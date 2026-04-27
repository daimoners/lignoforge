"""
Utility functions: graph ↔ SMILES/MOL conversion, visualisation,
random sampling helpers, and distance metrics.
"""

from __future__ import annotations

import os
import re
import platform
from collections import Counter
from typing import Optional, Tuple, TypeVar

import numpy as np
from numpy.random import RandomState
import networkx as nx
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import norm
from pysmiles import write_smiles, read_smiles, fill_valence
from rdkit import Chem
from rdkit.Chem import AllChem

from lignoforge.core.rules import (
    monomer_types,
    linkage_names,
    monomer_select_C1_C2,
    CHO,
    weight_CHO,
)

# Type aliases
nxgraph  = TypeVar("nxgraph")
molecule = TypeVar("molecule")
nparray  = TypeVar("nparray")

# Use non-interactive backend on Linux
if platform.system() == "Linux":
    matplotlib.use("Agg")


# ── Graph visualisation ───────────────────────────────────────────────────────

def draw_graph(
    G: nxgraph,
    node_labels: Optional[dict] = None,
    node_shape:  Optional[str]  = "o",
    node_size:   Optional[int]  = 500,
) -> None:
    """Draw a NetworkX graph with per-node colours."""
    plt.figure(figsize=(8, 8))
    nx.draw_networkx(
        G,
        with_labels=True,
        node_color=list(nx.get_node_attributes(G, "color").values()),
        labels=node_labels,
        node_shape=node_shape,
        node_size=node_size,
    )


def draw_big_graph(G: nxgraph) -> None:
    """Draw the coarse-grained monomer-level graph (hexagons)."""
    node_labels = {i: G.nodes[i]["mtype"] for i in range(len(G.nodes))}
    draw_graph(G, node_labels, node_shape="h", node_size=1000)


def draw_atomic_graph(G: nxgraph) -> None:
    """Draw the full atomistic graph (circles labelled by element)."""
    node_labels = {i: G.nodes[i]["element"] for i in range(len(G.nodes))}
    draw_graph(G, node_labels)


# ── SMILES / RDKit helpers ────────────────────────────────────────────────────

def graph_to_smile(G: nxgraph) -> str:
    """Convert an atomistic NetworkX graph to a SMILES string."""
    Gcopy = G.copy()
    fill_valence(Gcopy, respect_hcount=True)
    return write_smiles(Gcopy)


def graph_to_mol(
    G: nxgraph,
    save_mol:  Optional[bool] = False,
    name:      Optional[str]  = "lignoforge_mol",
    save_path: Optional[str]  = None,
    generate_3d: bool = False,
    random_seed: int = 42,
    optimize_3d: bool = True,
    max_uff_iterations: int = 500,
    add_hs_for_3d: bool = False,
    return_index_map: bool = False,
    chain_label: str = "",
) -> molecule:
    """
    Convert graph to an RDKit Mol object.

    When ``generate_3d`` is True, explicit Hs are added before ETKDGv3
    embedding so that (a) RDKit produces better conformers, (b) the
    "Molecule does not have explicit Hs" warnings are eliminated.
    Hs are stripped at the end if ``add_hs_for_3d`` is False.
    """
    if save_path is None:
        save_path = os.getcwd()
    mol, idx_to_node = _graph_to_ordered_rdkit_mol(G)

    # Fallback to SMILES path if direct construction fails
    if mol is None:
        smiles = graph_to_smile(G)
        mol = Chem.MolFromSmiles(smiles)
        idx_to_node = None

    if mol is None:
        return (None, None) if return_index_map else None

    if generate_3d:
        # ── helpers ──────────────────────────────────────────────────────────
        def _make_etkdg_params(seed, random_coords=False):
            params = AllChem.ETKDGv3()
            params.randomSeed = int(seed)
            # numThreads=1: process-level parallelism is already in use
            # (ProcessPoolExecutor); letting each worker grab all CPU threads
            # causes oversubscription and actually slows things down.
            params.numThreads = 1
            # useRandomCoords=True is more robust for macrocycles but much
            # slower.  Try standard DG first; fall back to random coords only
            # when the fast path fails.
            params.useRandomCoords = random_coords
            return params

        def _optimize_geometry(target_mol, max_iters):
            if not optimize_3d:
                return
            # Report progress at 0%, 25%, 50%, 75%, 100% of iterations.
            # Simple newline prints (no \r) so parallel workers don't overwrite
            # each other on the terminal.
            milestones = {max(1, int(max_iters * f)) for f in (0.25, 0.5, 0.75)}
            prefix = f"              {chain_label}" if chain_label else "             "

            try:
                mmff_props = AllChem.MMFFGetMoleculeProperties(target_mol)
            except Exception:
                mmff_props = None

            try:
                if mmff_props is not None:
                    ff = AllChem.MMFFGetMoleculeForceField(target_mol, mmff_props)
                    if ff is None:
                        raise RuntimeError("MMFF force field unavailable")
                    label = "MMFF"
                else:
                    ff = AllChem.UFFGetMoleculeForceField(target_mol)
                    if ff is None:
                        raise RuntimeError("UFF force field unavailable")
                    label = "UFF"

                total = int(max_iters)
                extra = max(100, total // 2)
                grand_total = total + extra
                done = 0
                converged = False
                chunk = 25   # iterations per ff.Minimize() call

                print(f"{prefix}• Optimizing geometry with {label} (max {grand_total} iter)")

                while done < grand_total and not converged:
                    step = min(chunk, grand_total - done)
                    result = ff.Minimize(maxIts=step)
                    done += step
                    if result == 0:
                        converged = True
                    elif done in milestones:
                        pct = int(done / total * 100)
                        print(f"{prefix}  {pct}% ({done}/{total} iter)...")

                status = "converged" if converged else "max iter reached"
                print(f"{prefix}  done ({done} iter, {status})")

            except Exception:
                # Fallback to single blocking call (e.g. force field init failed)
                try:
                    if mmff_props is not None:
                        print(f"{prefix}• Optimizing geometry with MMFF...")
                        AllChem.MMFFOptimizeMolecule(target_mol, maxIters=int(max_iters))
                    else:
                        print(f"{prefix}• Optimizing geometry with UFF...")
                        AllChem.UFFOptimizeMolecule(target_mol, maxIters=int(max_iters))
                except Exception:
                    pass

        # ── always embed with explicit Hs ─────────────────────────────────────
        prefix = f"              {chain_label}" if chain_label else "             "
        print(f"{prefix}• Embedding 3D coordinates (ETKDGv3)...")
        mol_with_hs = Chem.AddHs(mol)

        # Attempt 1: fast standard DG (no random coordinate initialisation)
        embed_status = AllChem.EmbedMolecule(mol_with_hs, _make_etkdg_params(random_seed, False))
        if embed_status != 0:
            # Attempt 2: robust random-coords init (handles complex ring closures)
            embed_status = AllChem.EmbedMolecule(
                mol_with_hs, _make_etkdg_params(random_seed + 1000, True)
            )

        if embed_status == 0:
            _optimize_geometry(mol_with_hs, max_uff_iterations)
            mol = mol_with_hs if add_hs_for_3d else Chem.RemoveHs(mol_with_hs)
        else:
            # Embedding failed completely – return mol without 3D coordinates
            mol = mol_with_hs if add_hs_for_3d else mol
    else:
        mol = Chem.RemoveHs(mol)

    if save_mol:
        os.makedirs(save_path, exist_ok=True)
        filename = os.path.join(save_path, name + ".png")
        Chem.Draw.MolToFile(mol, filename, size=(500, 500))
    if return_index_map:
        return mol, idx_to_node
    return mol


def _graph_to_ordered_rdkit_mol(G: nxgraph):
    """
    Build an RDKit mol from the atomistic graph preserving node index order.

    Returns
    -------
    mol : rdkit.Chem.Mol | None
    idx_to_node : dict[int, int] | None
        Mapping from RDKit atom indices (heavy atoms before AddHs) to
        original NetworkX node ids.
    """
    try:
        rw = Chem.RWMol()
        node_to_idx = {}
        idx_to_node = {}

        for node_id, attrs in sorted(G.nodes(data=True), key=lambda x: int(x[0])):
            atom = Chem.Atom(str(attrs.get("element", "C")))
            # Mark aromatic ring atoms so RDKit perceives the 6π system correctly.
            if attrs.get("aromatic", False):
                atom.SetIsAromatic(True)
            rd_idx = rw.AddAtom(atom)
            node_to_idx[int(node_id)] = int(rd_idx)
            idx_to_node[int(rd_idx)] = int(node_id)

        for n1, n2, attrs in G.edges(data=True):
            i1 = node_to_idx[int(n1)]
            i2 = node_to_idx[int(n2)]
            order = int(attrs.get("order", 1))

            n1_attrs = G.nodes[int(n1)]
            n2_attrs = G.nodes[int(n2)]
            n1_arom  = n1_attrs.get("aromatic", False)
            n2_arom  = n2_attrs.get("aromatic", False)
            # Intra-monomer ring bond: both atoms are part of the same aromatic
            # ring (same monomer index mi).  Inter-monomer bonds such as 5-5
            # biaryl also connect two aromatic carbons but from DIFFERENT
            # monomers (different mi) and must stay as single bonds.
            same_mi  = (n1_attrs.get("mi") is not None
                        and n1_attrs.get("mi") == n2_attrs.get("mi"))
            if n1_arom and n2_arom and same_mi:
                btype = Chem.BondType.AROMATIC
            elif order == 2:
                btype = Chem.BondType.DOUBLE
            else:
                btype = Chem.BondType.SINGLE

            if rw.GetBondBetweenAtoms(i1, i2) is None:
                rw.AddBond(i1, i2, btype)

        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        return mol, idx_to_node
    except Exception:
        return None, None


def smiles_to_formula(smiles: str) -> str:
    """Return the CₓHₓOₓ molecular formula for a SMILES string."""
    G_with_H = read_smiles(smiles, explicit_hydrogen=True)
    count_CHO = Counter(dict(G_with_H.nodes(data="element")).values())
    return "".join(f"{k}{count_CHO[k]}" for k in CHO)


def formula_to_MW(formula: str) -> float:
    """Compute molecular weight (g/mol) from a CₓHₓOₓ formula string."""
    counts = list(map(int, re.findall(r"\d+", formula)))
    return sum(weight_CHO[k] * counts[i] for i, k in enumerate(CHO))


# ── MW aggregation ────────────────────────────────────────────────────────────

def MW_array_to_number_average(MW: nparray) -> float:
    """Number-average molecular weight Mₙ."""
    return float(np.mean(MW))


def MW_array_to_weight_average(MW: nparray) -> float:
    """Weight-average molecular weight Mw."""
    MW = np.asarray(MW)
    return float(np.sum(MW ** 2) / np.sum(MW))


# ── Graph union / relabelling ─────────────────────────────────────────────────

def join_two(G1: nxgraph, G2: nxgraph) -> nxgraph:
    """Disjoint union of two graphs with contiguous node indices."""
    return nx.disjoint_union(G1, G2)


def adjust_indices(G: nxgraph) -> nxgraph:
    """Relabel nodes to a contiguous integer range after node removal."""
    mapping = {ni: ni_new for ni_new, ni in enumerate(G)}
    return nx.relabel_nodes(G, mapping)


# ── Bonding availability helpers ──────────────────────────────────────────────

def make_available(G: nxgraph, node_index: int) -> nxgraph:
    """Mark a carbon atom as available for inter-monomer bonding."""
    G.nodes[node_index]["bonding"] = True
    return G


def make_unavailable(G: nxgraph, node_index: int) -> nxgraph:
    """Mark a carbon atom as unavailable and reset neighbouring bond orders."""
    G.nodes[node_index]["bonding"] = False
    for u, v in G.edges(node_index):
        G.edges[u, v]["order"] = 1
    return G


def make_multi_available(G: nxgraph, monomer_type: str) -> nxgraph:
    """Mark all chemically valid bonding carbons as available in a fresh monomer."""
    if monomer_type not in monomer_types:
        raise ValueError(f"Unknown monomer type '{monomer_type}'. Must be H, G or S.")
    for Ci in monomer_select_C1_C2[monomer_type]:
        G = make_available(G, Ci - 1)  # convert 1-based chemical index to 0-based
    return G


# ── Random sampling ───────────────────────────────────────────────────────────

def select_one_from_many(many: list):
    """Uniformly sample one item from *many*."""
    return many[int(np.random.choice(len(many), 1)[0])]


def set_random_state(seed: Optional[int] = None) -> RandomState:
    """Create a NumPy RandomState, optionally seeded."""
    return np.random.RandomState(seed)


def generate_random_monomer(
    monomer_distribution: nparray,
    random_state: Optional[RandomState] = None,
) -> str:
    """Sample a monomer type from a discrete distribution [H, G, S]."""
    if random_state is None:
        random_state = np.random
    j = int(np.min(np.nonzero(random_state.rand() < np.cumsum(monomer_distribution))))
    return monomer_types[j]


def generate_random_linkage(
    linkage_distribution: nparray,
    random_state: Optional[RandomState] = None,
) -> str:
    """Sample a linkage name from a discrete distribution."""
    if random_state is None:
        random_state = np.random
    j = int(np.min(np.nonzero(random_state.rand() < np.cumsum(linkage_distribution))))
    return linkage_names[j]


def generate_random_branching_state(
    branching_propensity: float,
    random_state: Optional[RandomState] = None,
) -> bool:
    """Return True (allow branch) with probability *branching_propensity*."""
    if random_state is None:
        random_state = np.random
    dist = [branching_propensity, 1.0 - branching_propensity]
    j = int(np.min(np.nonzero(random_state.rand() < np.cumsum(dist))))
    return [True, False][j]


def generate_random_size_from_distribution(
    mean_size: float,
    max_size:  float,
    distribution_scaling: Optional[float] = 1.0,
    size_in_MW:           Optional[bool]   = False,
    random_state:         Optional[RandomState] = None,
) -> float:
    """Draw a polymer size from a truncated normal distribution."""
    if random_state is None:
        random_state = np.random
    min_size = 260.0 if size_in_MW else 2.0
    mu, sigma = mean_size, (max_size - min_size) / 1.96 / 2 * distribution_scaling
    dist = stats.truncnorm(
        (min_size - mu) / sigma, (max_size - mu) / sigma, loc=mu, scale=sigma
    )
    return dist.rvs(1, random_state=random_state)


# ── Metrics helpers ───────────────────────────────────────────────────────────

def cal_distance(
    target:  nparray,
    current: nparray,
    weights: Optional[nparray] = None,
) -> float:
    """Weighted Euclidean distance between two metrics arrays."""
    diff = np.array(target) - np.array(current)
    if weights is not None:
        diff = diff * np.array(weights)
    return float(np.sqrt(np.dot(diff, diff)))


def counts_to_metrics(
    counts:     nparray,
    additional: Optional[bool] = False,
) -> nparray:
    """Normalise raw counts array to a probability metrics array."""
    n_monomer  = len(monomer_types)
    n_linkage  = len(linkage_names)

    monomer_counts  = counts[:n_monomer]
    linkage_counts  = counts[n_monomer : n_monomer + n_linkage]

    monomer_sum = np.sum(monomer_counts)
    linkage_sum = np.sum(linkage_counts)

    m_dist = monomer_counts / monomer_sum if monomer_sum > 0 else np.zeros(n_monomer)
    l_dist = linkage_counts / linkage_sum if linkage_sum > 0 else np.zeros(n_linkage)

    metrics = np.concatenate([m_dist, l_dist])

    if additional and len(counts) > n_monomer + n_linkage:
        metrics = np.concatenate([metrics, counts[n_monomer + n_linkage:]])

    return metrics


def metrics_array_to_dict(
    metrics_array: nparray,
    metrics_names: list,
) -> dict:
    """Convert a flat metrics array to a labelled dictionary."""
    return {name: float(val) for name, val in zip(metrics_names, metrics_array)}


# ── Plotting helpers ──────────────────────────────────────────────────────────

def plot_single_distribution(
    ypred_all:   nparray,
    ypred:       Optional[float] = None,
    yobj:        Optional[float] = None,
    metric_name: Optional[str]   = "x",
    save_path:   Optional[str]   = None,
) -> None:
    """Plot a histogram + Gaussian fit for a single population metric."""
    if ypred_all.ndim > 1:
        ypred_all = ypred_all[:, 0]

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.hist(ypred_all, density=True, bins=10, alpha=0.5, color="steelblue", label="Population")
    if ypred is not None:
        ax.axvline(ypred, color="C0", linestyle="dotted", linewidth=2, label="Current")
    if yobj is not None:
        ax.axvline(yobj, color="C1", linestyle="--", linewidth=2, label="Target")
    mu, sigma = np.mean(ypred_all), np.std(ypred_all)
    x_norm = np.linspace(mu - 3 * sigma, mu + 3 * sigma, 100)
    ax.plot(x_norm, norm.pdf(x_norm, mu, sigma), color="r", label="Gaussian fit")
    xlabel = metric_name
    if metric_name in linkage_names or metric_name in monomer_types:
        xlabel += " (%)"
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Normalised Density")
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
    ax.set_title(rf"$\mu$={mu:.3f}  $\sigma$={sigma:.3f}")
    if save_path is None:
        save_path = os.getcwd()
    os.makedirs(save_path, exist_ok=True)
    fig.savefig(
        os.path.join(save_path, f"dist_{metric_name.replace(' ', '_')}.png"),
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_metrics(
    metrics_target:     dict,
    metrics_current:    dict,
    metrics_population: nparray,
    metrics_names:      list,
    save_path:          Optional[str] = None,
) -> None:
    """Plot distribution for each metric in the population."""
    for mi, name in enumerate(metrics_names):
        yobj  = metrics_target.get(name) if name != "branching_coeff" else None
        ypred = metrics_current.get(name)
        plot_single_distribution(
            metrics_population[:, mi], ypred, yobj,
            metric_name=name, save_path=save_path,
        )


def plot_distance_trajectory(
    distances:       list,
    simulation_name: Optional[str] = "x",
    distance_name:   Optional[str] = None,
    save_path:       Optional[str] = None,
) -> None:
    """Plot the MCMC distance (objective) trajectory."""
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot(np.arange(len(distances)), distances)
    ax.set_xlabel("Iterations")
    ylabel = f"{distance_name}_Distance" if distance_name else "Distance"
    ax.set_ylabel(ylabel)
    if save_path is None:
        save_path = os.getcwd()
    os.makedirs(save_path, exist_ok=True)
    fig.savefig(
        os.path.join(save_path, f"dist_{simulation_name.replace(' ', '_')}.png"),
        bbox_inches="tight",
    )
    plt.close(fig)
