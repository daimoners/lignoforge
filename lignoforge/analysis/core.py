"""
Trajectory access and monomer bookkeeping shared by all analyses.

A *run layout* describes which atoms of the simulated system belong to the
lignin molecule(s)::

    [{"chain_id": "c0001", "n_atoms": 117, "copies": 1}, ...]

(atomistic runs have a single entry; CG runs one entry per chain, each with
its number of copies, in the order of the topology).  Solvent follows.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np

__all__ = ["AnalysisError", "Selection", "load_universe", "locate_molecule",
           "select_molecule", "monomer_groups", "graph_endpoints", "AMU_A3_TO_G_CM3"]

AMU_A3_TO_G_CM3 = 1.66053907


class AnalysisError(RuntimeError):
    """Raised for missing files, missing optional dependencies or bad input."""


@dataclass
class Selection:
    """One lignin molecule inside a Universe."""
    atoms: "object"                    # MDAnalysis AtomGroup
    monomers: List[np.ndarray]         # per-monomer index arrays (into `atoms`)
    chain_id: str


def _mda():
    try:
        import MDAnalysis as mda
    except ImportError as e:             # pragma: no cover
        raise AnalysisError(
            "MDAnalysis is required for trajectory analysis: "
            "pip install 'lignoforge[analysis]'") from e
    return mda


def load_universe(tpr: str, xtc: Optional[str] = None):
    """Load a GROMACS run (``.tpr`` + optional ``.xtc``) as a MDAnalysis Universe."""
    mda = _mda()
    for f in (tpr, xtc):
        if f and not os.path.exists(f):
            raise AnalysisError(f"File not found: {f}")
    return mda.Universe(tpr, xtc) if xtc else mda.Universe(tpr)


def monomer_groups(chain_data: dict, coarse_grained: bool) -> List[np.ndarray]:
    """
    Indices (relative to the molecule's first atom) of the atoms of each
    monomer, in residue order.  For CG runs every atom is one monomer bead.
    """
    from lignoforge.forcefield import build_chain_topology

    topo = build_chain_topology(chain_data)
    by_res: dict = {}
    for i, a in enumerate(topo.atoms):
        by_res.setdefault(a.resid, []).append(i)
    if coarse_grained:
        return [np.array([r - 1]) for r in sorted(by_res)]
    return [np.array(by_res[r]) for r in sorted(by_res)]


def graph_endpoints(chain_data: dict) -> tuple:
    """
    Pair of monomers (1-based residue numbers) with the largest graph
    distance, i.e. the 'ends' of the chain (for branched chains, of its
    longest path).
    """
    import networkx as nx
    from lignoforge.forcefield import build_chain_topology

    topo = build_chain_topology(chain_data)
    g = nx.Graph()
    g.add_nodes_from({a.resid for a in topo.atoms})
    for lk in topo.linkages:
        g.add_edge(lk["residue_1"], lk["residue_2"])
    best, pair = -1, (1, 1)
    for u, dists in nx.all_pairs_shortest_path_length(g):
        for v, d in dists.items():
            if d > best:
                best, pair = d, (u, v)
    return pair


def locate_molecule(layout: Sequence[dict], molecule: int = 0) -> tuple:
    """Return ``(layout_entry, first_atom_index)`` of the *molecule*-th molecule."""
    start = 0
    count = 0
    for entry in layout:
        for _ in range(int(entry.get("copies", 1))):
            if count == molecule:
                return entry, start
            start += int(entry["n_atoms"])
            count += 1
    raise AnalysisError(f"Molecule index {molecule} out of range ({count} molecules)")


def select_molecule(universe, chain_data: dict, entry: dict, start: int,
                    coarse_grained: bool = False, unwrap: bool = True) -> Selection:
    """
    Return the molecule that starts at atom *start* and is described by
    *entry* / *chain_data*, made whole across periodic boundaries.
    """
    ag = universe.atoms[start:start + int(entry["n_atoms"])]
    if unwrap and ag.n_atoms and len(ag.bonds):
        from MDAnalysis import transformations as trans
        universe.trajectory.add_transformations(trans.unwrap(ag))
    return Selection(ag, monomer_groups(chain_data, coarse_grained),
                     entry["chain_id"])
