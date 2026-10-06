"""
Tests for the graph-based OPLS-AA topology generation
(:mod:`lignoforge.forcefield`) and for linear chain growth.

Run with ``python -m pytest tests/ -v``.  The grompp integration test is
skipped automatically when GROMACS is not installed.
"""

from __future__ import annotations

import collections
import copy
import io
import contextlib
import shutil
import subprocess

import numpy as np
import pytest

from lignoforge.cli.build_chain import _grow_chain_exact
from lignoforge.core.monomer import Monomer
from lignoforge.core.polymer import Polymer
from lignoforge.forcefield import TypingError, build_chain_topology, write_gromacs_system
from lignoforge.forcefield.charges import neutralize_residue
from lignoforge.forcefield.opls import linkage_position
from lignoforge.structure.generator import MolecularStructureGenerator

LINKAGES = ["beta-O-4", "alpha-O-4", "4-O-5", "5-5", "beta-5", "beta-beta"]
# order of lignoforge.core.rules.linkage_names
_LK_ORDER = ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta", "beta-1"]

_GEN = MolecularStructureGenerator()


def _atomistic(polymer, seed: int = 1) -> dict:
    with contextlib.redirect_stdout(io.StringIO()):
        return _GEN.polymer_atomistic_topology(
            polymer, random_seed=seed, optimize_3d=False
        )


def _grow(n, mono, linkages: dict, seed, branching=0.0):
    lk = np.array([linkages.get(k, 0.0) for k in _LK_ORDER], dtype=float)
    lk /= lk.sum()
    m = np.array(mono, dtype=float)
    m /= m.sum()
    with contextlib.redirect_stdout(io.StringIO()):
        return _grow_chain_exact(n, m.tolist(), lk.tolist(), seed, branching)


# ── Charge neutralisation ─────────────────────────────────────────────────────

@pytest.mark.parametrize("residual", [0.0, 0.2, -0.24, 0.085, 0.00031])
def test_neutralize_residue_sums_to_exactly_zero(residual):
    q = {f"A{i}": (0.1 if i % 2 else -0.1) for i in range(32)}
    q["A0"] += residual
    new, res, shift = neutralize_residue(q)
    assert res == pytest.approx(residual)
    assert sum(round(v * 10000) for v in new.values()) == 0
    assert shift <= abs(residual) / len(q) + 2e-4


# ── Linkage validation ────────────────────────────────────────────────────────

def test_linkage_position_labels_ring_closure():
    assert linkage_position("beta-O-4", "CB", "O4H") == "CB"
    assert linkage_position("beta-5", "CA", "O4H") == "CA_ringclose"
    assert linkage_position("beta-5", "C5", "CB") == "C5"
    assert linkage_position("beta-beta", "OG", "CA") == "OG_ringclose"


def test_linkage_position_rejects_invalid_pairs():
    with pytest.raises(TypingError):
        linkage_position("beta-O-4", "CB", "C5")
    with pytest.raises(TypingError):
        linkage_position("made-up", "CB", "O4H")


# ── Linear growth ─────────────────────────────────────────────────────────────

@pytest.mark.parametrize("seed", range(6))
def test_branching_zero_gives_linear_chain(seed):
    p = _grow(8, [0.2, 0.4, 0.4],
              {"4-O-5": 1, "beta-O-4": 1, "5-5": 1, "beta-5": 1}, seed)
    assert max(dict(p.bigG.degree()).values()) <= 2
    assert nx_is_tree(p.bigG)


def nx_is_tree(g) -> bool:
    import networkx as nx
    return nx.is_tree(g)


def test_branching_positive_can_branch():
    degs = []
    for seed in range(12):
        p = _grow(10, [0.3, 0.4, 0.3], {"4-O-5": 1, "beta-O-4": 1, "5-5": 1},
                  seed, branching=0.6)
        degs.append(max(dict(p.bigG.degree()).values()))
    assert max(degs) >= 3


# ── Topology construction ─────────────────────────────────────────────────────

@pytest.mark.parametrize("linkage", LINKAGES)
@pytest.mark.parametrize("mtype", ["G", "H"])
def test_dimer_topology_is_complete_and_neutral(linkage, mtype):
    p = Polymer(Monomer("G"), verbose=False)
    if not p.add_specific_monomer(mtype, linkage):
        pytest.skip(f"{linkage} G-{mtype} not constructible")
    chain = _atomistic(p)
    topo = build_chain_topology(chain)

    assert topo.net_charge == 0.0
    for res in topo.charge_report:
        assert res["max_atom_shift"] < 0.02
    assert {lk["linkage"] for lk in topo.linkages} <= {linkage}
    # every graph bond is present exactly once
    assert len(topo.bonds) == len(set(topo.bonds)) == len(chain["bonds"])
    # all atoms typed
    assert all(a.opls.startswith("opls_") for a in topo.atoms)


def test_aromatic_impropers_one_per_ring_carbon():
    p = _grow(4, [0, 1, 1], {"beta-O-4": 1, "5-5": 1}, 3)
    topo = build_chain_topology(_atomistic(p))
    assert len(topo.impropers) == 6 * 4


def test_topology_independent_of_residue_numbering():
    """Re-numbering the monomers must not change the chemistry."""
    p = _grow(7, [0.2, 0.4, 0.4],
              {"4-O-5": 1, "beta-O-4": 1, "beta-beta": 1}, 5)
    chain = _atomistic(p)
    ids = sorted({m["monomer_id"] for m in chain["monomers"]})
    perm = dict(zip(ids, reversed(ids)))

    shuffled = copy.deepcopy(chain)
    for m in shuffled["monomers"]:
        m["monomer_id"] = perm[m["monomer_id"]]
        for a in m["atoms"]:
            a["monomer_id"] = m["monomer_id"]
    for b in shuffled["bonds"]:
        b["monomer1_id"] = perm[b["monomer1_id"]]
        b["monomer2_id"] = perm[b["monomer2_id"]]

    t1, t2 = build_chain_topology(chain), build_chain_topology(shuffled)
    sig = lambda t: (  # noqa: E731
        collections.Counter(a.opls for a in t.atoms),
        len(t.bonds), len(t.angles), len(t.dihedrals), len(t.pairs),
        len(t.impropers), t.net_charge,
    )
    assert sig(t1) == sig(t2)


def test_every_inter_monomer_bond_has_a_linkage():
    p = _grow(6, [0.3, 0.4, 0.3], {"beta-O-4": 2, "4-O-5": 1, "beta-5": 1}, 2)
    chain = _atomistic(p)
    n_inter = sum(1 for b in chain["bonds"] if b["scope"] == "inter_monomer")
    assert len(build_chain_topology(chain).linkages) == n_inter


def test_acetal_alpha_carbon_is_rejected():
    """A Cα bearing both α-OH and an ether bond must raise, not be guessed."""
    from lignoforge.forcefield.opls import type_residue

    names = {"C1", "C2", "H2", "C3", "H3", "C4", "O4H", "C5", "H5", "C6", "H6",
             "CA", "HA", "OA", "HOA", "CB", "HB", "CG", "HG1", "HG2", "OG", "HOG"}
    with pytest.raises(TypingError):
        type_residue(names, {"CA": ("alpha-O-4", "CA"), "CB": ("beta-O-4", "CB")})


# ── Writers ───────────────────────────────────────────────────────────────────

def test_writers_produce_consistent_files(tmp_path):
    p = _grow(4, [0, 1, 1], {"beta-O-4": 1}, 0)
    chain = _atomistic(p)
    paths = write_gromacs_system(chain, str(tmp_path), name="lig")
    top = open(paths["top"]).read()
    gro = open(paths["gro"]).read().splitlines()
    assert '#include "oplsaa.ff/forcefield.itp"' in top
    assert int(gro[1]) == chain["n_atoms"]
    assert len(gro) == chain["n_atoms"] + 3


@pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS not installed")
@pytest.mark.parametrize("seed,linkages,branching", [
    (0, {"beta-O-4": 1}, 0.0),
    (1, {"beta-O-4": 1, "4-O-5": 1, "beta-beta": 1}, 0.0),
    (2, {"beta-O-4": 1, "alpha-O-4": 1, "beta-5": 1, "5-5": 1}, 0.0),
    (3, {"beta-O-4": 1, "4-O-5": 1, "5-5": 1}, 0.5),
])
def test_grompp_accepts_generated_topology(tmp_path, seed, linkages, branching):
    p = _grow(6, [0.2, 0.4, 0.4], linkages, seed, branching)
    write_gromacs_system(_atomistic(p), str(tmp_path), name="lig")
    r = subprocess.run(
        ["gmx", "grompp", "-f", "em.mdp", "-c", "lig.gro", "-p", "lig.top",
         "-o", "em.tpr"],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr[-1500:]
