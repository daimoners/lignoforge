"""Tests for virtual chemical derivatisation (OH masking) and its force-field typing."""

from __future__ import annotations

import collections
import shutil
import subprocess

import pytest
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from lignoforge.core.builder import normalise_derivatization, resolve_spec
from lignoforge.core.derivatization import derivatize_polymer, find_sites, parse_reagent
from lignoforge.core.utils import graph_to_mol
from lignoforge.forcefield import build_chain_topology, write_gromacs_system

from test_forcefield import _atomistic, _grow

GMX = shutil.which("gmx") is not None
BASE = dict(n=5, mono=[0.1, 0.5, 0.4],
            lk={"beta-O-4": 3, "4-O-5": 1, "beta-beta": 1, "5-5": 1}, seed=3)


@pytest.fixture(scope="module")
def base():
    return _grow(BASE["n"], BASE["mono"], BASE["lk"], BASE["seed"])


def formula(polymer):
    return CalcMolFormula(Chem.AddHs(graph_to_mol(polymer.G, generate_3d=False)))


def counts(f: str) -> collections.Counter:
    import re
    return collections.Counter({e: int(n or 1) for e, n in re.findall(r"([A-Z][a-z]?)(\d*)", f)})


def rule(site, reagent, degree=1.0):
    return {"site": site, "reagent": reagent, "degree": degree}


# ── input handling ────────────────────────────────────────────────────────────

def test_reagent_parsing():
    assert parse_reagent("methyl") == ("methyl", 0)
    assert parse_reagent("Acetyl") == ("acyl", 1)
    assert parse_reagent("acyl:3") == ("acyl", 3)
    for bad in ("acyl", "acyl:0", "acyl:99", "ethyl", ""):
        with pytest.raises(ValueError):
            parse_reagent(bad)


def test_rule_validation():
    assert normalise_derivatization([rule("gamma", "METHYL", 0.5)]) == [rule("gamma", "methyl", 0.5)]
    with pytest.raises(ValueError):
        normalise_derivatization([rule("ortho", "methyl")])
    with pytest.raises(ValueError):
        normalise_derivatization([rule("gamma", "methyl", 1.5)])
    with pytest.raises(ValueError):
        normalise_derivatization([rule("gamma", "methyl"), rule("gamma", "acetyl")])
    assert resolve_spec({"derivatization": [rule("phenolic", "acetyl")]})["derivatization"]


def test_sites_and_degree(base):
    for kind in ("phenolic", "benzylic", "gamma"):
        assert len(find_sites(base.G, kind)) > 0
    new, rep = derivatize_polymer(base, [rule("phenolic", "methyl", 0.5)], seed=1)
    r = rep["rules"][0]
    assert r["reacted_sites"] == round(0.5 * r["available_sites"])
    assert len(find_sites(new.G, "phenolic")) == r["available_sites"] - r["reacted_sites"]
    assert len(find_sites(base.G, "phenolic")) == r["available_sites"]      # input untouched


def test_reproducible_and_seed_dependent(base):
    a = derivatize_polymer(base, [rule("phenolic", "acetyl", 0.5)], seed=7)[1]["rules"][0]["nodes"]
    b = derivatize_polymer(base, [rule("phenolic", "acetyl", 0.5)], seed=7)[1]["rules"][0]["nodes"]
    assert a == b
    seen = {tuple(derivatize_polymer(base, [rule("phenolic", "acetyl", 0.5)], seed=s)[1]["rules"][0]["nodes"])
            for s in range(10)}
    assert len(seen) > 1


# ── chemistry ─────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("reagent,delta", [
    ("methyl", {"C": 1, "H": 2}),                 # O-H -> O-CH3
    ("acetyl", {"C": 2, "H": 2, "O": 1}),         # O-H -> O-C(=O)CH3
    ("acyl:3", {"C": 4, "H": 6, "O": 1}),         # butyryl
])
@pytest.mark.parametrize("site", ["phenolic", "benzylic", "gamma"])
def test_formula_change_per_site(base, site, reagent, delta):
    new, rep = derivatize_polymer(base, [rule(site, reagent)], seed=0)
    n = rep["rules"][0]["reacted_sites"]
    assert n > 0
    diff = counts(formula(new))
    diff.subtract(counts(formula(base)))
    assert {k: v for k, v in diff.items() if v} == {k: v * n for k, v in delta.items()}


def test_oh_hydrogens_are_consumed(base):
    def n_oh(poly):
        chain = _atomistic(poly)
        return sum(1 for m in chain["monomers"] for a in m["atoms"] if a["atom_name"] in ("HO4", "HOA", "HOG"))
    new, _ = derivatize_polymer(base, [rule(s, "acetyl") for s in ("phenolic", "benzylic", "gamma")], seed=0)
    assert n_oh(base) > 0 and n_oh(new) == 0


# ── force field ───────────────────────────────────────────────────────────────

ALL = ("phenolic", "benzylic", "gamma")


@pytest.mark.parametrize("reagent", ["methyl", "acetyl", "acyl:3"])
def test_topology_is_neutral_and_typed(base, reagent):
    new, rep = derivatize_polymer(base, [rule(s, reagent) for s in ALL], seed=0)
    topo = build_chain_topology(_atomistic(new))
    assert topo.net_charge == 0.0
    assert max(r["max_atom_shift"] for r in topo.charge_report) < 0.05   # local junction corrections only
    n_react = sum(r["reacted_sites"] for r in rep["rules"])
    types = collections.Counter(a.opls for a in topo.atoms)
    if reagent == "methyl":
        assert types["opls_181"] >= n_react
    else:
        assert types["opls_465"] == n_react and types["opls_466"] == n_react
        n_carbonyl_impropers = sum(1 for *_, m in topo.impropers if m == "improper_O_C_X_Y")
        assert n_carbonyl_impropers == n_react                          # one planarity term per ester


def test_aryl_ester_uses_phenyl_ester_types(base):
    new, _ = derivatize_polymer(base, [rule("phenolic", "acetyl")], seed=0)
    types = {a.name: a.opls for a in build_chain_topology(_atomistic(new)).atoms}
    assert types["C4"] == "opls_472" and types["O4H"] == "opls_473"


@pytest.mark.skipif(not GMX, reason="GROMACS not installed")
@pytest.mark.parametrize("rules", [
    [rule("phenolic", "methyl")],
    [rule(s, "acetyl") for s in ALL],
    [rule(s, "acyl:3") for s in ALL],
    [rule("phenolic", "acetyl", 0.5), rule("gamma", "methyl", 0.5)],
])
def test_grompp_accepts_derivatised_topologies(tmp_path, base, rules):
    new, _ = derivatize_polymer(base, rules, seed=2)
    write_gromacs_system(_atomistic(new), str(tmp_path), name="d")
    r = subprocess.run(["gmx", "grompp", "-f", "em.mdp", "-c", "d.gro", "-p", "d.top", "-o", "em.tpr"],
                       cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-1200:]


# ── integration with the project API ──────────────────────────────────────────

def test_project_build_with_derivatization(tmp_path):
    from lignoforge.api import Project
    p = Project.create(tmp_path / "p", "p")
    plain = p.build_chain({"n_monomers": 4, "monomers": {"G": 1}, "linkages": {"beta-O-4": 1}, "seed": 1})
    der = p.build_chain({"n_monomers": 4, "monomers": {"G": 1}, "linkages": {"beta-O-4": 1}, "seed": 1,
                         "derivatization": [rule("phenolic", "acetyl")]})
    assert der["status"] == "ready", der["error"]
    assert der["resolved"]["derivatization"][0]["reagent"] == "acetyl"
    rep = der["derivatization_report"]["rules"][0]
    assert rep["reacted_sites"] == rep["available_sites"] > 0
    assert der["stats"]["n_atoms"] == plain["stats"]["n_atoms"] + 5 * rep["reacted_sites"]   # C2H2O per site
    assert p.topology(der["id"])["net_charge"] == 0.0
