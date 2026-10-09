"""Tests for the MD workflow generation and the coarse-grained model."""

from __future__ import annotations

import json
import shutil
import subprocess

import numpy as np
import pytest

from lignoforge.cg import (
    BoltzmannFitter, CGParameters, default_parameters, map_topology,
    read_pdb_frames, write_cg_system,
)
from lignoforge.cg.boltzmann import KB, _harmonic_from_histogram
from lignoforge.forcefield import build_chain_topology
from lignoforge.md import mdp, write_atomistic_workflow

from test_forcefield import _atomistic, _grow

_LK = {"beta-O-4": 3, "4-O-5": 1, "beta-beta": 1, "5-5": 1}


@pytest.fixture(scope="module")
def chain():
    return _atomistic(_grow(7, [0.1, 0.5, 0.4], _LK, 4), seed=4)


# ── mdp ───────────────────────────────────────────────────────────────────────

def test_mdp_step_counts():
    assert "nsteps              = 50000" in mdp.nvt_mdp(100.0, 300.0)
    assert "ref-t               = 350.0" in mdp.prod_mdp(1000.0, 350.0)
    assert "pcoupl              = no" in mdp.nvt_mdp()
    assert "C-rescale" in mdp.npt_mdp()


# ── Atomistic workflow ────────────────────────────────────────────────────────

def test_workflow_solvated(tmp_path, chain):
    paths = write_atomistic_workflow(chain, str(tmp_path), name="lig")
    assert {"top", "gro", "em", "nvt", "npt", "prod", "run", "config"} <= set(paths)
    assert '#include "oplsaa.ff/tip3p.itp"' in open(paths["top"]).read()
    run = open(paths["run"]).read()
    assert "solvate" in run and "npt.mdp" in run
    assert json.load(open(paths["config"]))["net_charge"] == 0.0


def test_workflow_vacuum_has_no_npt_or_water(tmp_path, chain):
    paths = write_atomistic_workflow(chain, str(tmp_path), name="lig", solvent="none")
    assert "npt" not in paths
    assert "tip3p" not in open(paths["top"]).read()
    assert "solvate" not in open(paths["run"]).read()


def test_workflow_rejects_unknown_solvent(tmp_path, chain):
    with pytest.raises(ValueError):
        write_atomistic_workflow(chain, str(tmp_path), solvent="benzene")


# ── CG mapping ────────────────────────────────────────────────────────────────

def test_cg_mapping_matches_atomistic_graph(chain):
    topo = build_chain_topology(chain)
    cg = map_topology(topo)
    n_res = len({a.resid for a in topo.atoms})
    assert len(cg.beads) == n_res
    assert sum(len(b.atom_indices) for b in cg.beads) == len(topo.atoms)
    pairs = {tuple(sorted((l["residue_1"], l["residue_2"]))) for l in topo.linkages}
    assert {(i, j) for i, j, _ in cg.bonds} == pairs
    deg = {}
    for i, j, _ in cg.bonds:
        deg[i] = deg.get(i, 0) + 1
        deg[j] = deg.get(j, 0) + 1
    assert len(cg.angles) == sum(d * (d - 1) // 2 for d in deg.values())
    assert cg.ref_xyz.shape == (n_res, 3)


# ── CG parameters ─────────────────────────────────────────────────────────────

def test_parameter_lookup_and_roundtrip(tmp_path):
    p = default_parameters()
    assert p.bond("beta-O-4", "G", "S")["r0"] == pytest.approx(0.656)
    p.bonds["beta-O-4:G-S"] = {"r0": 0.7, "k": 1.0}
    assert p.bond("beta-O-4", "S", "G")["r0"] == 0.7           # type-specific wins
    assert p.angle("made-up", "also-made-up") == p.default_angle  # fallback
    path = str(tmp_path / "p.json")
    p.save(path)
    assert CGParameters.load(path).to_dict() == p.to_dict()


# ── Boltzmann inversion ───────────────────────────────────────────────────────

def test_boltzmann_recovers_known_harmonic_potentials():
    rng, kT = np.random.RandomState(0), KB * 300.0
    r0, k = 0.62, 3000.0
    r = np.linspace(0.3, 1.0, 4000)
    P = r ** 2 * np.exp(-k * (r - r0) ** 2 / 2 / kT)
    s = rng.choice(r, 150000, p=P / P.sum())
    got_r0, got_k = _harmonic_from_histogram(s, lambda c: c ** 2, kT, s.min(), s.max())
    assert got_r0 == pytest.approx(r0, abs=0.003)
    assert got_k == pytest.approx(k, rel=0.05)

    t0, ka = np.radians(140.0), 30.0
    t = np.linspace(0.3, np.pi, 4000)
    P = np.sin(t) * np.exp(-ka * (t - t0) ** 2 / 2 / kT)
    s = rng.choice(t, 150000, p=P / P.sum())
    got_t0, got_k = _harmonic_from_histogram(s, np.sin, kT, 0.2, np.pi)
    assert np.degrees(got_t0) == pytest.approx(140.0, abs=0.5)
    assert got_k == pytest.approx(ka, rel=0.05)


def test_fitter_keeps_defaults_when_undersampled(chain):
    topo = build_chain_topology(chain)
    cg = map_topology(topo)
    xyz = np.array([a.xyz for a in topo.atoms])[None] / 10.0
    fitter = BoltzmannFitter()
    fitter.add(cg, xyz)                       # a single frame
    params, report = fitter.fit(min_samples=200)
    assert report["fitted"] == {}
    assert params.to_dict()["bonds"] == default_parameters().to_dict()["bonds"]


def test_read_pdb_frames(tmp_path):
    lines = []
    for m in range(3):
        lines.append(f"MODEL {m + 1}")
        for a in range(4):
            lines.append(f"ATOM  {a + 1:5d}  C   LIG A   1    "
                         f"{10 * m + a:8.3f}{0:8.3f}{0:8.3f}")
        lines.append("ENDMDL")
    f = tmp_path / "f.pdb"
    f.write_text("\n".join(lines) + "\n")
    fr = read_pdb_frames(str(f), 4)
    assert fr.shape == (3, 4, 3)
    assert fr[2, 3, 0] == pytest.approx(2.3)   # (20+3) Å → nm


# ── CG system ─────────────────────────────────────────────────────────────────

def test_cg_system_files(tmp_path, chain):
    paths = write_cg_system([chain], str(tmp_path), name="cg", copies=3)
    top = open(paths["top"]).read()
    assert "[ atomtypes ]" in top and "cg_0  3" in top
    assert "insert-molecules" in open(paths["run"]).read()
    assert json.load(open(tmp_path / "workflow.json"))["copies"] == 3


@pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS not installed")
def test_cg_system_grompp(tmp_path, chain):
    write_cg_system([chain], str(tmp_path), name="cg", run_ns=0.1)
    r = subprocess.run(
        ["gmx", "grompp", "-f", "md.mdp", "-c", "cg_0.gro", "-p", "cg.top",
         "-o", "md.tpr"], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-1500:]


def test_angle_fit_is_constrained_to_physical_domain():
    """A distribution truncated at 180° must not give theta0 > 180°."""
    rng, kT = np.random.RandomState(1), KB * 300.0
    t = np.linspace(0.3, np.pi, 4000)
    P = np.sin(t) * np.exp(-25.0 * (t - np.radians(185.0)) ** 2 / 2 / kT)
    s = rng.choice(t, 200000, p=P / P.sum())
    t0, k = _harmonic_from_histogram(s, np.sin, kT, 0.2, np.pi, bounds=(0.0, np.pi))
    assert np.degrees(t0) == pytest.approx(180.0)
    assert k > 0


@pytest.mark.parametrize("out_ps", [0.1, 0.5, 1.0, 7.0, 10.0])
def test_output_intervals_are_multiples_of_nstcalcenergy(out_ps):
    """GROMACS refuses nstenergy that is not a multiple of nstcalcenergy."""
    txt = mdp.prod_mdp(100.0, 300.0, out_ps=out_ps)
    val = {l.split("=")[0].strip(): int(l.split("=")[1])
           for l in txt.splitlines()
           if l.split("=")[0].strip() in ("nstcalcenergy", "nstenergy", "nstlog")}
    assert val["nstenergy"] % val["nstcalcenergy"] == 0
    assert val["nstlog"] % val["nstcalcenergy"] == 0
