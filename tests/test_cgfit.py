"""Tests of the CG parametriser (:mod:`lignoforge.cgfit`)."""

from __future__ import annotations

import contextlib
import io
import shutil

import numpy as np
import pytest

from lignoforge.cgfit import CGParams, CGParamJob, map_chain, quick_config
from lignoforge.cgfit import fit as F
from lignoforge.cgfit import gmxio as G
from lignoforge.cgfit.sampling import SystemMap, ideal_pair_numbers, sample
from lignoforge.core.builder import build_chain_dict
from lignoforge.forcefield import build_chain_topology

needs_gmx = pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS not installed")

SPEC = {"n_monomers": 5, "monomers": {"H": 0.2, "G": 0.4, "S": 0.4},
        "linkages": {"beta-O-4": 0.5, "5-5": 0.2, "beta-5": 0.2, "4-O-5": 0.1}, "seed": 3}


@pytest.fixture(scope="module")
def topo():
    with contextlib.redirect_stdout(io.StringIO()):
        chain, _ = build_chain_dict(SPEC, optimize=False)
    return build_chain_topology(chain, name="L0")


# ── mapping ───────────────────────────────────────────────────────────────────

def test_ring_side_mapping_partitions_atoms_and_mass(topo):
    m = map_chain(topo, "ring-side")
    assert len(m.beads) == 2 * 5
    assert sorted(a for b in m.beads for a in b.atoms) == list(range(len(topo.atoms)))
    assert sum(m.masses) == pytest.approx(sum(m.atom_mass))
    side = [b for b in m.beads if b.type == "SC"]
    names = {topo.atoms[a].name for b in side for a in b.atoms}
    assert {"CA", "CB", "CG", "OG"} <= names and "C1" not in names and "HO4" not in names
    assert {t for t in m.types} <= {"RH", "RG", "RS", "SC"}


def test_monomer_mapping_one_bead_per_residue(topo):
    m = map_chain(topo, "monomer")
    assert len(m.beads) == 5 and all(b.type.startswith("M") for b in m.beads)


def test_bead_graph_is_connected_tree_for_linear_chain(topo):
    import networkx as nx
    m = map_chain(topo, "ring-side")
    g = nx.Graph([(i, j) for i, j, _ in m.bonds])
    assert g.number_of_nodes() == len(m.beads) and nx.is_connected(g)
    assert g.number_of_edges() == len(m.beads) - 1
    inter = [k for _, _, k in m.bonds if not k.startswith("intra")]
    assert len(inter) == 4
    # exclusions: every 1-2 and 1-3 pair, and nothing else
    pairs = set(m.exclusions)
    assert all((i, j) in pairs for i, j, _ in m.bonds)
    assert len(pairs) == len(m.bonds) + len(m.angles)


# ── parameters and topology writer ────────────────────────────────────────────

def _params(m, sigma=0.5, eps=2.0):
    pairs = sorted({"-".join(sorted((a, b))) for a in m.types for b in m.types})
    return CGParams("ring-side", 500.0, 1.0, {p: {"sigma": sigma, "epsilon": eps} for p in pairs})


def test_params_roundtrip(tmp_path, topo):
    p = _params(map_chain(topo))
    p.save(str(tmp_path / "p.json"))
    assert CGParams.load(str(tmp_path / "p.json")).nonbonded == p.nonbonded


def test_cg_topology_has_pairs_exclusions_and_valid_c6c12(topo):
    m = map_chain(topo)
    text = G.cg_top_text([m], [3], _params(m))
    assert "[ nonbond_params ]" in text and "[ exclusions ]" in text
    assert text.count(f"{m.name}  3") == 1
    c6, c12 = G._lj_c6c12(0.5, 2.0)
    assert c6 == pytest.approx(4 * 2.0 * 0.5 ** 6) and c12 == pytest.approx(4 * 2.0 * 0.5 ** 12)


@needs_gmx
def test_cg_topology_accepted_by_grompp(tmp_path, topo):
    m = map_chain(topo)
    p = _params(m)
    G.write_text(tmp_path / "cg.top", G.cg_top_text([m], [1], p))
    x = np.array([[0.4 * i, 0.1 * (i % 2), 0.0] for i in range(len(m.beads))]) + 2.0
    G.write_text(tmp_path / "s.gro", G.cg_gro_text(m.types, x, [4.0] * 3))
    G.write_text(tmp_path / "em.mdp", G.cg_mdp("em", 0, 500.0, 1.0))
    G.Gmx().grompp("em.mdp", "s.gro", "cg.top", "em.tpr", tmp_path)


# ── sampling ──────────────────────────────────────────────────────────────────

def test_sampling_pair_histogram_excludes_bonded_neighbours(topo):
    m = map_chain(topo)
    sm = SystemMap([m], [2])
    rng = np.random.default_rng(1)
    box = np.array([3.0, 3.0, 3.0])
    frames = [(rng.uniform(0, 3, (sm.n_beads, 3)), box) for _ in range(4)]
    st = sample(sm, iter(frames), mapped=True, rc=1.0)
    n_pairs = sm.n_beads * (sm.n_beads - 1) // 2
    n_excl = len(sm.excl_codes)
    # all counted pairs are non-excluded, and a random gas gives g(r) ~ 1
    assert 0 < st.pair_hist[0].sum() < n_pairs
    ideal = ideal_pair_numbers(sm)
    assert sum(ideal.values()) == pytest.approx(n_pairs)
    assert n_excl > 0 and st.n_frames == 4


# ── fitting ───────────────────────────────────────────────────────────────────

def test_lj_derivatives_match_finite_differences():
    r = np.linspace(0.3, 1.2, 50)
    s, e, rc, h = 0.5, 2.0, 1.0, 1e-6
    u, de, ds = F.lj(r, s, e, rc)
    ue = (F.lj(r, s, e + h, rc)[0] - F.lj(r, s, e - h, rc)[0]) / (2 * h)
    us = (F.lj(r, s + h, e, rc)[0] - F.lj(r, s - h, e, rc)[0]) / (2 * h)
    np.testing.assert_allclose(de, ue, atol=1e-5)
    np.testing.assert_allclose(ds, us, rtol=1e-4, atol=1e-4)
    assert np.all(u[r >= rc] == 0)


def test_mixing_lorentz_berthelot_and_jacobian():
    mix = F.Mixing(["A-A", "A-B", "B-B"], "lb")
    theta = np.array([2.0, 0.5, 3.0, 0.7])           # eps_A, sig_A, eps_B, sig_B
    d = mix.pairs_from(theta)
    assert d["A-B"]["epsilon"] == pytest.approx(np.sqrt(6.0))
    assert d["A-B"]["sigma"] == pytest.approx(0.6)
    J = mix.jacobian(theta)
    assert J.shape == (6, 4) and J[2, 0] == pytest.approx(0.5 * np.sqrt(3.0 / 2.0), rel=1e-3)


def test_pressure_step_direction():
    nb = {"A-A": {"sigma": 0.5, "epsilon": 2.0}}
    up, f_up = F.pressure_step(nb, p_cg=900.0, p_ref=1.0)
    down, f_down = F.pressure_step(nb, p_cg=-900.0, p_ref=1.0)
    assert f_up > 1 > f_down and up["A-A"]["epsilon"] > 2.0 > down["A-A"]["epsilon"]


def _sample_boltzmann(x, U, jac, n, rng):
    """Draw from P(x) ~ jac(x) exp(-U/kT) by inverse-CDF on a grid."""
    w = jac(x) * np.exp(-U)
    cdf = np.cumsum(w) / w.sum()
    return np.interp(rng.random(n), cdf, x)


def test_bonded_inversion_recovers_harmonic_parameters():
    """Samples drawn from the Boltzmann law of a known potential give that potential back."""
    from lignoforge.cgfit.sampling import Stats
    kT = 0.0083144626 * 300.0
    rng = np.random.default_rng(0)
    k, r0, ka, t0 = 4000.0, 0.45, 60.0, np.radians(125.0)
    r = np.linspace(0.2, 0.8, 4000)
    th = np.linspace(0.05, np.pi, 4000)
    st = Stats(["A-A"], np.linspace(0, 1, 3))
    st.bonds["intra:RG-SC"] = [_sample_boltzmann(r, 0.5 * k * (r - r0) ** 2 / kT, lambda x: x ** 2, 20000, rng)]
    st.angles["RG:intra:RG-SC|intra:RG-SC"] = [np.degrees(_sample_boltzmann(
        th, 0.5 * ka * (th - t0) ** 2 / kT, np.sin, 20000, rng))]
    bonds, angles, rep = F.fit_bonded(st, 300.0, min_samples=100)
    assert bonds["intra:RG-SC"]["r0"] == pytest.approx(r0, abs=0.003)
    assert bonds["intra:RG-SC"]["k"] == pytest.approx(k, rel=0.15)
    assert bonds["RG-SC"]["r0"] == pytest.approx(r0, abs=0.003)       # generic fallback
    a = list(angles.values())[0]
    assert a["theta0"] == pytest.approx(125.0, abs=1.5) and a["k"] == pytest.approx(ka, rel=0.2)


# ── monitor ───────────────────────────────────────────────────────────────────

def test_monitor_reports_progress_and_eta_from_gromacs_logs(tmp_path):
    import json
    import os
    import time
    from datetime import datetime, timedelta

    from lignoforge.cgfit.monitor import current_mdrun, fmt_hms, overview

    d = tmp_path / "reference" / "prod"
    d.mkdir(parents=True)
    (d / "prod.tpr").write_bytes(b"x")
    start = datetime.now() - timedelta(seconds=100)
    (d / "prod.log").write_text(
        "   dt                             = 0.002\n   nsteps                         = 10000\n"
        f"Started mdrun on rank 0 {start.strftime('%a %b %d %H:%M:%S %Y')}\n"
        "           Step           Time\n        2500        5.00000\n")
    cur = current_mdrun(tmp_path)
    assert cur["run"] == "reference/prod/prod.log" and cur["fraction"] == pytest.approx(0.25)
    assert 250 < cur["eta_s"] < 350                      # 25% in ~100 s -> ~300 s left
    assert cur["ns_per_day"] == pytest.approx(2500 / 100 * 0.002 * 86.4, rel=0.1)
    (tmp_path / "config.json").write_text(json.dumps(quick_config()))
    ov = overview(tmp_path)
    assert ov["current"]["step"] == 2500 and ov["stages"]["reference"] == "pending"
    # a finished or stale log is not "running"
    (d / "prod.log").write_text((d / "prod.log").read_text() + "Finished mdrun on rank 0\n")
    assert current_mdrun(tmp_path) is None
    assert fmt_hms(3725) == "1h 02m" and fmt_hms(75) == "1m 15s"


# ── end-to-end (needs GROMACS) ───────────────────────────────────────────────

@needs_gmx
def test_quick_pipeline_runs_all_stages(tmp_path):
    cfg = quick_config()
    cfg["atomistic"].update(compress_ps=10.0, equil_ps=10.0, prod_ps=20.0, out_ps=2.0)
    cfg["cg"].update(iterations=2, iter_ps=40.0, validate_ps=40.0, out_ps=1.0)
    job = CGParamJob(tmp_path, cfg)
    with contextlib.redirect_stdout(io.StringIO()):
        state = job.run()
    assert all(state[s]["status"] == "done" for s in
               ["build", "reference", "map", "iterate", "validate", "export"])
    assert (tmp_path / "result" / "cg_params.json").exists()
    assert (tmp_path / "result" / "cg_bulk.top").exists()
    # re-running skips finished stages
    again = CGParamJob(tmp_path)
    with contextlib.redirect_stdout(io.StringIO()):
        again.run()
    assert len((tmp_path / "cg" / "history.json").read_text()) > 10
    # the exported system is accepted by grompp
    G.Gmx().grompp("cg_npt.mdp", "cg_bulk.gro", "cg_bulk.top", "chk.tpr", tmp_path / "result")


@needs_gmx
@pytest.mark.parametrize("mixing", ["free", "lb"])
def test_optimiser_recovers_known_parameters(tmp_path, mixing):
    from lignoforge.cgfit.selftest import run_selftest
    with contextlib.redirect_stdout(io.StringIO()):
        r = run_selftest(tmp_path / "st", iterations=7, iter_ps=800.0, mixing=mixing)
    assert r["rmse_last"] < 0.5 * r["rmse_first"]
    assert r["rmse_last"] < 0.1
    assert r["max_rel_sigma_error"] < 0.1
