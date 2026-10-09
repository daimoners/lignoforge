"""Tests for the project API (chains, runs, analyses) and the HTTP layer."""

from __future__ import annotations

import json
import shutil
import time

import pytest

from lignoforge.api import NotFoundError, Project, ProjectError

GMX = shutil.which("gmx") is not None
needs_gmx = pytest.mark.skipif(not GMX, reason="GROMACS not installed")

SPEC = {"n_monomers": 5, "monomers": {"G": 1, "S": 1},
        "linkages": {"beta-O-4": 3, "4-O-5": 1, "beta-beta": 1}, "seed": 3}


def wait_run(p: Project, rid: str, timeout: float = 240):
    t0 = time.time()
    while time.time() - t0 < timeout:
        r = p.run(rid)
        if r["status"] != "running":
            return r
        time.sleep(0.5)
    raise TimeoutError(rid)


@pytest.fixture(scope="module")
def proj(tmp_path_factory):
    p = Project.create(tmp_path_factory.mktemp("ws") / "demo", "demo")
    p.build_chain(SPEC, name="five-mer")
    return p


# ── project / chains ──────────────────────────────────────────────────────────

def test_create_open_and_duplicate(tmp_path):
    p = Project.create(tmp_path / "x", "x")
    assert Project.open(tmp_path / "x").info()["name"] == "x"
    with pytest.raises(ProjectError):
        Project.create(tmp_path / "x")
    assert Project.create(tmp_path / "x", exist_ok=True).root == p.root
    with pytest.raises(ProjectError):
        Project.open(tmp_path / "missing")


def test_chain_record_is_complete_and_reproducible(proj):
    rec = proj.chain("c0001")
    assert rec["status"] == "ready" and rec["name"] == "five-mer"
    assert rec["stats"]["monomer_count"] == 5 and rec["stats"]["max_degree"] <= 2
    assert rec["resolved"]["seed"] == 3
    assert sum(rec["resolved"]["monomers"].values()) == pytest.approx(1.0)
    again = proj.build_chain(SPEC)
    assert again["stats"]["smiles"] == rec["stats"]["smiles"]       # same seed → same chain
    assert proj.structure_path("c0001").read_text().startswith("HEADER")
    g = proj.chain_graph("c0001")
    assert len(g["nodes"]) == 5 and len(g["edges"]) == 4


def test_invalid_specs_are_rejected(proj):
    with pytest.raises(ProjectError):
        proj.build_chain({"n_monomers": 0})
    with pytest.raises(ProjectError):
        proj.build_chain({"n_monomers": 5, "linkages": {"made-up": 1}})
    with pytest.raises(ProjectError):
        proj.build_chain({"n_monomers": 5, "monomers": {"X": 1}})


def test_beta_1_is_reported_not_silently_dropped(proj):
    rec = proj.build_chain({"n_monomers": 3, "linkages": {"beta-O-4": 1, "beta-1": 1}})
    assert any("beta-1" in w for w in rec["warnings"])
    assert rec["resolved"]["linkages"]["beta-1"] == 0.0


def test_experimental_input_drives_composition(proj):
    exp = {"material_origin": {"biomass_type": "softwood"},
           "extraction_process": {"process_type": "kraft"}}
    rec = proj.build_chain({"experimental_input": exp, "n_monomers": 4, "seed": 1})
    assert rec["status"] == "ready"
    m = rec["resolved"]["monomers"]
    assert m["G"] > m["S"]                       # softwood is guaiacyl-rich


def test_topology_and_charge_report(proj):
    t = proj.topology("c0001", include_terms=True)
    assert t["net_charge"] == 0.0 and t["n_atoms"] == len(t["atoms"])
    assert {"opls", "charge", "raw_charge"} <= set(t["atoms"][0])
    assert len(t["angles"]) == t["n_angles"]
    rep = proj.charge_report("c0001")
    assert rep["net_charge"] == 0.0 and rep["by_type"] and "text" in rep
    files = proj.write_topology("c0001")
    assert (proj.root / files["top"]).exists()


def test_unknown_ids_and_stale_builds(proj):
    with pytest.raises(NotFoundError):
        proj.chain("c9999")
    rec = proj.build_chain({"n_monomers": 2})
    p = proj.root / "chains" / rec["id"] / "record.json"
    d = json.loads(p.read_text())
    d["status"] = "building"                     # simulate a crash mid-build
    p.write_text(json.dumps(d))
    assert proj.chain(rec["id"])["status"] == "failed"


def test_async_build_reports_progress_states(proj):
    rec = proj.build_chain({"n_monomers": 3, "seed": 9}, wait=False)
    assert rec["status"] in ("building", "ready")
    for _ in range(120):
        if proj.chain(rec["id"])["status"] != "building":
            break
        time.sleep(0.5)
    assert proj.chain(rec["id"])["status"] == "ready"


# ── runs and analyses ─────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def md_run(proj):
    if not GMX:
        pytest.skip("GROMACS not installed")
    r = proj.start_run("c0001", "md", {"solvent": "none", "nvt_ps": 2, "prod_ns": 0.01,
                                       "traj_ps": 0.5, "threads": 2})
    r = wait_run(proj, r["id"])
    assert r["status"] == "done", proj.run_log(r["id"])
    return r


@needs_gmx
def test_run_record_and_files(proj, md_run):
    assert md_run["exit_code"] == 0 and md_run["layout"][0]["n_atoms"] > 50
    assert isinstance(md_run["params"]["seed"], int)                  # recorded for reproducibility
    names = {f["name"] for f in proj.run_files(md_run["id"])}
    assert {"prod.xtc", "prod.tpr", "prod.edr", "run.sh"} <= names
    assert "GROMACS" in proj.run_log(md_run["id"], 500)
    with pytest.raises(ProjectError):
        proj.run_file_path(md_run["id"], "../run.json")      # no path traversal


@needs_gmx
def test_run_status_survives_reopen(proj, md_run):
    again = Project.open(proj.root).run(md_run["id"])
    assert again["status"] == "done"


@needs_gmx
def test_invalid_run_requests(proj):
    before = len(proj.runs())
    with pytest.raises(ProjectError):
        proj.start_run("c0001", "nonsense")
    with pytest.raises(ProjectError):
        proj.start_run(["c0001", "c0002"], "md")
    with pytest.raises(ProjectError):
        proj.start_run("c0001", "md", {"bogus": 1})
    with pytest.raises(NotFoundError):
        proj.start_run("c9999", "em")
    assert len(proj.runs()) == before                         # failed prep leaves no folder


@needs_gmx
@pytest.mark.parametrize("obs,key", [
    ("rg", "Rg (Å)"), ("end_to_end", "End-to-end (Å)"), ("rmsd", "RMSD (Å)"),
    ("density", "Density (g/cm³)"), ("energy", "Temperature"),
])
def test_series_analyses(proj, md_run, obs, key):
    a = proj.analyze(md_run["id"], obs, {"stage": "prod"} if obs == "energy" else {})
    assert a["status"] == "done", a["error"]
    s = a["result"]["summary"][key]
    assert s["min"] - 1e-9 <= s["mean"] <= s["max"] + 1e-9
    assert len(a["result"]["time_ps"]) == len(a["result"]["series"][key]) > 5
    if obs == "energy":
        assert 250 < s["mean"] < 350                        # thermostat at 300 K


@needs_gmx
def test_contacts_and_rdf(proj, md_run):
    a = proj.analyze(md_run["id"], "contacts", {"cutoff": 4.0})
    m = a["result"]["extra"]["matrix"]
    assert len(m) == 5 and all(m[i][i] == 0 for i in range(5))
    assert all(abs(m[i][j] - m[j][i]) < 1e-12 for i in range(5) for j in range(5))
    r = proj.analyze(md_run["id"], "rdf", {"group_a": "solute", "group_b": "solute_O"})
    assert r["status"] == "done" and len(r["result"]["extra"]["g"]) == 150
    bad = proj.analyze(md_run["id"], "rdf")                 # needs water; none in vacuum
    assert bad["status"] == "failed" and "water" in bad["error"].lower()


@needs_gmx
def test_analysis_errors(proj, md_run):
    with pytest.raises(ProjectError):
        proj.analyze(md_run["id"], "nonsense")
    assert proj.analyze(md_run["id"], "rg", {"molecule": 5})["status"] == "failed"
    assert proj.analyze(md_run["id"], "rg", {"stage": "nope"})["status"] == "failed"


@needs_gmx
def test_cancel_run(proj):
    r = proj.start_run("c0001", "md", {"solvent": "none", "nvt_ps": 20000, "threads": 1})
    time.sleep(3)
    assert proj.cancel_run(r["id"])["status"] in ("running", "cancelled")
    r = wait_run(proj, r["id"], 30)
    assert r["status"] == "cancelled"
    proj.delete_run(r["id"])
    with pytest.raises(NotFoundError):
        proj.run(r["id"])


@needs_gmx
def test_cg_run_fit_and_analysis(proj, md_run):
    fit = proj.fit_cg_parameters([md_run["id"]], stage="prod", min_samples=5)
    assert fit["id"].startswith("p")
    assert fit["report"]["fitted"] or fit["report"]["kept_default"]   # short run: either is valid
    params = proj.cg_parameters(fit["id"])
    assert all(v["k"] > 0 for v in params["bonds"].values())
    r = wait_run(proj, proj.start_run(
        "c0001", "cg", {"copies": 2, "run_ns": 0.2, "density": 0.1,
                        "cg_parameters": fit["id"], "threads": 1})["id"])
    assert r["status"] == "done", proj.run_log(r["id"])
    assert r["layout"] == [{"chain_id": "c0001", "n_atoms": 5, "copies": 2}]
    a = proj.analyze(r["id"], "rg", {"molecule": 1})
    assert a["status"] == "done", a["error"]
    with pytest.raises(ProjectError):
        proj.fit_cg_parameters(["r9999"])


# ── HTTP layer ────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def client(tmp_path_factory):
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient
    from lignoforge.web.app import create_app
    return TestClient(create_app(str(tmp_path_factory.mktemp("web"))))


def test_http_system_catalog_and_resolve(client):
    s = client.get("/api/v1/system").json()
    assert s["features"]["build_structures"] and "models" in s
    cat = client.get("/api/v1/catalog").json()
    assert cat["monomers"] == ["H", "G", "S"] and "beta-1" not in cat["linkages"]
    assert {o["id"] for o in cat["observables"]} >= {"rg", "rdf", "energy"}
    r = client.post("/api/v1/resolve-spec", json={"monomers": {"S": 1, "G": 1, "H": 0},
                                                  "linkages": {"beta-O-4": 1}})
    assert r.status_code == 200 and r.json()["monomers"]["S"] == 0.5
    assert client.post("/api/v1/resolve-spec", json={"linkages": {"zzz": 1}}).status_code == 400


def test_http_full_chain_flow(client):
    assert client.post("/api/v1/projects", json={"name": "p one"}).status_code == 201
    assert client.post("/api/v1/projects", json={"name": "../evil"}).status_code == 400
    assert client.get("/api/v1/projects/p one").json()["n_chains"] == 0
    r = client.post("/api/v1/projects/p one/chains",
                    json={"n_monomers": 3, "monomers": {"G": 1}, "linkages": {"beta-O-4": 1}})
    assert r.status_code == 202
    cid = r.json()["id"]
    for _ in range(120):
        c = client.get(f"/api/v1/projects/p one/chains/{cid}").json()
        if c["status"] != "building":
            break
        time.sleep(0.5)
    assert c["status"] == "ready"
    base = f"/api/v1/projects/p one/chains/{cid}"
    assert client.get(base + "/structure").text.startswith("HEADER")
    assert client.get(base + "/topology").json()["net_charge"] == 0.0
    assert client.get(base + "/charge-report").json()["by_type"]
    assert client.get(base + "/graph").json()["edges"]
    files = client.post(base + "/topology/files").json()
    assert client.get(base + "/topology/files/lig.top").status_code == 200
    assert client.get(base + "/topology/files/..%2F..%2Frecord.json").status_code in (400, 404)
    assert client.delete(base).status_code == 204
    assert client.get(base).status_code == 404


def test_http_errors_map_to_status_codes(client):
    client.post("/api/v1/projects", json={"name": "errs"})
    assert client.get("/api/v1/projects/nope").status_code == 404
    assert client.post("/api/v1/projects/errs/chains", json={"n_monomers": 0}).status_code == 422
    assert client.post("/api/v1/projects/errs/runs",
                       json={"chain_ids": ["c0001"], "kind": "zzz"}).status_code == 422
    assert client.get("/api/v1/projects/errs/runs/r0001").status_code == 404
    assert client.post("/api/v1/projects/errs/analyses",
                       json={"run_id": "r1", "observable": "rg"}).status_code == 404


@needs_gmx
def test_http_run_and_analysis_flow(client):
    pytest.importorskip("MDAnalysis")
    client.post("/api/v1/projects", json={"name": "sim"})
    c = client.post("/api/v1/projects/sim/chains",
                    json={"n_monomers": 3, "monomers": {"G": 1}, "linkages": {"beta-O-4": 1}}).json()
    for _ in range(120):
        if client.get(f"/api/v1/projects/sim/chains/{c['id']}").json()["status"] != "building":
            break
        time.sleep(0.5)
    r = client.post("/api/v1/projects/sim/runs", json={
        "chain_ids": [c["id"]], "kind": "md",
        "params": {"solvent": "none", "nvt_ps": 1, "prod_ns": 0.004, "traj_ps": 0.5, "threads": 1}})
    assert r.status_code == 202
    rid = r.json()["id"]
    for _ in range(240):
        run = client.get(f"/api/v1/projects/sim/runs/{rid}").json()
        if run["status"] != "running":
            break
        time.sleep(0.5)
    assert run["status"] == "done"
    assert "GROMACS" in client.get(f"/api/v1/projects/sim/runs/{rid}/log?tail=50").text
    assert any(f["name"] == "prod.xtc" for f in client.get(f"/api/v1/projects/sim/runs/{rid}/files").json())
    assert client.get(f"/api/v1/projects/sim/runs/{rid}/files/..%2Frun.json").status_code in (400, 404)
    a = client.post("/api/v1/projects/sim/analyses", json={"run_id": rid, "observable": "rg"}).json()
    for _ in range(120):
        a = client.get(f"/api/v1/projects/sim/analyses/{a['id']}").json()
        if a["status"] != "running":
            break
        time.sleep(0.5)
    assert a["status"] == "done" and a["result"]["summary"]["Rg (Å)"]["mean"] > 3


def test_unknown_api_routes_are_404_not_the_app_shell(client):
    """With the front-end bundle installed, /api/... typos must not return index.html."""
    assert client.get("/api/v1/nope").status_code == 404
    assert client.get("/api/v1/projects/x/nope/y").status_code == 404
