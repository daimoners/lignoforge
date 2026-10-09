"""
FastAPI application serving the LignoForge web interface (local use only).

The server is a thin HTTP layer over :class:`lignoforge.api.Project`: every
endpoint maps to one method, so the GUI, the Python API and the CLI cannot
diverge.  It is meant to listen on ``127.0.0.1``; there is no authentication.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from lignoforge.api import NotFoundError, Project, ProjectError, system_check
from lignoforge.analysis import AnalysisError

API = "/api/v1"
STATIC_DIR = Path(__file__).parent / "static"
_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.\- ]{0,63}")


# ── request models (also drive the OpenAPI schema → TypeScript types) ─────────

class ProjectIn(BaseModel):
    name: str = Field(..., description="Project name (letters, digits, _ . - and spaces)")


class ChainSpecIn(BaseModel):
    name: Optional[str] = None
    n_monomers: Optional[int] = Field(None, ge=1, le=200)
    monomers: Dict[str, float] = Field(default_factory=dict, description="H/G/S fractions")
    linkages: Dict[str, float] = Field(default_factory=dict, description="linkage fractions")
    branching: float = Field(0.0, ge=0.0, le=1.0)
    seed: int = 42
    optimize: bool = True
    max_iter: int = Field(300, ge=10, le=5000)
    experimental_input: Optional[Dict[str, Any]] = None


class RunIn(BaseModel):
    chain_ids: List[str]
    kind: str = Field(..., pattern="^(em|md|cg)$")
    params: Dict[str, Any] = Field(default_factory=dict)
    name: Optional[str] = None


class AnalysisIn(BaseModel):
    run_id: str
    observable: str
    params: Dict[str, Any] = Field(default_factory=dict)


class CGParamsIn(BaseModel):
    name: str = "custom"
    parameters: Dict[str, Any]


class CGFitIn(BaseModel):
    run_ids: List[str]
    temperature: float = 300.0
    stage: str = "prod"
    stride: int = Field(1, ge=1)
    min_samples: int = Field(200, ge=2)
    name: str = "fitted"


# ── observable catalogue (drives the analysis forms in the GUI) ───────────────

OBSERVABLE_INFO = [
    {"id": "rg", "label": "Radius of gyration", "needs": "trajectory", "options": {}},
    {"id": "end_to_end", "label": "End-to-end distance", "needs": "trajectory", "options": {}},
    {"id": "rmsd", "label": "RMSD (heavy atoms)", "needs": "trajectory", "options": {}},
    {"id": "density", "label": "Box density", "needs": "trajectory", "options": {}},
    {"id": "contacts", "label": "Inter-monomer contacts", "needs": "trajectory",
     "options": {"cutoff": {"type": "number", "default": 4.5, "unit": "Å"}}},
    {"id": "rdf", "label": "Radial distribution function", "needs": "trajectory",
     "options": {"group_a": {"type": "string", "default": "solute"},
                 "group_b": {"type": "string", "default": "water_O"},
                 "rmax": {"type": "number", "default": 15.0, "unit": "Å"},
                 "nbins": {"type": "integer", "default": 150}}},
    {"id": "energy", "label": "Energies, T, P", "needs": "edr", "options": {}},
]
COMMON_OPTIONS = {
    "stage": {"type": "string", "default": None},
    "molecule": {"type": "integer", "default": 0},
    "start": {"type": "integer", "default": 0},
    "stop": {"type": "integer", "default": None},
    "step": {"type": "integer", "default": 1},
}


def create_app(workspace: Optional[str] = None, dev: bool = False) -> FastAPI:
    """Create the application serving projects stored under *workspace*."""
    ws = Path(workspace or os.environ.get("LIGNOFORGE_WORKSPACE")
              or Path.home() / "LignoForge").expanduser().resolve()
    ws.mkdir(parents=True, exist_ok=True)

    app = FastAPI(title="LignoForge", version="1", docs_url=f"{API}/docs",
                  openapi_url=f"{API}/openapi.json")
    app.state.workspace = ws
    if dev:
        app.add_middleware(CORSMiddleware, allow_origins=["http://localhost:5173"],
                           allow_methods=["*"], allow_headers=["*"])

    @app.exception_handler(NotFoundError)
    async def _nf(_, e):                     # noqa: ANN001
        return JSONResponse({"detail": str(e)}, status_code=404)

    @app.exception_handler(ProjectError)
    async def _pe(_, e):                     # noqa: ANN001
        return JSONResponse({"detail": str(e)}, status_code=400)

    @app.exception_handler(AnalysisError)
    async def _ae(_, e):                     # noqa: ANN001
        return JSONResponse({"detail": str(e)}, status_code=400)

    def project(name: str) -> Project:
        if not _NAME.fullmatch(name) or ".." in name:
            raise HTTPException(400, "Invalid project name")
        root = ws / name
        if not (root / "project.json").exists():
            raise HTTPException(404, f"Unknown project: {name}")
        return Project.open(root)

    # ── system ────────────────────────────────────────────────────────────────

    @app.get(f"{API}/system")
    def system() -> dict:
        return {**system_check(), "workspace": str(ws)}

    @app.get(f"{API}/catalog")
    def catalog() -> dict:
        from lignoforge.core.rules import linkage_names, monomer_types
        from lignoforge.core.polymer import ENABLE_BETA_1_LINKAGE
        schema = json.loads((Path(__file__).parents[1] / "io" / "lignin_info_schema.json").read_text())
        return {
            "monomers": list(monomer_types),
            "linkages": [l for l in linkage_names if l != "beta-1" or ENABLE_BETA_1_LINKAGE],
            "observables": OBSERVABLE_INFO, "common_options": COMMON_OPTIONS,
            "run_kinds": ["em", "md", "cg"], "input_schema": schema,
        }

    @app.post(f"{API}/resolve-spec")
    def resolve(spec: ChainSpecIn) -> dict:
        """Preview the explicit fractions a specification resolves to (no build)."""
        from lignoforge.core.builder import resolve_spec
        try:
            return resolve_spec(spec.model_dump(exclude_none=True))
        except Exception as e:                # noqa: BLE001
            raise HTTPException(400, str(e))

    # ── projects ──────────────────────────────────────────────────────────────

    @app.get(f"{API}/projects")
    def list_projects() -> List[dict]:
        out = []
        for d in sorted(ws.iterdir()):
            if (d / "project.json").exists():
                out.append(Project.open(d).info())
        return out

    @app.post(f"{API}/projects", status_code=201)
    def create_project(body: ProjectIn) -> dict:
        if not _NAME.fullmatch(body.name):
            raise HTTPException(400, "Invalid project name")
        return Project.create(ws / body.name, body.name).info()

    @app.get(f"{API}/projects/{{p}}")
    def get_project(p: str) -> dict:
        return project(p).info()

    # ── chains ────────────────────────────────────────────────────────────────

    @app.get(f"{API}/projects/{{p}}/chains")
    def list_chains(p: str) -> List[dict]:
        return project(p).chains()

    @app.post(f"{API}/projects/{{p}}/chains", status_code=202)
    def build_chain(p: str, body: ChainSpecIn) -> dict:
        spec = body.model_dump(exclude_none=True)
        return project(p).build_chain(spec, name=spec.pop("name", None), wait=False)

    @app.get(f"{API}/projects/{{p}}/chains/{{c}}")
    def get_chain(p: str, c: str) -> dict:
        return project(p).chain(c)

    @app.delete(f"{API}/projects/{{p}}/chains/{{c}}", status_code=204)
    def delete_chain(p: str, c: str) -> None:
        project(p).delete_chain(c)

    @app.get(f"{API}/projects/{{p}}/chains/{{c}}/structure")
    def chain_structure(p: str, c: str) -> FileResponse:
        return FileResponse(project(p).structure_path(c), media_type="chemical/x-pdb",
                            filename=f"{c}.pdb")

    @app.get(f"{API}/projects/{{p}}/chains/{{c}}/graph")
    def chain_graph(p: str, c: str) -> dict:
        return project(p).chain_graph(c)

    @app.get(f"{API}/projects/{{p}}/chains/{{c}}/topology")
    def chain_topology(p: str, c: str, terms: bool = False) -> dict:
        return project(p).topology(c, include_terms=terms)

    @app.get(f"{API}/projects/{{p}}/chains/{{c}}/charge-report")
    def chain_charges(p: str, c: str) -> dict:
        return project(p).charge_report(c)

    @app.post(f"{API}/projects/{{p}}/chains/{{c}}/topology/files")
    def write_topology(p: str, c: str) -> dict:
        return project(p).write_topology(c)

    @app.get(f"{API}/projects/{{p}}/chains/{{c}}/topology/files/{{name}}")
    def topology_file(p: str, c: str, name: str) -> FileResponse:
        pr = project(p)
        f = (pr.root / "chains" / c / "topology" / name).resolve()
        if f.parent != (pr.root / "chains" / c / "topology").resolve() or not f.is_file():
            raise HTTPException(404, "No such file")
        return FileResponse(f, filename=name)

    # ── runs ──────────────────────────────────────────────────────────────────

    @app.get(f"{API}/projects/{{p}}/runs")
    def list_runs(p: str) -> List[dict]:
        return project(p).runs()

    @app.post(f"{API}/projects/{{p}}/runs", status_code=202)
    def start_run(p: str, body: RunIn) -> dict:
        return project(p).start_run(body.chain_ids, body.kind, body.params, body.name)

    @app.get(f"{API}/projects/{{p}}/runs/{{r}}")
    def get_run(p: str, r: str) -> dict:
        return project(p).run(r)

    @app.get(f"{API}/projects/{{p}}/runs/{{r}}/log", response_class=PlainTextResponse)
    def run_log(p: str, r: str, tail: int = Query(200, ge=1, le=5000)) -> str:
        return project(p).run_log(r, tail)

    @app.get(f"{API}/projects/{{p}}/runs/{{r}}/files")
    def run_files(p: str, r: str) -> List[dict]:
        return project(p).run_files(r)

    @app.get(f"{API}/projects/{{p}}/runs/{{r}}/files/{{name}}")
    def run_file(p: str, r: str, name: str) -> FileResponse:
        return FileResponse(project(p).run_file_path(r, name), filename=name)

    @app.post(f"{API}/projects/{{p}}/runs/{{r}}/cancel")
    def cancel_run(p: str, r: str) -> dict:
        return project(p).cancel_run(r)

    @app.delete(f"{API}/projects/{{p}}/runs/{{r}}", status_code=204)
    def delete_run(p: str, r: str) -> None:
        project(p).delete_run(r)

    # ── analyses ──────────────────────────────────────────────────────────────

    @app.get(f"{API}/projects/{{p}}/analyses")
    def list_analyses(p: str) -> List[dict]:
        return project(p).analyses()

    @app.post(f"{API}/projects/{{p}}/analyses", status_code=202)
    def run_analysis(p: str, body: AnalysisIn) -> dict:
        return project(p).analyze(body.run_id, body.observable, body.params, wait=False)

    @app.get(f"{API}/projects/{{p}}/analyses/{{a}}")
    def get_analysis(p: str, a: str) -> dict:
        return project(p).analysis(a)

    @app.delete(f"{API}/projects/{{p}}/analyses/{{a}}", status_code=204)
    def delete_analysis(p: str, a: str) -> None:
        project(p).delete_analysis(a)

    # ── coarse-grained parameters ─────────────────────────────────────────────

    @app.get(f"{API}/projects/{{p}}/cg")
    def list_cg(p: str) -> List[dict]:
        return project(p).cg_parameter_sets()

    @app.get(f"{API}/cg/defaults")
    def cg_defaults() -> dict:
        from lignoforge.cg import default_parameters
        return default_parameters().to_dict()

    @app.post(f"{API}/projects/{{p}}/cg", status_code=201)
    def save_cg(p: str, body: CGParamsIn) -> dict:
        return project(p).save_cg_parameters(body.parameters, name=body.name)

    @app.post(f"{API}/projects/{{p}}/cg/fit", status_code=201)
    def fit_cg(p: str, body: CGFitIn) -> dict:
        return project(p).fit_cg_parameters(
            body.run_ids, body.temperature, body.stage, body.stride,
            body.min_samples, body.name)

    @app.get(f"{API}/projects/{{p}}/cg/{{pid}}")
    def get_cg(p: str, pid: str) -> dict:
        return project(p).cg_parameters(pid)

    # ── static front-end (built bundle) ───────────────────────────────────────

    if (STATIC_DIR / "index.html").exists():
        app.mount("/assets", StaticFiles(directory=STATIC_DIR / "assets"), name="assets") \
            if (STATIC_DIR / "assets").is_dir() else None

        @app.get("/{path:path}", include_in_schema=False)
        def spa(path: str) -> FileResponse:      # client-side routing fallback
            if path.startswith("api/"):          # unknown API route: a real 404, not the app shell
                raise HTTPException(404, "Not found")
            f = (STATIC_DIR / path).resolve()
            if path and f.is_file() and STATIC_DIR in f.parents:
                return FileResponse(f)
            return FileResponse(STATIC_DIR / "index.html")
    else:
        @app.get("/", include_in_schema=False)
        def landing() -> PlainTextResponse:
            return PlainTextResponse(
                "LignoForge API is running. The web interface bundle is not installed "
                f"in this build.\nAPI documentation: {API}/docs\n")

    return app
