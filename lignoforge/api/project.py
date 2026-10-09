"""
Project: the single, UI-independent entry point to LignoForge.

A project is a folder::

    project.json            name, schema version, creation info
    chains/<id>/            record.json, chain.json, structure.pdb, graph.json,
                            topology/ (GROMACS files + charge report)
    runs/<id>/              run.json, log.txt, work/ (inputs and outputs)
    analyses/<id>/          analysis.json (parameters + results)
    cg/<id>/                parameters.json (+ fit report)

Everything the command line, notebooks and the web interface do goes through
this class, so all three behave identically and every result is reproducible
from the records stored next to it.
"""

from __future__ import annotations

import random
import shutil
import threading
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from lignoforge.api import runner
from lignoforge.api.store import (
    SCHEMA_VERSION, lignoforge_version, next_id, now, read_json, write_json,
)

__all__ = ["Project", "ProjectError", "NotFoundError"]

_POOL = ThreadPoolExecutor(max_workers=2, thread_name_prefix="lignoforge")
_LOCK = threading.RLock()
_ACTIVE: set = set()      # ids of chain builds / analyses currently in flight

RUN_KINDS = ("em", "md", "cg")
OBSERVABLES_NEEDING_TRAJECTORY = {"rg", "end_to_end", "density", "rdf", "contacts", "rmsd"}


class ProjectError(RuntimeError):
    """User-facing error (missing prerequisite, invalid parameters)."""


class NotFoundError(ProjectError):
    """An id (chain, run, analysis, ...) does not exist in the project."""


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple, set)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


class Project:
    # ── lifecycle ─────────────────────────────────────────────────────────────

    def __init__(self, root: Path):
        self.root = Path(root).resolve()
        if not (self.root / "project.json").exists():
            raise ProjectError(f"Not a LignoForge project: {self.root}")

    @classmethod
    def create(cls, root, name: Optional[str] = None, exist_ok: bool = False) -> "Project":
        root = Path(root)
        if (root / "project.json").exists():
            if exist_ok:
                return cls(root)
            raise ProjectError(f"Project already exists: {root}")
        root.mkdir(parents=True, exist_ok=True)
        write_json(root / "project.json", {
            "schema_version": SCHEMA_VERSION, "name": name or root.name,
            "created": now(), "lignoforge_version": lignoforge_version(),
        })
        for sub in ("chains", "runs", "analyses", "cg"):
            (root / sub).mkdir(exist_ok=True)
        return cls(root)

    @classmethod
    def open(cls, root) -> "Project":
        return cls(Path(root))

    def info(self) -> dict:
        d = read_json(self.root / "project.json")
        d.update(path=str(self.root),
                 n_chains=len(self._ids("chains")), n_runs=len(self._ids("runs")),
                 n_analyses=len(self._ids("analyses")), n_cg=len(self._ids("cg")))
        return d

    def _ids(self, kind: str) -> List[str]:
        return sorted(p.name for p in (self.root / kind).iterdir() if p.is_dir())

    def _dir(self, kind: str, rid: str) -> Path:
        d = self.root / kind / rid
        if not d.is_dir():
            raise NotFoundError(f"Unknown {kind[:-1]} id: {rid}")
        return d

    # ── chains ────────────────────────────────────────────────────────────────

    def build_chain(self, spec: dict, name: Optional[str] = None,
                    wait: bool = True) -> dict:
        """
        Build a chain from *spec* (see :func:`lignoforge.core.builder.resolve_spec`;
        extra keys: ``optimize`` (bool, default True), ``max_iter`` (int, 300)).

        With ``wait=False`` the record is returned immediately with
        ``status == "building"``; poll :meth:`chain` until ``ready``/``failed``.
        """
        from lignoforge.core.builder import resolve_spec

        try:
            resolved = resolve_spec(spec)
        except Exception as e:                       # noqa: BLE001
            raise ProjectError(f"Invalid chain specification: {e}") from e
        if resolved["n_monomers"] < 1 or resolved["n_monomers"] > 200:
            raise ProjectError("n_monomers must be between 1 and 200")
        cid = next_id(self.root / "chains", "c")
        rec = {
            "id": cid, "name": name or spec.get("name") or cid, "status": "building",
            "created": now(), "lignoforge_version": lignoforge_version(),
            "spec": _jsonable({k: v for k, v in spec.items()}), "resolved": resolved,
            "warnings": list(resolved["warnings"]), "error": None, "stats": {},
        }
        write_json(self.root / "chains" / cid / "record.json", rec)
        _ACTIVE.add(cid)
        fut = _POOL.submit(self._build_chain_worker, cid, spec, resolved)
        if wait:
            fut.result()
        return self.chain(cid)

    def _build_chain_worker(self, cid: str, spec: dict, resolved: dict) -> None:
        try:
            self._build_chain_inner(cid, spec, resolved)
        finally:
            _ACTIVE.discard(cid)

    def _build_chain_inner(self, cid: str, spec: dict, resolved: dict) -> None:
        from lignoforge.core.builder import grow_chain_exact
        from lignoforge.core.characterization import Characterize
        from lignoforge.core.rules import linkage_names, monomer_types
        from lignoforge.core.utils import graph_to_smile
        from lignoforge.structure.generator import MolecularStructureGenerator
        from lignoforge.structure.pdb import PDBStructureWriter

        d = self.root / "chains" / cid
        rec = read_json(d / "record.json")
        try:
            # (stdout is deliberately not redirected: redirect_stdout is
            # process-global and would swallow output of other threads.)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                polymer = grow_chain_exact(
                    resolved["n_monomers"],
                    [resolved["monomers"][k] for k in monomer_types],
                    [resolved["linkages"][k] for k in linkage_names],
                    resolved["seed"], resolved["branching"])
                chain = MolecularStructureGenerator().polymer_atomistic_topology(
                    polymer, random_seed=resolved["seed"],
                    optimize_3d=bool(spec.get("optimize", True)),
                    max_uff_iterations=int(spec.get("max_iter", 300)))
            rec["warnings"] += [str(w.message) for w in caught]
            write_json(d / "chain.json", chain, indent=None)
            PDBStructureWriter()._write_atomistic_topology_as_pdb(chain, str(d / "structure.pdb"))
            graph = {
                "nodes": [{"id": n + 1, "type": a["mtype"]} for n, a in polymer.bigG.nodes(data=True)],
                "edges": [{"source": u + 1, "target": v + 1, "linkage": a.get("btype")}
                          for u, v, a in polymer.bigG.edges(data=True)],
            }
            write_json(d / "graph.json", graph)
            stats = _jsonable(Characterize(polymer).summary())
            stats.update(n_atoms=chain["n_atoms"], smiles=graph_to_smile(polymer.G),
                         max_degree=max(dict(polymer.bigG.degree()).values()))
            rec.update(status="ready", stats=stats)
            if stats.get("monomer_count", resolved["n_monomers"]) < resolved["n_monomers"]:
                rec["warnings"].append(
                    f"Requested {resolved['n_monomers']} monomers, built "
                    f"{stats.get('monomer_count')}.")
        except Exception as e:                         # noqa: BLE001
            rec.update(status="failed", error=f"{type(e).__name__}: {e}")
        write_json(d / "record.json", rec)

    def chains(self) -> List[dict]:
        return [self.chain(c) for c in self._ids("chains")]

    def chain(self, cid: str) -> dict:
        rec = read_json(self._dir("chains", cid) / "record.json")
        if rec["status"] == "building" and cid not in _ACTIVE:
            rec["status"] = "failed"
            rec["error"] = "interrupted (application was restarted during the build)"
            write_json(self.root / "chains" / cid / "record.json", rec)
        return rec

    def chain_data(self, cid: str) -> dict:
        self._require_ready(cid)
        return read_json(self.root / "chains" / cid / "chain.json")

    def chain_graph(self, cid: str) -> dict:
        self._require_ready(cid)
        return read_json(self.root / "chains" / cid / "graph.json")

    def structure_path(self, cid: str) -> Path:
        self._require_ready(cid)
        return self.root / "chains" / cid / "structure.pdb"

    def _require_ready(self, cid: str) -> None:
        st = self.chain(cid)["status"]
        if st != "ready":
            raise ProjectError(f"Chain {cid} is not ready (status: {st})")

    def delete_chain(self, cid: str) -> None:
        shutil.rmtree(self._dir("chains", cid))

    # ── topology / force field ────────────────────────────────────────────────

    def topology(self, cid: str, include_terms: bool = False) -> dict:
        """Typed OPLS-AA topology (atoms with type/charge/raw charge, bonds, terms)."""
        from lignoforge.forcefield import build_chain_topology
        return build_chain_topology(self.chain_data(cid), name="lig").to_dict(include_terms)

    def write_topology(self, cid: str) -> Dict[str, str]:
        """Write the GROMACS files and the charge report into ``chains/<id>/topology``."""
        from lignoforge.forcefield import write_gromacs_system
        out = self.root / "chains" / cid / "topology"
        paths = write_gromacs_system(self.chain_data(cid), str(out), name="lig")
        return {k: str(Path(v).relative_to(self.root)) for k, v in paths.items()}

    def charge_report(self, cid: str) -> dict:
        from lignoforge.forcefield import build_chain_topology
        t = build_chain_topology(self.chain_data(cid), name="lig")
        return {"text": t.charge_report_text(), "residues": t.charge_report,
                "by_type": t.charge_adjustments_by_type(),
                "atoms": t.charge_adjustments(), "net_charge": t.net_charge}

    # ── simulation runs ───────────────────────────────────────────────────────

    def start_run(self, chain_ids: Sequence[str] | str, kind: str,
                  params: Optional[dict] = None, name: Optional[str] = None) -> dict:
        """
        Prepare and launch a local GROMACS run.

        ``kind``: ``"em"`` (minimisation), ``"md"`` (full atomistic workflow) or
        ``"cg"`` (coarse-grained system of all given chains).
        Common params: ``threads`` (0 = all), ``gmx`` (binary).  ``md``:
        ``temperature, solvent (tip3p|spce|none), nvt_ps, npt_ps, prod_ns,
        traj_ps, box_distance_nm``.  ``cg``: ``copies, density, temperature,
        run_ns, cg_parameters`` (id of a stored CG parameter set).
        """
        params = dict(params or {})
        if kind not in RUN_KINDS:
            raise ProjectError(f"Unknown run kind '{kind}' (use {RUN_KINDS})")
        ids = [chain_ids] if isinstance(chain_ids, str) else list(chain_ids)
        if not ids:
            raise ProjectError("At least one chain is required")
        if kind != "cg" and len(ids) != 1:
            raise ProjectError(f"Run kind '{kind}' takes exactly one chain")
        gmx = params.pop("gmx", "gmx")
        if kind in ("md", "cg"):
            # explicit seed → the run can be repeated exactly (given the same GROMACS build)
            params["seed"] = int(params.get("seed") or random.randrange(1, 2**31 - 1))
        if shutil.which(gmx) is None:
            raise ProjectError("GROMACS ('gmx') was not found on PATH; "
                               "install it to run simulations.")
        chains = {c: self.chain_data(c) for c in ids}

        rid = next_id(self.root / "runs", "r")
        run_dir = self.root / "runs" / rid
        work = run_dir / "work"
        work.mkdir()
        threads = int(params.pop("threads", 0))
        try:
            layout = self._write_run_inputs(kind, ids, chains, work, params)
        except Exception as e:                       # noqa: BLE001
            shutil.rmtree(run_dir, ignore_errors=True)
            raise ProjectError(f"Could not prepare run: {e}") from e

        rec = {
            "id": rid, "name": name or f"{kind}-{ids[0]}", "kind": kind,
            "chain_ids": ids, "params": _jsonable(params), "threads": threads,
            "gmx": gmx, "layout": layout, "status": "running", "created": now(),
            "started": now(), "pid": None, "exit_code": None, "error": None,
            "cancel_requested": False, "lignoforge_version": lignoforge_version(),
        }
        rec["pid"] = runner.start(work, run_dir / "log.txt", threads, gmx)
        write_json(run_dir / "run.json", rec)
        return self.run(rid)

    def _write_run_inputs(self, kind, ids, chains, work: Path, params: dict) -> list:
        from lignoforge.md import write_atomistic_workflow

        if kind == "em":
            from lignoforge.forcefield import write_gromacs_system
            write_gromacs_system(chains[ids[0]], str(work), name="system")
            (work / "run.sh").write_text(
                "set -euo pipefail\n"
                "$GMX grompp -f em.mdp -c system.gro -p system.top -o em.tpr\n"
                '$GMX mdrun -nt "$NT" -deffnm em\n')
            n = len(chains[ids[0]]["monomers"]) and sum(
                m["n_atoms"] for m in chains[ids[0]]["monomers"])
            return [{"chain_id": ids[0], "n_atoms": n, "copies": 1}]
        if kind == "md":
            allowed = {"temperature", "solvent", "nvt_ps", "npt_ps", "prod_ns",
                       "traj_ps", "box_distance_nm", "seed"}
            bad = set(params) - allowed
            if bad:
                raise ValueError(f"unknown md parameter(s): {sorted(bad)}")
            write_atomistic_workflow(chains[ids[0]], str(work), name="system", **params)
            n = sum(m["n_atoms"] for m in chains[ids[0]]["monomers"])
            return [{"chain_id": ids[0], "n_atoms": n, "copies": 1}]
        # cg
        from lignoforge.cg import CGParameters, write_cg_system
        allowed = {"copies", "density", "temperature", "run_ns", "cg_parameters", "seed"}
        bad = set(params) - allowed
        if bad:
            raise ValueError(f"unknown cg parameter(s): {sorted(bad)}")
        cgp = None
        if params.get("cg_parameters"):
            cgp = CGParameters.load(str(self._dir("cg", params["cg_parameters"]) / "parameters.json"))
        write_cg_system(
            [chains[c] for c in ids], str(work), name="cg", parameters=cgp,
            copies=int(params.get("copies", 1)),
            density_g_cm3=float(params.get("density", 0.3)),
            temperature=float(params.get("temperature", 300.0)),
            run_ns=float(params.get("run_ns", 100.0)), seed=int(params["seed"]))
        return [{"chain_id": c, "n_atoms": len(chains[c]["monomers"]),
                 "copies": int(params.get("copies", 1))} for c in ids]

    def runs(self) -> List[dict]:
        return [self.run(r) for r in self._ids("runs")]

    def run(self, rid: str) -> dict:
        d = self._dir("runs", rid)
        with _LOCK:
            rec = read_json(d / "run.json")
            before = rec["status"]
            rec = runner.resolve_status(rec, d / "work")
            if rec["status"] != before:
                rec["finished"] = now()
                write_json(d / "run.json", rec)
        rec["progress"] = runner.progress(d / "work") if rec["status"] == "running" else None
        return rec

    def run_log(self, rid: str, tail: int = 200) -> str:
        log = self._dir("runs", rid) / "log.txt"
        if not log.exists():
            return ""
        return "\n".join(log.read_text(errors="replace").splitlines()[-tail:])

    def run_files(self, rid: str) -> List[dict]:
        work = self._dir("runs", rid) / "work"
        return [{"name": p.name, "size": p.stat().st_size}
                for p in sorted(work.iterdir()) if p.is_file() and not p.name.startswith(".")]

    def run_file_path(self, rid: str, name: str) -> Path:
        p = (self._dir("runs", rid) / "work" / name).resolve()
        if p.parent != (self._dir("runs", rid) / "work").resolve() or not p.is_file():
            raise ProjectError(f"No such file in run {rid}: {name}")
        return p

    def cancel_run(self, rid: str) -> dict:
        d = self._dir("runs", rid)
        with _LOCK:
            rec = read_json(d / "run.json")
            if rec["status"] in runner.TERMINAL:
                return self.run(rid)
            rec["cancel_requested"] = True
            write_json(d / "run.json", rec)
        runner.cancel(rec.get("pid"))
        return self.run(rid)

    def delete_run(self, rid: str) -> None:
        if self.run(rid)["status"] == "running":
            raise ProjectError("Cancel the run before deleting it")
        shutil.rmtree(self._dir("runs", rid))

    # ── analysis ──────────────────────────────────────────────────────────────

    def analyze(self, run_id: str, observable: str, params: Optional[dict] = None,
                wait: bool = True) -> dict:
        """
        Compute *observable* on a finished run.

        ``params``: ``stage`` (default ``prod`` for md, ``md`` for cg, ``em``),
        ``molecule`` (index among the lignin molecules, default 0), ``start``,
        ``stop``, ``step`` (frames) plus observable-specific options
        (``cutoff``, ``group_a``/``group_b``/``rmax``/``nbins``, ``terms``).
        """
        from lignoforge.analysis import OBSERVABLES
        if observable not in OBSERVABLES:
            raise ProjectError(f"Unknown observable '{observable}' "
                               f"(available: {sorted(OBSERVABLES)})")
        run = self.run(run_id)
        if run["status"] not in ("done", "running"):
            raise ProjectError(f"Run {run_id} is {run['status']}; nothing to analyse")
        aid = next_id(self.root / "analyses", "a")
        rec = {"id": aid, "run_id": run_id, "observable": observable,
               "params": _jsonable(params or {}), "status": "running",
               "created": now(), "error": None, "result": None,
               "lignoforge_version": lignoforge_version()}
        write_json(self.root / "analyses" / aid / "analysis.json", rec)
        _ACTIVE.add(aid)
        fut = _POOL.submit(self._analysis_worker, aid)
        if wait:
            fut.result()
        return self.analysis(aid)

    def _analysis_worker(self, aid: str) -> None:
        try:
            self._analysis_inner(aid)
        finally:
            _ACTIVE.discard(aid)

    def _analysis_inner(self, aid: str) -> None:
        from lignoforge.analysis import (
            AnalysisError, OBSERVABLES, load_universe, locate_molecule, select_molecule)

        d = self.root / "analyses" / aid
        rec = read_json(d / "analysis.json")
        try:
            run = self.run(rec["run_id"])
            p = dict(rec["params"])
            cg = run["kind"] == "cg"
            stage = p.pop("stage", {"md": "prod", "cg": "md", "em": "em"}[run["kind"]])
            work = self.root / "runs" / run["id"] / "work"
            obs = rec["observable"]
            kw = {k: p[k] for k in p if k not in ("molecule",)}
            if obs == "energy":
                edr = work / f"{stage}.edr"
                if not edr.exists():
                    raise AnalysisError(f"No energy file for stage '{stage}'")
                result = OBSERVABLES[obs](edr=str(edr), gmx=run["gmx"], **kw)
            else:
                tpr, xtc = work / f"{stage}.tpr", work / f"{stage}.xtc"
                if not xtc.exists():
                    raise AnalysisError(
                        f"No trajectory for stage '{stage}' (is the run finished?)")
                u = load_universe(str(tpr), str(xtc))
                entry, start = locate_molecule(run["layout"], int(p.get("molecule", 0)))
                chain_data = self.chain_data(entry["chain_id"])
                sel = select_molecule(u, chain_data, entry, start, coarse_grained=cg)
                result = OBSERVABLES[obs](u=u, sel=sel, chain_data=chain_data, **kw)
                result["extra"]["n_frames_total"] = int(u.trajectory.n_frames)
            rec.update(status="done", result=_jsonable(result))
        except Exception as e:                           # noqa: BLE001
            rec.update(status="failed", error=f"{type(e).__name__}: {e}")
        write_json(d / "analysis.json", rec)

    def analyses(self) -> List[dict]:
        return [self.analysis(a) for a in self._ids("analyses")]

    def analysis(self, aid: str) -> dict:
        rec = read_json(self._dir("analyses", aid) / "analysis.json")
        if rec["status"] == "running" and aid not in _ACTIVE:
            rec.update(status="failed", error="interrupted")
        return rec

    def delete_analysis(self, aid: str) -> None:
        shutil.rmtree(self._dir("analyses", aid))

    # ── coarse-grained parameters ─────────────────────────────────────────────

    def cg_parameter_sets(self) -> List[dict]:
        out = []
        for c in self._ids("cg"):
            rec = read_json(self.root / "cg" / c / "record.json")
            out.append(rec)
        return out

    def save_cg_parameters(self, params: dict, name: str = "custom",
                           source: str = "manual", report: Optional[dict] = None) -> dict:
        from lignoforge.cg import CGParameters
        cp = CGParameters.from_dict(params)
        cid = next_id(self.root / "cg", "p")
        d = self.root / "cg" / cid
        cp.save(str(d / "parameters.json"))
        rec = {"id": cid, "name": name, "source": source, "created": now(),
               "provenance": cp.provenance, "report": _jsonable(report or {})}
        write_json(d / "record.json", rec)
        return rec

    def cg_parameters(self, pid: str) -> dict:
        from lignoforge.cg import CGParameters
        return CGParameters.load(str(self._dir("cg", pid) / "parameters.json")).to_dict()

    def default_cg_parameters(self) -> dict:
        from lignoforge.cg import default_parameters
        return default_parameters().to_dict()

    def fit_cg_parameters(self, run_ids: Sequence[str], temperature: float = 300.0,
                          stage: str = "prod", stride: int = 1, min_samples: int = 200,
                          name: str = "fitted") -> dict:
        """
        Boltzmann-invert the CG bond/angle distributions of finished atomistic
        runs (pooled over all runs) and store the result as a CG parameter set.
        """
        from lignoforge.analysis import load_universe
        from lignoforge.cg import BoltzmannFitter, map_topology
        from lignoforge.forcefield import build_chain_topology

        fitter = BoltzmannFitter()
        for rid in run_ids:
            run = self.run(rid)
            if run["kind"] != "md" or run["status"] != "done":
                raise ProjectError(f"Run {rid} must be a finished atomistic md run")
            work = self.root / "runs" / rid / "work"
            u = load_universe(str(work / f"{stage}.tpr"), str(work / f"{stage}.xtc"))
            entry = run["layout"][0]
            chain = self.chain_data(entry["chain_id"])
            topo = build_chain_topology(chain)
            ag = u.atoms[:entry["n_atoms"]]
            from MDAnalysis import transformations as trans
            u.trajectory.add_transformations(trans.unwrap(ag))
            frames = np.array([ag.positions / 10.0 for _ in u.trajectory[::stride]])
            fitter.add(map_topology(topo), frames)
        params, report = fitter.fit(temperature, min_samples=min_samples)
        return self.save_cg_parameters(params.to_dict(), name=name,
                                       source=f"boltzmann-inversion:{','.join(run_ids)}",
                                       report=report)
