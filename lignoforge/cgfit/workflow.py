"""
CG parametrisation workflow for bulk lignin.

Stages (each is resumable; finished work is skipped on re-run)::

    build      atomistic chains of the ensemble  → chains/chain_*.json
    reference  pack, minimise, compress, equilibrate and run the atomistic bulk (OPLS-AA)
    map        project the reference on the beads; Boltzmann-invert bonded terms;
               reference g(r); initial non-bonded parameters
    iterate    CG MD at the reference density ↔ relative-entropy / pressure updates
    validate   CG run at constant pressure; density, g(r), chain size vs reference
    export     final parameters, CG bulk input, report

Everything is written below the work directory (``config.json``, ``state.json``,
``progress.json``, ``log.txt``); nothing is written to the repository.
"""

from __future__ import annotations

import copy
import json
import pickle
import time
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np

from lignoforge.cgfit import fit as F
from lignoforge.cgfit import gmxio as G
from lignoforge.cgfit.model import MASS, CGParams, SCHEMES, map_chain
from lignoforge.cgfit.sampling import SystemMap, Stats, ideal_pair_numbers, iter_frames, sample

STAGES = ["build", "reference", "map", "iterate", "validate", "export"]

DEFAULT_CONFIG: dict = {
    "name": "lignin-bulk",
    "scheme": "ring-side",
    "temperature": 500.0,                 # K, reference and CG temperature (melt)
    "seed": 1,
    "gmx": None,                          # path of the gmx binary (None: $GMX or PATH)
    "ensemble": {
        "spec": {"n_monomers": 6, "monomers": {"H": 0.05, "G": 0.55, "S": 0.40},
                 "linkages": {"beta-O-4": 0.5, "5-5": 0.1, "beta-5": 0.15, "beta-beta": 0.1,
                              "4-O-5": 0.1, "alpha-O-4": 0.05},
                 "seed": 11},
        "n_distinct": 4,                  # distinct chain sequences (seed, seed+1, ...)
        "copies": 6,                      # copies of each
        "chains": None,                   # or an explicit list of chain specs (overrides spec)
    },
    "atomistic": {
        "start_density": 0.4,             # g/cm3 at packing
        "compress_ps": 500.0, "compress_pressure": 500.0,
        "equil_ps": 2000.0, "prod_ps": 3000.0, "out_ps": 10.0, "dt": 0.002,
        "mdrun_args": [],
    },
    "cg": {
        "cutoff": 1.2, "dt": 0.005, "iterations": 15, "iter_ps": 2000.0, "discard_frac": 0.25,
        "out_ps": 5.0, "tau_t": 1.0, "mixing": "free", "pressure_gain": 0.05,
        "trust_eps": 0.15, "trust_sigma": 0.04, "damping": 0.3, "min_samples": 300,
        "validate_ps": 5000.0, "mdrun_args": [],
        "rdf_tol": 0.05, "pressure_tol": 200.0,           # convergence: overall RMSE(g), |P-Pref| bar
    },
    "tolerances": {"density_rel": 0.05, "rdf_rmse": 0.25, "rg_rel": 0.10},
}


def _merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def quick_config() -> dict:
    """Seconds-to-minutes settings for tool validation (not for production)."""
    return _merge(DEFAULT_CONFIG, {
        "ensemble": {"spec": {"n_monomers": 4}, "n_distinct": 2, "copies": 8},
        "atomistic": {"compress_ps": 20.0, "equil_ps": 20.0, "prod_ps": 40.0, "out_ps": 2.0,
                      "start_density": 0.35},
        "cg": {"cutoff": 0.8, "iterations": 3, "iter_ps": 100.0, "out_ps": 1.0, "validate_ps": 100.0,
               "min_samples": 50},
    })


class CGParamJob:
    def __init__(self, workdir, config: Optional[dict] = None, gmx: Optional[str] = None,
                 dry_run: bool = False, log: Optional[Callable[[str], None]] = None):
        self.dir = Path(workdir)
        self.dir.mkdir(parents=True, exist_ok=True)
        cfg_path = self.dir / "config.json"
        if config is not None:
            self.cfg = _merge(DEFAULT_CONFIG, config)
            cfg_path.write_text(json.dumps(self.cfg, indent=2))
        elif cfg_path.exists():
            self.cfg = _merge(DEFAULT_CONFIG, json.loads(cfg_path.read_text()))
        else:
            self.cfg = copy.deepcopy(DEFAULT_CONFIG)
            cfg_path.write_text(json.dumps(self.cfg, indent=2))
        if self.cfg["scheme"] not in SCHEMES:
            raise ValueError(f"Unknown scheme '{self.cfg['scheme']}' (available: {sorted(SCHEMES)})")
        gmx = gmx or self.cfg.get("gmx")
        self.gmx_atom = G.Gmx(gmx, self.cfg["atomistic"]["mdrun_args"], dry_run)
        self.gmx_cg = G.Gmx(gmx, self.cfg["cg"]["mdrun_args"], dry_run)
        self.dry_run = dry_run
        self._log = log or (lambda s: None)
        self.state = self._read("state.json", {})

    # ── bookkeeping ───────────────────────────────────────────────────────────
    def _read(self, name, default=None):
        p = self.dir / name
        return json.loads(p.read_text()) if p.exists() else default

    def _write(self, name, obj):
        (self.dir / name).write_text(json.dumps(obj, indent=2, default=float))

    def log(self, msg: str) -> None:
        with open(self.dir / "log.txt", "a") as fh:
            fh.write(f"{time.strftime('%H:%M:%S')} {msg}\n")
        self._log(msg)

    def progress(self, stage: str, message: str = "", fraction: Optional[float] = None) -> None:
        self._write("progress.json", {"stage": stage, "message": message, "fraction": fraction,
                                      "time": time.time()})

    def _mark(self, stage, status, **extra):
        self.state[stage] = {"status": status, **extra}
        self._write("state.json", self.state)

    # ── chains / mapping ──────────────────────────────────────────────────────
    def _chain_specs(self) -> List[dict]:
        e = self.cfg["ensemble"]
        if e.get("chains"):
            return list(e["chains"])
        specs = []
        for i in range(int(e["n_distinct"])):
            s = copy.deepcopy(e["spec"])
            s["seed"] = int(s.get("seed", 1)) + i
            specs.append(s)
        return specs

    def _copies(self, n: int) -> List[int]:
        c = self.cfg["ensemble"]["copies"]
        return list(c) if isinstance(c, (list, tuple)) else [int(c)] * n

    def load_topologies(self):
        from lignoforge.forcefield import build_chain_topology
        paths = sorted((self.dir / "chains").glob("chain_*.json"))
        if not paths:
            raise RuntimeError("no chains: run stage 'build' first")
        topos = []
        for n, p in enumerate(paths):
            topos.append(build_chain_topology(json.loads(p.read_text()), name=f"L{n}"))
        return topos

    def system_map(self) -> SystemMap:
        topos = self.load_topologies()
        mols = [map_chain(t, self.cfg["scheme"]) for t in topos]
        return SystemMap(mols, self._copies(len(mols)))

    # ── driver ────────────────────────────────────────────────────────────────
    def run(self, stages: Optional[Sequence[str]] = None, force: bool = False) -> dict:
        stages = list(stages or STAGES)
        for s in stages:
            if s not in STAGES:
                raise ValueError(f"unknown stage '{s}' (available: {STAGES})")
        for s in STAGES:
            if s not in stages:
                continue
            if self.state.get(s, {}).get("status") == "done" and not force:
                self.log(f"stage {s}: already done (use --force to redo)")
                continue
            self.log(f"stage {s}: start")
            self.progress(s, "running")
            self._mark(s, "running")
            try:
                getattr(self, f"stage_{s}")()
            except Exception as e:                       # noqa: BLE001
                self._mark(s, "failed", error=f"{type(e).__name__}: {e}")
                self.progress(s, f"failed: {e}")
                self.log(f"stage {s}: FAILED {type(e).__name__}: {e}")
                raise
            self._mark(s, "done")
            self.log(f"stage {s}: done")
            if self.dry_run and s == "reference":
                break
        self.progress(stages[-1], "finished", 1.0)
        return self.state

    # ── stage: build ──────────────────────────────────────────────────────────
    def stage_build(self):
        from lignoforge.core.builder import build_chain_dict
        from lignoforge.forcefield import build_chain_topology
        d = self.dir / "chains"
        d.mkdir(exist_ok=True)
        specs = self._chain_specs()
        reports = []
        for n, spec in enumerate(specs):
            chain, rep = build_chain_dict(spec, optimize=True)
            (d / f"chain_{n}.json").write_text(json.dumps(chain))
            topo = build_chain_topology(chain, name=f"L{n}")
            topo.write_gro(str(d / f"chain_{n}.gro"), margin_nm=0.1)
            reports.append({"chain": n, "n_atoms": len(topo.atoms), "net_charge": topo.net_charge,
                            "resolved": rep["resolved"], "derivatization": rep["derivatization"],
                            "warnings": rep["warnings"]})
            self.log(f"chain {n}: {len(topo.atoms)} atoms")
        (d / "ensemble.json").write_text(json.dumps(reports, indent=2, default=float))
        sm = self.system_map()
        self.log(f"scheme {self.cfg['scheme']}: {sm.n_beads} beads, types {sm.types}")

    # ── stage: reference ──────────────────────────────────────────────────────
    def stage_reference(self):
        a, T = self.cfg["atomistic"], self.cfg["temperature"]
        topos = self.load_topologies()
        copies = self._copies(len(topos))
        ref = self.dir / "reference"
        ref.mkdir(exist_ok=True)
        g = self.gmx_atom
        G.write_atomistic_top(topos, copies, ref / "bulk.top")
        mass = sum(c * sum(MASS[x.element] for x in t.atoms) for t, c in zip(topos, copies))
        n_atoms = sum(c * len(t.atoms) for t, c in zip(topos, copies))
        final_edge = G.box_edge_for_density(mass, 1.1)
        if final_edge < 2.3:
            raise RuntimeError(
                f"the box would shrink to about {final_edge:.2f} nm at melt density, too small for the "
                "1.0 nm atomistic cut-off (needs about 2.3 nm or more): increase ensemble.copies or "
                "ensemble.n_distinct (production runs need roughly 25 chains of 6-10 monomers)")
        packed = None
        for dens in (a["start_density"], a["start_density"] * 0.7, a["start_density"] * 0.5):
            edge = G.box_edge_for_density(mass, dens)
            (ref / "packed.gro").unlink(missing_ok=True)
            for n, (t, c) in enumerate(zip(topos, copies)):
                src = str((self.dir / "chains" / f"chain_{n}.gro").resolve())
                args = ["insert-molecules", "-ci", src, "-nmol", str(c), "-o", "packed.gro",
                        "-rot", "xyz", "-try", "300", "-seed", str(self.cfg["seed"] + n)]
                args += ["-f", "packed.gro"] if n else ["-box", *[f"{edge:.4f}"] * 3]
                g.run(args, ref, log=f"insert_{n}.log")
            if g.dry_run:
                packed = True
                break
            if len(G.read_gro_xyz(ref / "packed.gro")) == n_atoms:
                packed = True
                self.log(f"packed {n_atoms} atoms at {dens:.2f} g/cm3 (box {edge:.2f} nm)")
                break
            self.log(f"packing incomplete at {dens:.2f} g/cm3, retrying at lower density")
        if not packed:
            raise RuntimeError("could not pack the chains; lower atomistic.start_density or copies")

        stages = [
            ("em", "em", 0.0, None), ("compress", "npt", a["compress_ps"], a["compress_pressure"]),
            ("equil", "npt", a["equil_ps"], 1.0), ("prod", "npt", a["prod_ps"], 1.0)]
        prev = "packed.gro"
        for k, (name, kind, ps, press) in enumerate(stages):
            d = ref / name
            d.mkdir(exist_ok=True)
            if (d / f"{name}.gro").exists():
                prev = f"{name}/{name}.gro"
                continue
            self.progress("reference", name, k / len(stages))
            G.write_text(d / f"{name}.mdp", G.atomistic_mdp(
                kind, ps, T, a["dt"], pressure=press or 1.0, out_ps=a["out_ps"],
                seed=self.cfg["seed"], continuation=name not in ("compress",)))
            g.grompp(f"{name}.mdp", f"../{prev}", "../bulk.top", f"{name}.tpr", d, maxwarn=1)
            g.mdrun(name, d)
            prev = f"{name}/{name}.gro"
        if g.dry_run:
            return
        prod = ref / "prod"
        g.run(["trjconv", "-s", "prod.tpr", "-f", "prod.xtc", "-pbc", "mol", "-o", "whole.xtc"],
              prod, stdin="0\n", log="trjconv.log")
        half = a["prod_ps"] / 2
        summary = {"pressure_bar": G.energy_average(g, "prod.edr", "Pressure", prod, half),
                   "density_g_cm3": G.energy_average(g, "prod.edr", "Density", prod, half),
                   "temperature": T, "n_atoms": n_atoms}
        self._write("reference/summary.json", summary)
        self.log(f"reference: density {summary['density_g_cm3']} g/cm3, P {summary['pressure_bar']} bar")

    # ── stage: map ────────────────────────────────────────────────────────────
    def stage_map(self):
        cg, T = self.cfg["cg"], self.cfg["temperature"]
        sm = self.system_map()
        prod = self.dir / "reference" / "prod"
        out = self.dir / "cg"
        out.mkdir(exist_ok=True)
        frames = iter_frames(str(prod / "prod.gro"), str(prod / "whole.xtc"))
        st = sample(sm, frames, mapped=False, rc=cg["cutoff"])
        self.log(f"reference mapped: {st.n_frames} frames, {sm.n_beads} beads")
        with open(out / "ref_stats.pkl", "wb") as fh:
            pickle.dump(st, fh)
        bonds, angles, rep = F.fit_bonded(st, T, cg["min_samples"])
        self._write("cg/bonded_report.json", rep)
        gr = st.gr(ideal_pair_numbers(sm))
        nb = F.initial_nonbonded(gr, st.r_centers)
        edge = float(np.mean(st.volumes) ** (1 / 3))
        if 2.2 * cg["cutoff"] > edge:
            raise RuntimeError(f"cutoff {cg['cutoff']} nm is too long for the box ({edge:.2f} nm; "
                               "needs box > 2.2 x cutoff): "
                               "increase the system size or lower cg.cutoff")
        params = CGParams(
            scheme=self.cfg["scheme"], temperature=T, cutoff=cg["cutoff"], nonbonded=nb,
            bonds=bonds, angles=angles,
            provenance=f"cgfit: bonded by Boltzmann inversion of the mapped OPLS-AA reference at {T} K; "
                       f"non-bonded by relative-entropy minimisation (initial guess).")
        params.save(str(out / "params_000.json"))
        # starting CG configuration: last mapped frame, in the mean box
        last = None
        for last in iter_frames(str(prod / "prod.gro"), str(prod / "whole.xtc")):
            pass
        x = sm.map_positions(last[0])
        types = [sm.types[t] for t in sm.bead_type]
        G.write_text(out / "start.gro", G.cg_gro_text(types, x, last[1], "cg start"))
        self._write("cg/reference_info.json", {"box_nm": last[1].tolist(),
                                              "mean_edge_nm": edge, "n_frames": st.n_frames})

    # ── stage: iterate ────────────────────────────────────────────────────────
    def _cg_run(self, d: Path, params: CGParams, gro_in: Path, kind: str, ps: float, seed: int,
                pressure: float = 1.0):
        sm = self.system_map()
        cg = self.cfg["cg"]
        d.mkdir(exist_ok=True)
        G.write_text(d / "cg.top", G.cg_top_text(sm.mols, sm.copies, params))
        G.write_text(d / "em.mdp", G.cg_mdp("em", 0, params.temperature, cg["cutoff"]))
        G.write_text(d / "md.mdp", G.cg_mdp(kind, ps, params.temperature, cg["cutoff"], cg["dt"],
                                            cg["out_ps"], cg["tau_t"], pressure, seed))
        g = self.gmx_cg
        g.grompp("em.mdp", str(gro_in.resolve()), "cg.top", "em.tpr", d)
        g.mdrun("em", d)
        g.grompp("md.mdp", "em.gro", "cg.top", "md.tpr", d)
        g.mdrun("md", d)

    def stage_iterate(self):
        cg, T = self.cfg["cg"], self.cfg["temperature"]
        sm = self.system_map()
        out = self.dir / "cg"
        with open(out / "ref_stats.pkl", "rb") as fh:
            ref: Stats = pickle.load(fh)
        info = self._read("cg/reference_info.json")
        p_ref = (self._read("reference/summary.json") or {}).get("pressure_bar") or 1.0
        ideal = ideal_pair_numbers(sm)
        ref_gr = ref.gr(ideal)
        weights = {p: ideal[p] for p in ref_gr}
        mixing = F.Mixing(sm.pair_names, cg["mixing"])
        history = self._read("cg/history.json", [])
        done = len(history)
        params = CGParams.load(str(out / f"params_{done:03d}.json"))
        gro = out / "start.gro" if done == 0 else out / f"iter_{done:02d}" / "md.gro"
        # the reference box is fixed for the NVT iterations
        box = info["box_nm"]
        self._set_box(gro, box)
        good = 0
        for it in range(done + 1, int(cg["iterations"]) + 1):
            self.progress("iterate", f"iteration {it}/{cg['iterations']}", (it - 1) / cg["iterations"])
            d = out / f"iter_{it:02d}"
            self._cg_run(d, params, gro, "nvt", cg["iter_ps"], self.cfg["seed"] + it)
            frames = iter_frames(str(d / "md.gro"), str(d / "md.xtc"),
                                 begin_ps=cg["discard_frac"] * cg["iter_ps"])
            st = sample(sm, frames, mapped=True, rc=cg["cutoff"])
            gr = st.gr(ideal)
            err = F.rdf_rmse(ref_gr, gr, ref.r_centers, weights)
            p_cg = G.energy_average(self.gmx_cg, "md.edr", "Pressure", d, cg["discard_frac"] * cg["iter_ps"])
            rec = {"iteration": it, "rdf_rmse": err, "pressure_bar": p_cg, "pressure_ref_bar": p_ref,
                   "nonbonded": copy.deepcopy(params.nonbonded)}
            self.log(f"iter {it}: RMSE(g)={err['overall']:.4f}  P={p_cg} bar (ref {p_ref:.0f})")
            converged = (err["overall"] < cg["rdf_tol"] and p_cg is not None
                         and abs(p_cg - p_ref) < cg["pressure_tol"])
            good = good + 1 if converged else 0
            new_nb, step = F.re_step(ref, st, params.nonbonded, mixing, cg["cutoff"], T,
                                     cg["trust_eps"], cg["trust_sigma"], cg["damping"])
            rec.update(step)
            if cg["pressure_gain"] and p_cg is not None:
                new_nb, f = F.pressure_step(new_nb, p_cg, p_ref, cg["pressure_gain"])
                rec["pressure_scale"] = f
            history.append(rec)
            self._write("cg/history.json", history)
            params = CGParams.from_dict({**params.to_dict(), "nonbonded": new_nb})
            params.save(str(out / f"params_{it:03d}.json"))
            gro = d / "md.gro"
            if good >= 2:
                self.log(f"converged after {it} iterations")
                break

    @staticmethod
    def _set_box(gro: Path, box) -> None:
        lines = Path(gro).read_text().splitlines()
        lines[-1] = "".join(f"{b:10.5f}" for b in box[:3])
        Path(gro).write_text("\n".join(lines) + "\n")

    # ── stage: validate ───────────────────────────────────────────────────────
    def stage_validate(self):
        cg, T, tol = self.cfg["cg"], self.cfg["temperature"], self.cfg["tolerances"]
        sm = self.system_map()
        out = self.dir / "cg"
        hist = self._read("cg/history.json", [])
        final = out / f"params_{len(hist):03d}.json"
        params = CGParams.load(str(final))
        with open(out / "ref_stats.pkl", "rb") as fh:
            ref: Stats = pickle.load(fh)
        ideal = ideal_pair_numbers(sm)
        d = self.dir / "validation"
        gro = out / "start.gro" if not hist else out / f"iter_{len(hist):02d}" / "md.gro"
        self._cg_run(d, params, gro, "npt", cg["validate_ps"], self.cfg["seed"] + 999, 1.0)
        st = sample(sm, iter_frames(str(d / "md.gro"), str(d / "md.xtc"),
                                    begin_ps=0.4 * cg["validate_ps"]), mapped=True, rc=cg["cutoff"])
        ref_gr, gr = ref.gr(ideal), st.gr(ideal)
        err = F.rdf_rmse(ref_gr, gr, ref.r_centers, {p: ideal[p] for p in ref_gr})
        rsum = self._read("reference/summary.json") or {}
        dens = G.energy_average(self.gmx_cg, "md.edr", "Density", d, 0.4 * cg["validate_ps"])
        rg = {n: {"reference": float(np.mean(ref.rg[n])), "cg": float(np.mean(st.rg[n]))}
              for n in ref.rg}
        for v in rg.values():
            v["rel_err"] = (v["cg"] - v["reference"]) / v["reference"]
        dens_ref = rsum.get("density_g_cm3")
        checks = {
            "density_rel_err": None if not (dens and dens_ref) else (dens - dens_ref) / dens_ref,
            "rdf_rmse": err["overall"],
            "rg_max_rel_err": max(abs(v["rel_err"]) for v in rg.values()),
        }
        checks["pass"] = {
            "density": checks["density_rel_err"] is not None and abs(checks["density_rel_err"]) <= tol["density_rel"],
            "rdf": checks["rdf_rmse"] <= tol["rdf_rmse"],
            "rg": checks["rg_max_rel_err"] <= tol["rg_rel"],
        }
        res = {"density_cg": dens, "density_reference": dens_ref, "rg": rg, "rdf_rmse": err,
               "checks": checks, "tolerances": tol}
        self._write("validation/validation.json", res)
        self.log(f"validation: {checks}")

    # ── stage: export ─────────────────────────────────────────────────────────
    def stage_export(self):
        out = self.dir / "result"
        out.mkdir(exist_ok=True)
        hist = self._read("cg/history.json", [])
        params = CGParams.load(str(self.dir / "cg" / f"params_{len(hist):03d}.json"))
        params.provenance += f" {len(hist)} relative-entropy iterations."
        params.save(str(out / "cg_params.json"))
        sm = self.system_map()
        G.write_text(out / "cg_bulk.top", G.cg_top_text(sm.mols, sm.copies, params))
        d = self.dir / "validation"
        if (d / "md.gro").exists():
            import shutil
            shutil.copy(d / "md.gro", out / "cg_bulk.gro")
            shutil.copy(d / "md.mdp", out / "cg_npt.mdp")
        val = self._read("validation/validation.json") or {}
        lines = [f"# CG parametrisation report — {self.cfg['name']}", "",
                 f"scheme `{self.cfg['scheme']}`, T = {self.cfg['temperature']} K, "
                 f"{len(hist)} iterations.", ""]
        if hist:
            last = hist[-1]
            lines += [f"Final RMSE(g) = {last['rdf_rmse']['overall']:.4f}; "
                      f"P = {last['pressure_bar']} bar (reference {last['pressure_ref_bar']}).", ""]
        if val:
            lines += ["## Validation (NPT)", "", "```", json.dumps(val["checks"], indent=2), "```", ""]
        (out / "report.md").write_text("\n".join(lines))
        self.log(f"exported to {out}")
