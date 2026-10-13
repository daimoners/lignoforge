"""
Progress and time-remaining estimates for a CG parametrisation work directory,
computed from the files GROMACS and the workflow leave on disk (so it works for
a job running in another terminal, under Slurm, or from the web interface).

Estimates are rough: the current ``mdrun`` is extrapolated from its own speed;
the rest of the reference stage reuses that speed (same system); CG iterations
use the mean duration of the iterations already finished.
"""

from __future__ import annotations

import json
import re
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

_STARTED = re.compile(r"Started mdrun on rank 0 (.+)")
_STEP = re.compile(r"^\s+Step\s+Time\s*\n\s+(\d+)\s+([\d.]+)", re.M)
_NSTEPS = re.compile(r"^\s*nsteps\s*=\s*(\d+)", re.M)
_DT = re.compile(r"^\s*dt\s*=\s*([\d.eE+-]+)", re.M)
_SKIP = (".mdrun.log", ".grompp.log")
FRESH_S = 900          # a log untouched for this long is not "running"
REF_STAGES = ["em", "compress", "equil", "prod"]


def fmt_hms(seconds: Optional[float]) -> str:
    if seconds is None or seconds != seconds:
        return "?"
    s = int(max(0, seconds))
    d, s = divmod(s, 86400)
    h, s = divmod(s, 3600)
    m, s = divmod(s, 60)
    return (f"{d}d " if d else "") + f"{h:d}h {m:02d}m" if (d or h) else f"{m:d}m {s:02d}s"


def _tail(path: Path, nbytes: int = 20000) -> str:
    with open(path, "rb") as fh:
        fh.seek(0, 2)
        fh.seek(max(0, fh.tell() - nbytes))
        return fh.read().decode(errors="replace")


def current_mdrun(workdir, fresh_s: float = FRESH_S) -> Optional[dict]:
    """The mdrun that is being written right now (newest fresh ``*.log`` of a md run), if any."""
    wd = Path(workdir)
    best = None
    for p in wd.rglob("*.log"):
        if p.name.endswith(_SKIP) or not p.with_suffix(".tpr").exists():
            continue
        if best is None or p.stat().st_mtime > best.stat().st_mtime:
            best = p
    if best is None or time.time() - best.stat().st_mtime > fresh_s:
        return None
    head = best.read_text(errors="replace")[:40000]
    tail = _tail(best)
    if "Finished mdrun" in tail:
        return None
    ns, dt = _NSTEPS.search(head), _DT.search(head)
    st = _STEP.findall(tail)
    out = {"run": str(best.relative_to(wd)), "deffnm": best.stem, "step": None, "nsteps": None,
           "time_ps": None, "fraction": None, "elapsed_s": None, "eta_s": None, "ns_per_day": None}
    if ns:
        out["nsteps"] = int(ns.group(1))
    if st:
        out["step"], out["time_ps"] = int(st[-1][0]), float(st[-1][1])
    sm = _STARTED.search(head) or _STARTED.search(tail)
    if sm:
        try:
            t0 = datetime.strptime(" ".join(sm.group(1).split()), "%a %b %d %H:%M:%S %Y").timestamp()
            out["elapsed_s"] = max(0.0, best.stat().st_mtime - t0)
        except ValueError:
            pass
    if out["step"] and out["nsteps"]:
        out["fraction"] = min(1.0, out["step"] / out["nsteps"])
        if out["elapsed_s"] and out["step"] > 0:
            rate = out["step"] / out["elapsed_s"]
            out["eta_s"] = (out["nsteps"] - out["step"]) / rate
            if dt:
                out["ns_per_day"] = rate * float(dt.group(1)) * 86400 / 1000.0
    return out


def _mdp_ps(path: Path) -> Optional[float]:
    if not path.exists():
        return None
    t = path.read_text()
    ns, dt = _NSTEPS.search(t), _DT.search(t)
    return int(ns.group(1)) * float(dt.group(1)) if ns and dt else None


def overview(workdir) -> dict:
    """Stage table, current mdrun and an overall estimate for one work directory."""
    wd = Path(workdir)
    state = json.loads((wd / "state.json").read_text()) if (wd / "state.json").exists() else {}
    cfg = json.loads((wd / "config.json").read_text()) if (wd / "config.json").exists() else {}
    hist = json.loads((wd / "cg" / "history.json").read_text()) if (wd / "cg" / "history.json").exists() else []
    cur = current_mdrun(wd)
    res = {"stages": {s: state.get(s, {}).get("status", "pending")
                      for s in ["build", "reference", "map", "iterate", "validate", "export"]},
           "current": cur, "iterations_done": len(hist), "iterations_planned": None,
           "last_rmse": hist[-1]["rdf_rmse"]["overall"] if hist else None,
           "eta_s": None, "detail": ""}
    if not cfg:
        return res
    at, cg = cfg["atomistic"], cfg["cg"]
    res["iterations_planned"] = cg["iterations"]

    eta = 0.0
    known = True
    if res["stages"]["reference"] != "done":
        plan = {"em": 0.0, "compress": at["compress_ps"], "equil": at["equil_ps"], "prod": at["prod_ps"]}
        done_ps = sum(ps for s, ps in plan.items() if (wd / "reference" / s / f"{s}.gro").exists())
        total_ps = sum(plan.values())
        cur_ps = cur["time_ps"] if cur and cur["run"].startswith("reference") and cur["time_ps"] else 0.0
        if cur and cur["run"].startswith("reference") and cur["ns_per_day"]:
            speed = cur["ns_per_day"] * 1000.0 / 86400.0            # ps / s
            eta += max(0.0, total_ps - done_ps - cur_ps) / speed
            res["detail"] = (f"reference: {done_ps + cur_ps:.0f}/{total_ps:.0f} ps "
                             f"at {cur['ns_per_day']:.0f} ns/day")
        else:
            known = False
    n_left = max(0, cg["iterations"] - len(hist)) if res["stages"]["iterate"] != "done" else 0
    if n_left:
        dirs = sorted((wd / "cg").glob("iter_*/md.gro"), key=lambda p: p.stat().st_mtime)
        if len(dirs) >= 2:
            per = (dirs[-1].stat().st_mtime - dirs[0].stat().st_mtime) / (len(dirs) - 1)
        elif cur and cur["run"].startswith("cg/iter") and cur["fraction"] and cur["elapsed_s"]:
            per = cur["elapsed_s"] / cur["fraction"]
        else:
            per = None
        if per is not None and known:
            eta += per * n_left - (cur["elapsed_s"] if cur and cur["run"].startswith("cg/iter") and cur["elapsed_s"] else 0)
            res["detail"] = (res["detail"] + "; " if res["detail"] else "") + \
                f"CG: {len(hist)}/{cg['iterations']} iterations (≤ {fmt_hms(per)} each; may stop earlier if converged)"
        else:
            known = False
            res["detail"] = (res["detail"] + "; " if res["detail"] else "") + \
                f"CG part ({cg['iterations']} iterations) not estimated until the first iteration finishes"
    res["eta_s"] = eta if known and (n_left or res["stages"]["reference"] != "done") else None
    return res


def format_overview(workdir) -> str:
    o = overview(workdir)
    lines = [f"{Path(workdir)}"]
    lines.append("  " + "  ".join(f"{s}:{v}" for s, v in o["stages"].items()))
    c = o["current"]
    if c:
        lines.append(f"  running: {c['run']}  step {c['step']}/{c['nsteps']}"
                     + (f" ({100 * c['fraction']:.0f}%)" if c["fraction"] is not None else "")
                     + (f"  ETA {fmt_hms(c['eta_s'])}" if c["eta_s"] is not None else "")
                     + (f"  {c['ns_per_day']:.0f} ns/day" if c["ns_per_day"] else ""))
    if o["iterations_planned"] and o["iterations_done"]:
        lines.append(f"  iterations: {o['iterations_done']}/{o['iterations_planned']}"
                     + (f"  last RMSE(g) = {o['last_rmse']:.4f}" if o["last_rmse"] is not None else ""))
    if o["detail"]:
        lines.append("  " + o["detail"])
    if o["eta_s"] is not None:
        lines.append(f"  estimated time remaining (whole job): {fmt_hms(o['eta_s'])}")
    return "\n".join(lines)
