"""
Local execution of GROMACS jobs as detached subprocesses.

State lives on disk (``run.json`` plus a ``.exit_code`` file written by the
wrapper), so a run survives a restart of the application: status is derived
on every read from the exit-code file and the liveness of the process.
"""

from __future__ import annotations

import os
import re
import signal
import subprocess
from pathlib import Path
from typing import Optional

TERMINAL = {"done", "failed", "cancelled"}


def start(work: Path, log: Path, threads: int, gmx: str) -> int:
    """Launch ``run.sh`` in *work*; returns the PID of the process group leader."""
    work = Path(work)
    for stale in (work / ".exit_code",):
        if stale.exists():
            stale.unlink()
    env = dict(os.environ, GMX=gmx, NT=str(threads),
               OMP_NUM_THREADS=str(max(threads, 1)))
    with open(log, "ab") as fh:
        proc = subprocess.Popen(
            ["bash", "-c", "bash run.sh; echo $? > .exit_code"],
            cwd=work, stdout=fh, stderr=subprocess.STDOUT, env=env,
            start_new_session=True,
        )
    return proc.pid


def alive(pid: Optional[int]) -> bool:
    if not pid:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    # A zombie child of this process is not alive for our purposes.
    try:
        with open(f"/proc/{pid}/stat") as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return True


def cancel(pid: Optional[int]) -> None:
    if pid and alive(pid):
        try:
            os.killpg(os.getpgid(pid), signal.SIGTERM)
        except ProcessLookupError:
            pass


def resolve_status(record: dict, work: Path) -> dict:
    """Return *record* updated with the current status / exit code."""
    if record["status"] in TERMINAL:
        return record
    code_file = Path(work) / ".exit_code"
    if code_file.exists():
        try:
            code = int(code_file.read_text().strip())
        except ValueError:
            code = -1
        record["exit_code"] = code
        record["status"] = ("cancelled" if record.get("cancel_requested")
                            else "done" if code == 0 else "failed")
    elif not alive(record.get("pid")):
        record["status"] = "cancelled" if record.get("cancel_requested") else "failed"
        record["exit_code"] = record.get("exit_code")
        record["error"] = record.get("error") or "process ended without an exit code"
    else:
        record["status"] = "running"
    return record


_STEP = re.compile(r"Step\s+Time\s*\n\s*(\d+)\s+([\d.]+)")
_NSTEPS = re.compile(r"^\s*nsteps\s*=\s*(\d+)", re.M)
_DT = re.compile(r"^\s*dt\s*=\s*([\d.eE+-]+)", re.M)


def progress(work: Path) -> Optional[dict]:
    """Stage / step / fraction of the stage currently being written, if any."""
    logs = sorted(Path(work).glob("*.log"), key=lambda p: p.stat().st_mtime)
    if not logs:
        return None
    log = logs[-1]
    with open(log, "rb") as fh:
        fh.seek(0, 2)
        fh.seek(max(0, fh.tell() - 8192))
        tail = fh.read().decode(errors="replace")
    steps = _STEP.findall(tail)
    out = {"stage": log.stem, "step": None, "time_ps": None, "fraction": None}
    if steps:
        step, time = int(steps[-1][0]), float(steps[-1][1])
        out.update(step=step, time_ps=time)
        mdp = Path(work) / f"{log.stem}.mdp"
        if mdp.exists():
            txt = mdp.read_text()
            m = _NSTEPS.search(txt)
            if m and int(m.group(1)) > 0:
                out["fraction"] = min(1.0, step / int(m.group(1)))
    return out
