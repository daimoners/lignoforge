"""
Machine-aware configuration for the CG parametriser: detects GROMACS (and its GPU
support) and writes a complete, consistent ``config.json`` from a few options.
"""

from __future__ import annotations

import copy
import os
import re
import shutil
import subprocess
from typing import Dict, Optional

from lignoforge.cgfit.gmxio import box_edge_for_density
from lignoforge.cgfit.workflow import DEFAULT_CONFIG, _merge, quick_config

SIZES = {"small": (4, 4), "medium": (4, 6), "large": (4, 10)}     # (distinct sequences, copies)
_ATOMS_PER_MONOMER, _DALTON_PER_ATOM = 25.5, 7.3                    # measured on built chains


def detect_gmx(gmx: Optional[str] = None) -> dict:
    """``{"path", "version", "gpu", "gpu_kind"}``; ``path`` is None if GROMACS is missing."""
    path = gmx or os.environ.get("GMX") or shutil.which("gmx")
    if not path or not (shutil.which(path) or os.path.exists(path)):
        return {"path": None, "version": None, "gpu": False, "gpu_kind": None}
    try:
        out = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return {"path": path, "version": None, "gpu": False, "gpu_kind": None}
    ver = re.search(r"GROMACS version:\s*(\S+)", out)
    gpu = re.search(r"GPU support:\s*(\S+)", out)
    kind = gpu.group(1) if gpu else None
    return {"path": path, "version": ver.group(1) if ver else None,
            "gpu": bool(kind) and kind.lower() not in ("disabled", "none"), "gpu_kind": kind}


def _mdrun_args(gpu: bool, threads: int):
    atom = ["-ntmpi", "1", "-ntomp", str(threads)]
    if gpu:
        atom += ["-nb", "gpu", "-pme", "gpu", "-update", "gpu"]
    return atom, ["-ntmpi", "1", "-ntomp", str(threads)]


def make_config(*, quick: bool = False, scheme: Optional[str] = None,
                temperature: Optional[float] = None, mixing: Optional[str] = None,
                size: Optional[str] = None, n_distinct: Optional[int] = None,
                copies: Optional[int] = None, monomers: Optional[int] = None,
                fractions: Optional[Dict[str, float]] = None, seed: Optional[int] = None,
                gpu: str = "auto", threads: Optional[int] = None, gmx: Optional[str] = None,
                compress_ns: Optional[float] = None, equil_ns: Optional[float] = None,
                prod_ns: Optional[float] = None, iterations: Optional[int] = None,
                iter_ns: Optional[float] = None, cutoff: Optional[float] = None) -> tuple:
    """Return ``(config, info)``; *info* has the machine detection and the system estimate."""
    cfg = copy.deepcopy(quick_config() if quick else DEFAULT_CONFIG)
    if scheme:
        cfg["scheme"] = scheme
    if temperature:
        cfg["temperature"] = float(temperature)
    if mixing:
        cfg["cg"]["mixing"] = mixing
    if seed is not None:
        cfg["seed"] = seed
        cfg["ensemble"]["spec"]["seed"] = seed
    e = cfg["ensemble"]
    if size and not quick:
        e["n_distinct"], e["copies"] = SIZES[size]
    if n_distinct:
        e["n_distinct"] = n_distinct
    if copies:
        e["copies"] = copies
    if monomers:
        e["spec"]["n_monomers"] = monomers
    if fractions:
        e["spec"]["monomers"] = {k: float(fractions.get(k, 0.0)) for k in "HGS"}
    a, c = cfg["atomistic"], cfg["cg"]
    if compress_ns is not None:
        a["compress_ps"] = compress_ns * 1000
    if equil_ns is not None:
        a["equil_ps"] = equil_ns * 1000
    if prod_ns is not None:
        a["prod_ps"] = prod_ns * 1000
    if iterations:
        c["iterations"] = iterations
    if iter_ns is not None:
        c["iter_ps"] = iter_ns * 1000
    if cutoff:
        c["cutoff"] = cutoff

    g = detect_gmx(gmx)
    if g["path"]:
        cfg["gmx"] = os.path.abspath(shutil.which(g["path"]) or g["path"])
    use_gpu = g["gpu"] if gpu == "auto" else gpu == "yes"
    n_thr = threads or os.cpu_count() or 4
    if g["path"] or gpu != "auto":
        a["mdrun_args"], c["mdrun_args"] = _mdrun_args(use_gpu, n_thr)

    n_chains = int(e["n_distinct"]) * int(e["copies"]) if not e.get("chains") else \
        len(e["chains"]) * int(e["copies"])
    atoms = int(n_chains * e["spec"]["n_monomers"] * _ATOMS_PER_MONOMER)
    edge = box_edge_for_density(atoms * _DALTON_PER_ATOM, 1.1)
    warnings = []
    if edge < 2.3:
        warnings.append(f"box ≈ {edge:.2f} nm at melt density is too small for the atomistic cut-off: more chains needed")
    elif edge < 2.2 * c["cutoff"]:
        warnings.append(f"box ≈ {edge:.2f} nm is too small for the CG cut-off {c['cutoff']} nm "
                        f"(needs > {2.2 * c['cutoff']:.2f}): more chains or a shorter cut-off")
    info = {"gmx": g, "gpu_used": use_gpu, "threads": n_thr, "chains": n_chains,
            "atoms_estimate": atoms, "box_nm_estimate": round(edge, 2), "warnings": warnings}
    return cfg, info
