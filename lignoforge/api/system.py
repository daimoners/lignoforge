"""Environment check and model-validation status reported to users."""

from __future__ import annotations

import platform
import re
import shutil
import subprocess
from typing import Optional

from lignoforge.api.store import lignoforge_version

# Honest, user-visible status of every scientific component.  The GUI shows
# these so that provisional models are never mistaken for validated ones.
MODEL_STATUS = {
    "structure_generation": {
        "status": "validated-internally",
        "detail": "Chemistry sanity-checked (RDKit); linkage/composition statistics "
                  "not yet compared with experimental NMR/GPC data."},
    "opls_types": {
        "status": "validated-internally",
        "detail": "All bonds/angles/torsions of generated chains resolve in OPLS-AA "
                  "(grompp); two torsions are analogues (see docs)."},
    "opls_charges": {
        "status": "provisional",
        "detail": "OPLS-AA fragment charges with per-residue neutralisation "
                  "(shifts ≤ ~0.013 e); not compared with RESP/DFT."},
    "cg_model": {
        "status": "provisional",
        "detail": "Bond r0 from MMFF ensembles; force constants, bead sizes and "
                  "well depths are placeholders. Fit from atomistic MD before use."},
}


def _version(cmd: list) -> Optional[str]:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        m = re.search(r"GROMACS version:\s*(\S+)", out.stdout + out.stderr)
        return m.group(1) if m else None
    except Exception:                       # noqa: BLE001
        return None


def system_check(gmx: str = "gmx") -> dict:
    """Report installed dependencies and what features are therefore available."""
    def mod(name):
        try:
            m = __import__(name)
            return getattr(m, "__version__", "installed")
        except Exception:                   # noqa: BLE001
            return None

    gmx_path = shutil.which(gmx)
    info = {
        "lignoforge": lignoforge_version(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "gromacs": {"path": gmx_path, "version": _version([gmx, "--version"]) if gmx_path else None},
        "rdkit": mod("rdkit"),
        "mdanalysis": mod("MDAnalysis"),
        "networkx": mod("networkx"),
    }
    gmx_ok = bool(info["gromacs"]["version"])
    info["features"] = {
        "build_structures": bool(info["rdkit"]),
        "topology": bool(info["rdkit"]),
        "simulations": gmx_ok,
        "trajectory_analysis": bool(info["mdanalysis"]),
        "energy_analysis": gmx_ok,
    }
    info["hints"] = [
        h for h in (
            None if info["rdkit"] else "RDKit missing: conda install -c conda-forge rdkit",
            None if gmx_ok else "GROMACS not found on PATH: conda install -c conda-forge gromacs",
            None if info["mdanalysis"] else "MDAnalysis missing: pip install 'lignoforge[analysis]'",
        ) if h
    ]
    info["models"] = MODEL_STATUS
    return info
