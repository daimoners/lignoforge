"""
Trajectory observables.  Every function returns a JSON-serialisable dict::

    {"observable": name, "time_ps": [...], "series": {label: [...]},
     "summary": {label: {"mean", "std", ...}}, "extra": {...}}
"""

from __future__ import annotations

import re
import subprocess
import tempfile
from typing import Dict, List, Optional

import numpy as np

from lignoforge.analysis.core import (
    AMU_A3_TO_G_CM3, AnalysisError, Selection, graph_endpoints,
)

__all__ = ["radius_of_gyration", "end_to_end", "rdf", "density", "contacts",
           "rmsd", "energy", "OBSERVABLES"]


def _summ(x) -> dict:
    x = np.asarray(x, float)
    n = len(x)
    half = x[n // 2:] if n >= 4 else x
    return {"mean": float(x.mean()), "std": float(x.std()),
            "min": float(x.min()), "max": float(x.max()),
            "mean_second_half": float(half.mean())}


def _frames(u, start=0, stop=None, step=1):
    return u.trajectory[start:stop:step]


def _result(name, time, series, **extra) -> dict:
    return {
        "observable": name,
        "time_ps": [float(t) for t in time],
        "series": {k: [float(v) for v in vals] for k, vals in series.items()},
        "summary": {k: _summ(v) for k, v in series.items() if len(v)},
        "extra": extra,
    }


def _com(ag, groups) -> np.ndarray:
    """Centre of mass of every monomer, shape (n_monomers, 3)."""
    m = ag.masses
    pos = ag.positions
    return np.array([(pos[g] * m[g, None]).sum(0) / m[g].sum() for g in groups])


# ── observables ───────────────────────────────────────────────────────────────

def radius_of_gyration(u, sel: Selection, start=0, stop=None, step=1, **_) -> dict:
    """Mass-weighted Rg and asphericity of the lignin molecule."""
    t, rg, asph = [], [], []
    for ts in _frames(u, start, stop, step):
        t.append(ts.time)
        rg.append(sel.atoms.radius_of_gyration())
        asph.append(sel.atoms.asphericity())
    return _result("rg", t, {"Rg (Å)": rg, "asphericity": asph})


def end_to_end(u, sel: Selection, chain_data: dict, start=0, stop=None, step=1, **_) -> dict:
    """Distance between the centres of mass of the two most distant monomers."""
    a, b = graph_endpoints(chain_data)
    t, d = [], []
    for ts in _frames(u, start, stop, step):
        com = _com(sel.atoms, sel.monomers)
        t.append(ts.time)
        d.append(float(np.linalg.norm(com[a - 1] - com[b - 1])))
    return _result("end_to_end", t, {"End-to-end (Å)": d}, monomers=[a, b])


def density(u, start=0, stop=None, step=1, **_) -> dict:
    """Mass density of the whole simulation box (g/cm³)."""
    m = float(u.atoms.masses.sum())
    t, rho, vol = [], [], []
    for ts in _frames(u, start, stop, step):
        v = float(ts.volume)
        if v <= 0:
            raise AnalysisError("Trajectory has no box (volume 0)")
        t.append(ts.time)
        rho.append(m / v * AMU_A3_TO_G_CM3)
        vol.append(v / 1000.0)
    return _result("density", t, {"Density (g/cm³)": rho, "Volume (nm³)": vol})


def rdf(u, sel: Selection, group_a: str = "solute", group_b: str = "water_O",
        rmax: float = 15.0, nbins: int = 150, start=0, stop=None, step=1, **_) -> dict:
    """
    Radial distribution function between two groups.

    ``group_a`` / ``group_b`` are ``"solute"``, ``"solute_heavy"``,
    ``"solute_O"``, ``"water_O"`` or any MDAnalysis selection string.
    """
    from MDAnalysis.analysis.rdf import InterRDF

    def pick(name):
        solute_idx = sel.atoms.indices
        if name == "solute":
            return sel.atoms
        if name == "solute_heavy":
            return sel.atoms[[i for i, a in enumerate(sel.atoms) if not a.name.startswith("H")]]
        if name == "solute_O":
            return sel.atoms[[i for i, a in enumerate(sel.atoms) if a.name.startswith("O")]]
        if name == "water_O":
            g = u.select_atoms("resname SOL and name OW")
            if not g.n_atoms:
                raise AnalysisError("No water oxygens (resname SOL, name OW) in this run")
            return g
        g = u.select_atoms(name)
        if not g.n_atoms:
            raise AnalysisError(f"Selection '{name}' matched no atoms")
        return g

    ga, gb = pick(group_a), pick(group_b)
    if np.all(u.dimensions[:3] <= 0):
        raise AnalysisError("RDF needs a periodic box")
    exclusion = (1, 1) if group_a == group_b else None
    r = InterRDF(ga, gb, nbins=nbins, range=(0.0, rmax), exclusion_block=exclusion)
    r.run(start=start, stop=stop, step=step)
    return {
        "observable": "rdf", "time_ps": [], "series": {}, "summary": {},
        "extra": {"r": [float(x) for x in r.results.bins],
                  "g": [float(x) for x in r.results.rdf],
                  "group_a": group_a, "group_b": group_b,
                  "n_a": int(ga.n_atoms), "n_b": int(gb.n_atoms)},
    }


def contacts(u, sel: Selection, cutoff: float = 4.5, start=0, stop=None, step=1, **_) -> dict:
    """
    Inter-monomer contact map: fraction of frames in which the minimum
    heavy-atom distance between two monomers is below *cutoff* (Å).
    """
    from MDAnalysis.lib.distances import distance_array

    heavy = [i for i, a in enumerate(sel.atoms) if not a.name.startswith("H")]
    owner = np.full(len(sel.atoms), -1)
    for k, g in enumerate(sel.monomers):
        owner[g] = k
    heavy = np.array(heavy)
    order = np.argsort(owner[heavy], kind="stable")
    heavy = heavy[order]
    res_of = owner[heavy]
    bounds = np.r_[0, np.flatnonzero(np.diff(res_of)) + 1]
    n = len(sel.monomers)

    acc = np.zeros((n, n))
    mind = np.zeros((n, n))
    nfr = 0
    for ts in _frames(u, start, stop, step):
        pos = sel.atoms.positions[heavy]
        D = distance_array(pos, pos)
        M = np.minimum.reduceat(np.minimum.reduceat(D, bounds, axis=0), bounds, axis=1)
        acc += (M < cutoff)
        mind += M
        nfr += 1
    if not nfr:
        raise AnalysisError("No frames selected")
    frac = acc / nfr
    np.fill_diagonal(frac, 0.0)
    return {
        "observable": "contacts", "time_ps": [], "series": {}, "summary": {},
        "extra": {"matrix": frac.tolist(), "mean_min_distance": (mind / nfr).tolist(),
                  "cutoff": cutoff, "n_frames": nfr,
                  "labels": [f"{k + 1}" for k in range(n)]},
    }


def rmsd(u, sel: Selection, start=0, stop=None, step=1, **_) -> dict:
    """Heavy-atom RMSD of the molecule to its first analysed frame (Å), after fitting."""
    from MDAnalysis.analysis import rms

    heavy = sel.atoms[[i for i, a in enumerate(sel.atoms) if not a.name.startswith("H")]]
    u.trajectory[start]
    r = rms.RMSD(heavy, heavy, ref_frame=start)
    r.run(start=start, stop=stop, step=step)
    res = r.results.rmsd
    return _result("rmsd", res[:, 1], {"RMSD (Å)": res[:, 2]})


def energy(edr: str, terms: Optional[List[str]] = None, gmx: str = "gmx", **_) -> dict:
    """
    Energy-file terms (e.g. Potential, Temperature, Pressure, Density) read
    with ``gmx energy``.  Terms absent from the file are reported in
    ``extra['missing']``.
    """
    terms = terms or ["Potential", "Kinetic-En.", "Total-Energy", "Temperature",
                      "Pressure", "Density"]
    with tempfile.TemporaryDirectory() as tmp:
        menu = subprocess.run([gmx, "energy", "-f", edr, "-o", f"{tmp}/q.xvg"],
                              input="0\n", capture_output=True, text=True, cwd=tmp)
        text = menu.stdout + menu.stderr
        avail = {name: int(num) for num, name in
                 re.findall(r"(\d+)\s+([A-Za-z][\w.\-()+]*)", text.split("Select the terms")[-1])}
        wanted = [t for t in terms if t in avail]
        missing = [t for t in terms if t not in avail]
        if not wanted:
            raise AnalysisError(f"None of {terms} found in {edr}; available: {sorted(avail)}")
        stdin = "\n".join(str(avail[t]) for t in wanted) + "\n\n"
        r = subprocess.run([gmx, "energy", "-f", edr, "-o", f"{tmp}/e.xvg"],
                           input=stdin, capture_output=True, text=True, cwd=tmp)
        try:
            data = np.loadtxt(f"{tmp}/e.xvg", comments=("#", "@"), ndmin=2)
        except Exception as e:           # noqa: BLE001
            raise AnalysisError(f"gmx energy failed: {r.stderr[-300:]}") from e
    series = {t: data[:, i + 1] for i, t in enumerate(wanted)}
    return _result("energy", data[:, 0], series, missing=missing)


OBSERVABLES = {
    "rg": radius_of_gyration, "end_to_end": end_to_end, "density": density,
    "rdf": rdf, "contacts": contacts, "rmsd": rmsd, "energy": energy,
}
