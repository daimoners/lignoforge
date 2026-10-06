"""
Boltzmann inversion of CG bond and angle distributions.

Atomistic frames are mapped onto the monomer beads (centre of mass) and the
distributions of bond lengths and angles are inverted,

    U(r)     = -kT ln[ P(r) / r^2 ]        (bonds)
    U(theta) = -kT ln[ P(theta) / sin(theta) ]   (angles)

then fitted with a harmonic well around the distribution minimum, which gives
``r0`` / ``theta0`` and ``k``.  Samples are pooled by linkage (bonds) and by
linkage pair (angles) over all chains supplied, so many atomistic chains of
different sequence contribute to the same parameter.
"""

from __future__ import annotations

import collections
from typing import Dict, List, Tuple

import numpy as np

from lignoforge.cg.mapping import CGChain, bead_positions
from lignoforge.cg.parameters import CGParameters, default_parameters, _angle_key

__all__ = ["BoltzmannFitter", "read_pdb_frames", "KB"]

KB = 0.0083144626  # kJ/mol/K

_MIN_SAMPLES = 200


def read_pdb_frames(path: str, n_atoms: int) -> np.ndarray:
    """
    Read a multi-model PDB (e.g. ``gmx trjconv -pbc mol -o frames.pdb``).

    Only the first *n_atoms* atoms of each model are kept (the solute comes
    first in a GROMACS system).  Returns ``(n_frames, n_atoms, 3)`` in nm.
    """
    frames: List[np.ndarray] = []
    cur: List[Tuple[float, float, float]] = []
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM"):
                if len(cur) < n_atoms:
                    cur.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            elif rec in ("ENDMDL", "END") and cur:
                if len(cur) != n_atoms:
                    raise ValueError(f"Model has {len(cur)} atoms, expected {n_atoms}")
                frames.append(np.array(cur) / 10.0)
                cur = []
    if cur and len(cur) == n_atoms:
        frames.append(np.array(cur) / 10.0)
    if not frames:
        raise ValueError(f"No frames read from {path}")
    return np.array(frames)


def _harmonic_from_histogram(x, weights_fn, kT, lo, hi, bins=60, bounds=None):
    """
    Fit U(x)=0.5 k (x-x0)^2 to -kT ln(P/jacobian) around its minimum.

    ``bounds=(xmin, xmax)`` is the physical domain of x.  If the unconstrained
    vertex x0 falls outside it (typical for angles near 180 degrees, where the
    distribution is truncated), the fit is repeated with x0 fixed at the
    nearest bound.
    """
    hist, edges = np.histogram(x, bins=bins, range=(lo, hi), density=True)
    c = 0.5 * (edges[1:] + edges[:-1])
    jac = weights_fn(c)
    ok = hist > 0
    with np.errstate(divide="ignore", invalid="ignore"):
        U = -kT * np.log(hist / jac)
    U[~ok] = np.nan
    # restrict to the well: contiguous bins within 2.5 kT of the minimum
    imin = int(np.nanargmin(U))
    sel = np.zeros_like(ok)
    for step in (-1, 1):
        i = imin
        while 0 <= i < len(U) and ok[i] and U[i] - U[imin] < 2.5 * kT:
            sel[i] = True
            i += step
    if sel.sum() < 4:
        return None
    w = np.sqrt(hist[sel])
    a, b, _ = np.polyfit(c[sel], U[sel], 2, w=w)
    if a <= 0:
        return None
    x0 = -b / (2 * a)
    if bounds is not None and not (bounds[0] <= x0 <= bounds[1]):
        x0 = min(max(x0, bounds[0]), bounds[1])
        slope, _ = np.polyfit((c[sel] - x0) ** 2, U[sel], 1, w=w)
        if slope <= 0:
            return None
        return x0, 2 * slope
    return x0, 2 * a


class BoltzmannFitter:
    """Accumulate atomistic frames of any number of chains, then fit."""

    def __init__(self) -> None:
        self._bonds: Dict[str, List[float]] = collections.defaultdict(list)
        self._angles: Dict[str, List[float]] = collections.defaultdict(list)

    def add(self, cg: CGChain, atom_frames_nm: np.ndarray) -> int:
        """Add frames ``(n_frames, n_atoms, 3)`` of one chain; returns n_frames."""
        for frame in atom_frames_nm:
            x = bead_positions(cg.beads, frame)
            for i, j, link in cg.bonds:
                self._bonds[link].append(float(np.linalg.norm(x[i - 1] - x[j - 1])))
            for i, j, k, l1, l2 in cg.angles:
                a, b = x[i - 1] - x[j - 1], x[k - 1] - x[j - 1]
                cos = a @ b / (np.linalg.norm(a) * np.linalg.norm(b))
                self._angles[_angle_key(l1, l2)].append(
                    float(np.degrees(np.arccos(np.clip(cos, -1, 1)))))
        return len(atom_frames_nm)

    def fit(
        self,
        temperature: float = 300.0,
        base: CGParameters | None = None,
        min_samples: int = _MIN_SAMPLES,
    ) -> Tuple[CGParameters, dict]:
        """
        Return ``(parameters, report)``.  Entries with fewer than
        ``min_samples`` samples, or whose well cannot be fitted, keep the
        values of *base* (default: provisional defaults) and are listed in
        ``report["kept_default"]``.
        """
        kT = KB * temperature
        params = base or default_parameters()
        params = CGParameters.from_dict(params.to_dict())
        report = {"temperature": temperature, "fitted": {}, "kept_default": []}

        for link, r in self._bonds.items():
            r = np.asarray(r)
            res = None
            if len(r) >= min_samples:
                res = _harmonic_from_histogram(
                    r, lambda c: c ** 2, kT, r.min(), r.max())
            if res is None:
                report["kept_default"].append(f"bond:{link} (n={len(r)})")
                continue
            r0, k = res
            params.bonds[link] = {"r0": round(float(r0), 4), "k": round(float(k), 1)}
            report["fitted"][f"bond:{link}"] = {**params.bonds[link], "n": len(r)}

        for key, th in self._angles.items():
            th = np.asarray(th)
            res = None
            if len(th) >= min_samples:
                res = _harmonic_from_histogram(
                    np.radians(th), np.sin, kT, 0.2, np.pi,
                    bounds=(0.0, np.pi))
            if res is None:
                report["kept_default"].append(f"angle:{key} (n={len(th)})")
                continue
            t0, k = res
            params.angles[key] = {"theta0": round(float(np.degrees(t0)), 2),
                                  "k": round(float(k), 2)}
            report["fitted"][f"angle:{key}"] = {**params.angles[key], "n": len(th)}

        params.provenance = (
            f"Boltzmann inversion of atomistic MD at {temperature} K; entries "
            f"not in report['fitted'] keep provisional defaults."
        )
        return params, report
