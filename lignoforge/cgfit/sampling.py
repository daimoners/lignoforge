"""
Statistics of a bead system from trajectories: bonded samples, non-excluded
pair-distance histograms (per bead-type pair, per frame) and chain radii of
gyration.  The same code runs on the mapped atomistic reference and on CG
trajectories, so the two are directly comparable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

from lignoforge.cgfit.model import CGMolecule, all_type_pairs, pair_key

__all__ = ["SystemMap", "Stats", "iter_frames", "sample"]


class SystemMap:
    """Global (whole-system) index arrays for a list of molecules with copies."""

    def __init__(self, mols: Sequence[CGMolecule], copies: Sequence[int]):
        self.mols, self.copies = list(mols), list(copies)
        atom_bead, atom_mass, btypes, bmass = [], [], [], []
        bonds, angles, excl, mol_of_bead, mol_slices = [], [], [], [], []
        boff = aoff = 0
        for m, c in zip(mols, copies):
            for _ in range(c):
                atom_bead.append(m.atom_bead + boff)
                atom_mass.append(m.atom_mass)
                btypes += m.types
                bmass.append(m.masses)
                bonds += [(i - 1 + boff, j - 1 + boff, k) for i, j, k in m.bonds]
                angles += [(i - 1 + boff, j - 1 + boff, k - 1 + boff, s, g)
                           for i, j, k, s, g in m.angles]
                excl += [(i - 1 + boff, j - 1 + boff) for i, j in m.exclusions]
                mol_slices.append((boff, boff + len(m.beads), m.name))
                boff += len(m.beads)
                aoff += m.n_atoms
        self.n_atoms, self.n_beads = aoff, boff
        self.atom_bead = np.concatenate(atom_bead)
        self.atom_mass = np.concatenate(atom_mass)
        self.bead_mass = np.concatenate(bmass)
        self.types = sorted(set(btypes))
        self.bead_type = np.array([self.types.index(t) for t in btypes])
        self.mol_slices = mol_slices

        self.bond_ij = np.array([(i, j) for i, j, _ in bonds], dtype=int).reshape(-1, 2)
        self.bond_keys = sorted({k for *_, k in bonds})
        kid = {k: n for n, k in enumerate(self.bond_keys)}
        self.bond_kid = np.array([kid[k] for *_, k in bonds], dtype=int)
        self.ang_ijk = np.array([(i, j, k) for i, j, k, *_ in angles], dtype=int).reshape(-1, 3)
        self.ang_keys = sorted({s for *_, s, g in angles} | {g for *_, s, g in angles})
        akid = {k: n for n, k in enumerate(self.ang_keys)}
        self.ang_spec = np.array([akid[s] for *_, s, g in angles], dtype=int)
        self.ang_gen = np.array([akid[g] for *_, s, g in angles], dtype=int)

        self.pair_names = all_type_pairs(self.types)
        pid = {p: n for n, p in enumerate(self.pair_names)}
        nt = len(self.types)
        self.pair_index = np.zeros((nt, nt), dtype=int)
        for a in range(nt):
            for b in range(nt):
                self.pair_index[a, b] = pid[pair_key(self.types[a], self.types[b])]
        n = self.n_beads
        e = np.array([(min(i, j), max(i, j)) for i, j in excl], dtype=np.int64).reshape(-1, 2)
        self.excl_codes = np.unique(e[:, 0] * n + e[:, 1])
        self.type_counts = np.bincount(self.bead_type, minlength=nt)

    def map_positions(self, atom_xyz_nm: np.ndarray) -> np.ndarray:
        """Centre-of-mass bead positions for one atomistic frame (molecules whole)."""
        w = self.atom_mass[:, None]
        num = np.zeros((self.n_beads, 3))
        np.add.at(num, self.atom_bead, atom_xyz_nm * w)
        return num / self.bead_mass[:, None]


@dataclass
class Stats:
    pair_names: List[str]
    r_edges: np.ndarray
    pair_hist: List[np.ndarray] = field(default_factory=list)      # per frame (n_pairs, n_bins)
    volumes: List[float] = field(default_factory=list)
    bonds: Dict[str, List[np.ndarray]] = field(default_factory=dict)
    angles: Dict[str, List[np.ndarray]] = field(default_factory=dict)   # degrees
    rg: Dict[str, List[float]] = field(default_factory=dict)
    type_counts: Optional[np.ndarray] = None
    types: Optional[List[str]] = None
    n_frames: int = 0

    @property
    def r_centers(self) -> np.ndarray:
        return 0.5 * (self.r_edges[1:] + self.r_edges[:-1])

    def mean_hist(self) -> np.ndarray:
        return np.mean(self.pair_hist, axis=0)

    def gr(self, sm_pairs_n: Dict[str, float]) -> Dict[str, np.ndarray]:
        """g(r) per pair: counts / expected counts of an ideal gas at the mean density."""
        V = float(np.mean(self.volumes))
        shell = 4 * np.pi / 3 * (self.r_edges[1:] ** 3 - self.r_edges[:-1] ** 3)
        h = self.mean_hist()
        return {p: h[n] / (sm_pairs_n[p] / V * shell) for n, p in enumerate(self.pair_names)}


def ideal_pair_numbers(sm: SystemMap) -> Dict[str, float]:
    """Number of unordered bead pairs per type pair (ideal-gas normalisation)."""
    out = {}
    for a, ta in enumerate(sm.types):
        for b in range(a, len(sm.types)):
            na, nb = sm.type_counts[a], sm.type_counts[b]
            out[pair_key(ta, sm.types[b])] = na * (na - 1) / 2 if a == b else na * nb
    return out


def iter_frames(gro: str, xtc: str, begin_ps: float = 0.0, stride: int = 1
                ) -> Iterator[Tuple[np.ndarray, np.ndarray]]:
    """Yield ``(xyz_nm, box_nm)`` for each frame at or after *begin_ps*."""
    import MDAnalysis as mda
    u = mda.Universe(gro, xtc)
    for n, ts in enumerate(u.trajectory):
        if ts.time < begin_ps - 1e-6 or n % stride:
            continue
        yield ts.positions.astype(float) / 10.0, ts.dimensions[:3].astype(float) / 10.0


def _mic(d: np.ndarray, box: np.ndarray) -> np.ndarray:
    return d - box * np.round(d / box)


def sample(sm: SystemMap, frames: Iterator[Tuple[np.ndarray, np.ndarray]], mapped: bool,
           rc: float, dr: float = 0.01) -> Stats:
    """
    Accumulate statistics.  ``mapped=False``: frames are atomistic and are
    projected on the beads; ``True``: frames already hold bead positions.
    """
    edges = np.arange(0.0, rc + 1e-9, dr)
    nb = len(edges) - 1
    st = Stats(sm.pair_names, edges, type_counts=sm.type_counts, types=sm.types)
    n = sm.n_beads
    bkeys = [np.where(sm.bond_kid == k)[0] for k in range(len(sm.bond_keys))]
    for xyz, box in frames:
        x = xyz if mapped else sm.map_positions(xyz)
        if len(x) != n:
            raise ValueError(f"frame has {len(x)} particles, expected {n}")
        st.n_frames += 1
        st.volumes.append(float(np.prod(box)))

        if len(sm.bond_ij):
            d = np.linalg.norm(_mic(x[sm.bond_ij[:, 0]] - x[sm.bond_ij[:, 1]], box), axis=1)
            for k, idx in enumerate(bkeys):
                st.bonds.setdefault(sm.bond_keys[k], []).append(d[idx])
        if len(sm.ang_ijk):
            a = _mic(x[sm.ang_ijk[:, 0]] - x[sm.ang_ijk[:, 1]], box)
            b = _mic(x[sm.ang_ijk[:, 2]] - x[sm.ang_ijk[:, 1]], box)
            cos = (a * b).sum(1) / np.linalg.norm(a, axis=1) / np.linalg.norm(b, axis=1)
            th = np.degrees(np.arccos(np.clip(cos, -1, 1)))
            for ids in (sm.ang_spec, sm.ang_gen):
                for k in np.unique(ids):
                    st.angles.setdefault(sm.ang_keys[k], []).append(th[ids == k])

        tree = cKDTree(np.mod(x, box), boxsize=box)
        pr = tree.query_pairs(rc, output_type="ndarray")
        if len(pr):
            i, j = pr[:, 0].astype(np.int64), pr[:, 1].astype(np.int64)
            lo, hi = np.minimum(i, j), np.maximum(i, j)
            keep = ~np.isin(lo * n + hi, sm.excl_codes, assume_unique=False)
            i, j = i[keep], j[keep]
            dist = np.linalg.norm(_mic(x[i] - x[j], box), axis=1)
            cls = sm.pair_index[sm.bead_type[i], sm.bead_type[j]]
            b = np.minimum((dist / dr).astype(int), nb - 1)
            h = np.bincount(cls * nb + b, minlength=len(sm.pair_names) * nb)
            st.pair_hist.append(h.reshape(len(sm.pair_names), nb).astype(np.float64))
        else:
            st.pair_hist.append(np.zeros((len(sm.pair_names), nb)))

        for s, e, name in sm.mol_slices:
            xm = x[s:e]
            xm = xm[0] + _mic(xm - xm[0], box)
            m = sm.bead_mass[s:e]
            com = (xm * m[:, None]).sum(0) / m.sum()
            st.rg.setdefault(name, []).append(
                float(np.sqrt((m * ((xm - com) ** 2).sum(1)).sum() / m.sum())))
    if not st.n_frames:
        raise ValueError("no frames to analyse")
    return st
