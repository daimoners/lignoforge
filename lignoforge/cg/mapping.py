"""
Coarse-grained mapping: one bead per monomer, placed at its centre of mass.

The CG graph is exactly the monomer graph (``bigG``): beads are monomers,
CG bonds are inter-monomer linkages, CG angles are pairs of bonds sharing a
bead.  Mapping is defined on a :class:`~lignoforge.forcefield.ChainTopology`
so that atom indices in an atomistic MD trajectory can be mapped onto the same
beads (see :mod:`lignoforge.cg.boltzmann`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

from lignoforge.forcefield.topology import ChainTopology

__all__ = ["CGBead", "CGChain", "map_topology", "bead_positions"]

_MASS = {"C": 12.011, "H": 1.008, "O": 15.999}
_RES_TO_TYPE = {"HPU": "H", "GYU": "G", "SYU": "S"}


@dataclass
class CGBead:
    index: int                 # 1-based
    mtype: str                 # 'H' | 'G' | 'S'
    resname: str
    mass: float                # g/mol, sum of the atoms of the monomer
    atom_indices: List[int]    # 0-based indices into ChainTopology.atoms
    elements: List[str]        # element of each of those atoms


@dataclass
class CGChain:
    name: str
    beads: List[CGBead]
    bonds: List[Tuple[int, int, str]]                 # (i, j, linkage), i<j, 1-based
    angles: List[Tuple[int, int, int, str, str]]      # (i, j, k, link_ij, link_jk)
    ref_xyz: np.ndarray                               # (n_beads, 3) nm

    def mass_array(self) -> np.ndarray:
        return np.array([b.mass for b in self.beads])


def bead_positions(beads: Sequence[CGBead], atom_xyz_nm: np.ndarray) -> np.ndarray:
    """Centre-of-mass positions (nm) of every bead for one frame ``(n_atoms, 3)``."""
    out = np.empty((len(beads), 3))
    for b, bead in enumerate(beads):
        m = np.array([_MASS[a] for a in bead.elements])
        out[b] = (atom_xyz_nm[bead.atom_indices] * m[:, None]).sum(0) / m.sum()
    return out


def map_topology(topo: ChainTopology) -> CGChain:
    """Build the CG chain (beads, bonds, angles, reference positions)."""
    by_res: Dict[int, List[int]] = {}
    for i, a in enumerate(topo.atoms):
        by_res.setdefault(a.resid, []).append(i)

    beads: List[CGBead] = []
    for resid in sorted(by_res):
        idx = by_res[resid]
        resname = topo.atoms[idx[0]].resname
        if resname not in _RES_TO_TYPE:
            raise ValueError(f"Unknown residue name '{resname}' for CG mapping")
        bead = CGBead(
            index=resid, mtype=_RES_TO_TYPE[resname], resname=resname,
            mass=sum(_MASS[topo.atoms[i].element] for i in idx), atom_indices=idx,
            elements=[topo.atoms[i].element for i in idx],
        )
        beads.append(bead)

    seen: Dict[Tuple[int, int], str] = {}
    for lk in topo.linkages:
        key = tuple(sorted((lk["residue_1"], lk["residue_2"])))
        seen.setdefault(key, lk["linkage"])
    bonds = [(i, j, l) for (i, j), l in sorted(seen.items())]

    nbrs: Dict[int, Dict[int, str]] = {}
    for i, j, l in bonds:
        nbrs.setdefault(i, {})[j] = l
        nbrs.setdefault(j, {})[i] = l
    angles = []
    for j in sorted(nbrs):
        ks = sorted(nbrs[j])
        for a in range(len(ks)):
            for b in range(a + 1, len(ks)):
                angles.append((ks[a], j, ks[b], nbrs[j][ks[a]], nbrs[j][ks[b]]))

    xyz_nm = np.array([a.xyz for a in topo.atoms]) / 10.0
    return CGChain(
        name=topo.name, beads=beads, bonds=bonds, angles=angles,
        ref_xyz=bead_positions(beads, xyz_nm),
    )
