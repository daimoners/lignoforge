"""
Coarse-grained model definition for the parametriser.

* **Mapping schemes** assign every atom of a chain to a bead:

  ``monomer``    one bead per monomer (types ``MH``, ``MG``, ``MS``);
  ``ring-side``  two beads per monomer: the aromatic ring with its methoxyls
                 and phenolic oxygen (``RH``/``RG``/``RS``) and the propyl side
                 chain with its oxygens (``SC``).

* :class:`CGMolecule` is the bead graph of one chain: beads, bonds (labelled by
  linkage class), angles and 1-2/1-3 exclusions, plus the atom→bead map used
  to project atomistic trajectories onto the beads.
* :class:`CGParams` holds everything that is fitted: Lennard-Jones ``sigma`` and
  ``epsilon`` per bead-type pair, harmonic bonds and angles.  It is a plain,
  JSON-serialisable object.
"""

from __future__ import annotations

import itertools
import json
import re
from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

import numpy as np

from lignoforge.forcefield.topology import ChainTopology

__all__ = ["SCHEMES", "CGMolecule", "CGParams", "map_chain", "pair_key", "MASS"]

MASS = {"C": 12.011, "H": 1.008, "O": 15.999}
_RES_LETTER = {"HPU": "H", "GYU": "G", "SYU": "S"}
_SIDE = re.compile(r"^(C|O|H|HO)[ABG]\d*$")          # CA CB CG OA OB OG HA HB1 HG2 HOG ...

SCHEMES: Dict[str, str] = {
    "monomer": "one bead per monomer (centre of mass)",
    "ring-side": "two beads per monomer: aromatic ring + OMe/OH, and propyl side chain",
}


def pair_key(a: str, b: str) -> str:
    """Canonical key of an unordered bead-type pair."""
    return "-".join(sorted((a, b)))


def _bead_of_atom(scheme: str, name: str, resname: str) -> Tuple[str, str]:
    """Return ``(local_key, bead_type)`` of one atom."""
    letter = _RES_LETTER.get(resname)
    if letter is None:
        raise ValueError(f"Unknown residue name '{resname}' for CG mapping")
    if scheme == "monomer":
        return "M", f"M{letter}"
    if scheme == "ring-side":
        return ("S", "SC") if _SIDE.match(name) else ("R", f"R{letter}")
    raise ValueError(f"Unknown CG scheme '{scheme}' (available: {sorted(SCHEMES)})")


@dataclass
class Bead:
    index: int                 # 1-based within the molecule
    type: str
    resid: int
    key: str                   # local key within the residue ('M', 'R', 'S')
    mass: float
    atoms: List[int]           # 0-based atom indices within the molecule


@dataclass
class CGMolecule:
    name: str
    scheme: str
    beads: List[Bead]
    bonds: List[Tuple[int, int, str]]                  # (i, j, bond key), 1-based, i<j
    angles: List[Tuple[int, int, int, str, str]]       # (i, j, k, specific key, generic key)
    exclusions: List[Tuple[int, int]]                  # 1-2 and 1-3 pairs, 1-based, i<j
    n_atoms: int
    atom_bead: np.ndarray = field(repr=False, default=None)   # (n_atoms,) bead index 0-based
    atom_mass: np.ndarray = field(repr=False, default=None)

    @property
    def types(self) -> List[str]:
        return [b.type for b in self.beads]

    @property
    def masses(self) -> np.ndarray:
        return np.array([b.mass for b in self.beads])


def _bond_key(cls: str, ta: str, tb: str) -> str:
    return f"{cls}:{pair_key(ta, tb)}"


def map_chain(topo: ChainTopology, scheme: str = "ring-side") -> CGMolecule:
    """Build the bead graph of one chain from its atomistic topology."""
    keys: Dict[Tuple[int, str], str] = {}                      # (resid, key) → type
    atom_keys: List[Tuple[int, str]] = []
    for a in topo.atoms:
        k, t = _bead_of_atom(scheme, a.name, a.resname)
        keys[(a.resid, k)] = t
        atom_keys.append((a.resid, k))
    order = sorted(keys, key=lambda rk: (rk[0], "RSM".index(rk[1])))
    index = {rk: n for n, rk in enumerate(order)}
    atom_bead = np.array([index[rk] for rk in atom_keys])
    atom_mass = np.array([MASS[a.element] for a in topo.atoms])
    beads = []
    for n, (resid, k) in enumerate(order):
        members = np.where(atom_bead == n)[0]
        beads.append(Bead(n + 1, keys[(resid, k)], resid, k,
                          float(atom_mass[members].sum()), members.tolist()))

    link_of = {}
    for lk in topo.linkages:
        link_of[tuple(sorted((lk["residue_1"], lk["residue_2"])))] = lk["linkage"]

    pairs: Dict[Tuple[int, int], str] = {}
    for i, j in topo.bonds:
        bi, bj = atom_bead[i - 1], atom_bead[j - 1]
        if bi == bj:
            continue
        a, b = sorted((int(bi), int(bj)))
        if (a, b) in pairs:
            continue
        ra, rb = beads[a].resid, beads[b].resid
        cls = "intra" if ra == rb else link_of.get(tuple(sorted((ra, rb))), "inter")
        pairs[(a, b)] = _bond_key(cls, beads[a].type, beads[b].type)
    bonds = [(a + 1, b + 1, k) for (a, b), k in sorted(pairs.items())]

    nbrs: Dict[int, Dict[int, str]] = {}
    for a, b, k in bonds:
        nbrs.setdefault(a, {})[b] = k
        nbrs.setdefault(b, {})[a] = k
    angles = []
    for j in sorted(nbrs):
        for a, b in itertools.combinations(sorted(nbrs[j]), 2):
            ka, kb = nbrs[j][a], nbrs[j][b]
            tj = beads[j - 1].type
            specific = f"{tj}:" + "|".join(sorted((ka, kb)))
            generic = f"{pair_key(beads[a - 1].type, beads[b - 1].type)}@{tj}"
            angles.append((a, j, b, specific, generic))

    excl = set()
    for a, b, _ in bonds:
        excl.add((a, b))
    for j in nbrs:
        for a, b in itertools.combinations(sorted(nbrs[j]), 2):
            excl.add((a, b))
    return CGMolecule(topo.name, scheme, beads, bonds, angles, sorted(excl),
                      len(topo.atoms), atom_bead, atom_mass)


# ── parameters ────────────────────────────────────────────────────────────────

@dataclass
class CGParams:
    """
    Fitted parameters of the CG model.

    nonbonded
        ``{"RG-SC": {"sigma": nm, "epsilon": kJ/mol}}`` (12-6 Lennard-Jones,
        potential shifted at ``cutoff``); every pair of bead types present.
    bonds
        ``{"beta-O-4:RG-SC": {"r0": nm, "k": kJ/mol/nm2}}``; generic fallbacks
        are keyed by the plain type pair ``"RG-SC"``.
    angles
        ``{"RG:intra:RG-SC|beta-O-4:RG-SC": {"theta0": deg, "k": kJ/mol/rad2}}``;
        generic fallbacks keyed ``"SC-SC@RG"``.
    """

    scheme: str
    temperature: float
    cutoff: float
    nonbonded: Dict[str, Dict[str, float]] = field(default_factory=dict)
    bonds: Dict[str, Dict[str, float]] = field(default_factory=dict)
    angles: Dict[str, Dict[str, float]] = field(default_factory=dict)
    default_bond: Dict[str, float] = field(default_factory=lambda: {"r0": 0.45, "k": 3000.0})
    default_angle: Dict[str, float] = field(default_factory=lambda: {"theta0": 130.0, "k": 40.0})
    provenance: str = ""

    def bond(self, key: str) -> Dict[str, float]:
        if key in self.bonds:
            return self.bonds[key]
        generic = key.split(":", 1)[1]
        return self.bonds.get(generic, self.default_bond)

    def angle(self, specific: str, generic: str) -> Dict[str, float]:
        return self.angles.get(specific) or self.angles.get(generic) or self.default_angle

    def to_dict(self) -> dict:
        return {"format": "lignoforge-cgfit-1", **self.__dict__}

    @classmethod
    def from_dict(cls, d: dict) -> "CGParams":
        d = {k: v for k, v in d.items() if k != "format"}
        return cls(**d)

    def save(self, path: str) -> None:
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=2)

    @classmethod
    def load(cls, path: str) -> "CGParams":
        with open(path) as fh:
            return cls.from_dict(json.load(fh))


def all_type_pairs(types: Sequence[str]) -> List[str]:
    ts = sorted(set(types))
    return [pair_key(a, b) for a, b in itertools.combinations_with_replacement(ts, 2)]
