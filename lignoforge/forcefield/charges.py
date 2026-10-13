"""
Partial-charge correction by junction.

Fragment charges from :mod:`lignoforge.forcefield.opls` are the OPLS-AA charges of
the reference groups.  A residue is *not* a neutral entity: where a bond joins two
residues (an inter-monomer linkage) or a group is attached (a derivative), the
neutral OPLS group is split across the two sides, so individual residues carry
charges such as ±0.2 e that are physical, not errors.

What really goes wrong is smaller and has one cause each: every junction type leaves
a fixed excess charge ``e_J`` because the OPLS fragments on its two sides do not add
up to a neutral whole.  This excess is **additive** over junctions and depends only on
the junction type (verified over random chains in the test-suite), so

* the chain total of the raw charges is exactly ``sum(e_J)`` over its junctions;
* the correction ``-e_J`` is applied **locally**, spread uniformly (the minimum-norm
  correction) over the atoms of the junction and their directly bonded neighbours;
* no residue is forced to be neutral.

Forcing each residue to zero by a uniform shift was tried first and rejected: it moves
the ±0.2 e of the split ether group onto the whole residue, which adds a spurious
dipole of 4-6 D across a β-O-4 dimer and doubles the electrostatic-potential error with
respect to HF/6-31G* (see ``docs`` → force field → partial charges).

The excess values were obtained by least squares of the raw chain charge on the
junction counts (60 random chains, zero residual) and are checked by
``tests/test_forcefield.py::test_junction_excess_is_additive``.  If a chain still has an
unexplained remainder (a chemistry not in the tables) it is spread over all junction
atoms and reported as ``unexplained_residual``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Sequence, Set, Tuple

__all__ = [
    "LINKAGE_EXCESS", "DERIVATIVE_EXCESS", "DEFAULT_DECIMALS", "Junction",
    "correct_charges", "junction_atoms",
]

DEFAULT_DECIMALS = 4

# Excess charge (e) of the raw OPLS fragments per inter-monomer bond of the graph
# (ring linkages have several bonds: β-5 two, β-β three).
LINKAGE_EXCESS: Dict[str, float] = {
    "beta-O-4": 0.0,
    "alpha-O-4": -0.025,
    "4-O-5": -0.115,
    "5-5": 0.0,
    "beta-5": -0.0125,
    "beta-beta": 0.103333333,
    "beta-1": 0.0,          # disabled in the builder
}

# Excess (e) per derivatised hydroxyl, keyed by (site letter, group kind):
# P phenolic, A benzylic (α), G γ; kind 'methyl' or 'ester'.
DERIVATIVE_EXCESS: Dict[Tuple[str, str], float] = {
    ("P", "methyl"): 0.0, ("P", "ester"): 0.0,
    ("A", "methyl"): -0.025, ("A", "ester"): -0.025,
    ("G", "methyl"): 0.0, ("G", "ester"): 0.06,
}


@dataclass
class Junction:
    """A place where the raw charges are corrected."""
    kind: str                       # 'linkage' | 'derivative'
    label: str                      # e.g. 'beta-O-4' or 'G ester'
    core: List[int]                 # atom indices (1-based) of the bond that defines it
    excess: float                   # raw excess charge of this junction (e)
    atoms: List[int] = field(default_factory=list)   # atoms receiving the correction
    shift_per_atom: float = 0.0


def junction_atoms(core: Iterable[int], adj: Dict[int, Set[int]]) -> List[int]:
    """The core atoms and every atom directly bonded to them."""
    out: Set[int] = set(core)
    for c in list(out):
        out |= adj.get(c, set())
    return sorted(out)


def correct_charges(
    raw: Dict[int, float],
    junctions: Sequence[Junction],
    adj: Dict[int, Set[int]],
    decimals: int = DEFAULT_DECIMALS,
) -> Tuple[Dict[int, float], float]:
    """
    Return ``(charges, unexplained)``: charges (rounded to *decimals*) summing to
    exactly zero, and the part of the chain excess not accounted for by the junction
    tables.  The *junctions* are updated in place with their atoms and shifts.
    """
    scale = 10 ** decimals
    delta: Dict[int, float] = {i: 0.0 for i in raw}
    total = sum(raw.values())
    explained = 0.0
    for j in junctions:
        j.atoms = junction_atoms(j.core, adj)
        j.shift_per_atom = -j.excess / len(j.atoms)
        for a in j.atoms:
            delta[a] += j.shift_per_atom
        explained += j.excess
    unexplained = total - explained
    touched = sorted({a for j in junctions for a in j.atoms})
    if abs(unexplained) > 0.5 / scale:
        pool = touched or sorted(raw)
        for a in pool:
            delta[a] -= unexplained / len(pool)
    else:
        unexplained = 0.0

    units = {i: round((raw[i] + delta[i]) * scale) for i in raw}
    leftover = -sum(units.values())
    pool = touched or sorted(raw)
    order = sorted(pool, key=lambda a: abs(units[a]))      # least polar atoms first
    step = 1 if leftover > 0 else -1
    for k in range(abs(leftover)):
        units[order[k % len(order)]] += step
    return {i: units[i] / scale for i in raw}, unexplained
