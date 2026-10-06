"""
Partial-charge neutralisation.

Fragment charges from :mod:`lignoforge.forcefield.opls` are exact for the
reference fragments but a residue that carries one or more linkages can end
up with a small net charge (for example a β-O-4 acceptor-only end unit).
Each residue is made exactly neutral by spreading its residual uniformly over
*all* of its atoms.  Because every residue is neutralised on its own, the
chain total is zero by construction, whatever the topology (linear, branched
or ring-closed) and however the residues are numbered.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

__all__ = ["neutralize_residue", "DEFAULT_DECIMALS"]

DEFAULT_DECIMALS = 4


def neutralize_residue(
    charges: Dict[str, float],
    decimals: int = DEFAULT_DECIMALS,
) -> Tuple[Dict[str, float], float, float]:
    """
    Return charges rounded to *decimals* places that sum to exactly zero.

    Parameters
    ----------
    charges
        ``{atom_name: raw_charge}`` for one residue.

    Returns
    -------
    (new_charges, residual, max_shift)
        ``residual`` is the net charge before neutralisation and
        ``max_shift`` the largest absolute change applied to any atom.
    """
    if not charges:
        return {}, 0.0, 0.0

    scale = 10 ** decimals
    names: List[str] = list(charges)
    residual = sum(charges.values())
    per_atom = residual / len(names)

    # Work in integer units of 10**-decimals so the sum is exactly zero.
    units = {k: round((charges[k] - per_atom) * scale) for k in names}
    leftover = -sum(units.values())
    # Give the rounding leftover (a few units at most) to the atoms with the
    # smallest |q|, i.e. the least polar ones.
    for k in sorted(names, key=lambda a: abs(units[a]))[: abs(leftover)]:
        units[k] += 1 if leftover > 0 else -1

    new = {k: units[k] / scale for k in names}
    max_shift = max(abs(new[k] - charges[k]) for k in names)
    return new, residual, max_shift
