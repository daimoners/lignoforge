"""
Coarse-grained force-field parameters (one bead per monomer).

The parameter set is a plain, serialisable object so that parameters derived
from atomistic simulations (:mod:`lignoforge.cg.boltzmann`) can replace the
defaults without code changes.

Defaults are **provisional** and intended as a starting point:

* bond ``r0`` per linkage: mean centre-of-mass distance over 60 MMFF-relaxed
  6-mers (n = 14-149 bonds per linkage);
* angle ``theta0`` for the linkage pairs sampled at least 10 times in that
  ensemble; all other pairs use the generic value;
* force constants, bead sizes and well depths are *not* derived from data.
  Replace them with Boltzmann-inversion results from atomistic MD before
  using the model quantitatively.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Dict, Optional

__all__ = ["CGParameters", "default_parameters"]

PROVENANCE_DEFAULT = (
    "LignoForge provisional CG defaults: r0/theta0 from MMFF-relaxed 6-mer "
    "ensemble; k, sigma, epsilon are generic placeholders."
)


def _angle_key(l1: str, l2: str) -> str:
    return "|".join(sorted((l1, l2)))


@dataclass
class CGParameters:
    """
    Parameters of the monomer-bead model.

    Attributes
    ----------
    beads
        ``{"H"|"G"|"S": {"sigma": nm, "epsilon": kJ/mol}}``.
    bonds
        ``{key: {"r0": nm, "k": kJ/mol/nm^2}}``.  Key is the linkage name, or
        ``"<linkage>:<A>-<B>"`` (bead types sorted) for a type-specific entry.
    angles
        ``{key: {"theta0": deg, "k": kJ/mol/rad^2}}``.  Key is the two
        linkage names sorted and joined by ``|``.
    default_bond / default_angle
        Fallbacks used when a key is not found.
    """

    beads: Dict[str, Dict[str, float]]
    bonds: Dict[str, Dict[str, float]]
    angles: Dict[str, Dict[str, float]]
    default_bond: Dict[str, float] = field(
        default_factory=lambda: {"r0": 0.62, "k": 2500.0})
    default_angle: Dict[str, float] = field(
        default_factory=lambda: {"theta0": 140.0, "k": 30.0})
    provenance: str = PROVENANCE_DEFAULT

    def bond(self, linkage: str, type_a: str, type_b: str) -> Dict[str, float]:
        specific = f"{linkage}:{'-'.join(sorted((type_a, type_b)))}"
        return self.bonds.get(specific) or self.bonds.get(linkage) or self.default_bond

    def angle(self, link_1: str, link_2: str) -> Dict[str, float]:
        return self.angles.get(_angle_key(link_1, link_2)) or self.default_angle

    # ── (de)serialisation ─────────────────────────────────────────────────────
    def to_dict(self) -> dict:
        return {
            "provenance": self.provenance, "beads": self.beads,
            "bonds": self.bonds, "angles": self.angles,
            "default_bond": self.default_bond, "default_angle": self.default_angle,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "CGParameters":
        return cls(
            beads=d["beads"], bonds=d.get("bonds", {}), angles=d.get("angles", {}),
            default_bond=d.get("default_bond", {"r0": 0.62, "k": 2500.0}),
            default_angle=d.get("default_angle", {"theta0": 140.0, "k": 30.0}),
            provenance=d.get("provenance", ""),
        )

    def save(self, path: str) -> None:
        with open(path, "w") as fh:
            json.dump(self.to_dict(), fh, indent=2)

    @classmethod
    def load(cls, path: str) -> "CGParameters":
        with open(path) as fh:
            return cls.from_dict(json.load(fh))


def default_parameters() -> CGParameters:
    """Return a fresh copy of the provisional default parameter set."""
    k_b, k_a = 2500.0, 30.0
    r0 = {"4-O-5": 0.604, "5-5": 0.574, "alpha-O-4": 0.633,
          "beta-5": 0.608, "beta-O-4": 0.656, "beta-beta": 0.663}
    theta0 = {"alpha-O-4|beta-O-4": 147.8, "beta-O-4|beta-O-4": 142.0,
              "beta-5|beta-O-4": 153.8, "beta-O-4|beta-beta": 156.0,
              "4-O-5|beta-O-4": 109.5, "5-5|beta-O-4": 92.5,
              "alpha-O-4|beta-beta": 151.8}
    return CGParameters(
        beads={
            "H": {"sigma": 0.62, "epsilon": 3.0},
            "G": {"sigma": 0.66, "epsilon": 3.0},
            "S": {"sigma": 0.70, "epsilon": 3.0},
        },
        bonds={l: {"r0": r, "k": k_b} for l, r in r0.items()},
        angles={l: {"theta0": t, "k": k_a} for l, t in theta0.items()},
    )
