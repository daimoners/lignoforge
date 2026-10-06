"""
OPLS-AA atom typing for lignin chains.

Typing is driven by two pieces of information only, both available from the
graph that LignoForge generates:

* the set of atom names present in a residue (RTP-style names produced by
  :mod:`lignoforge.structure.generator`), and
* the *linkage context* of the residue: which atoms take part in an
  inter-monomer bond, and of which linkage type.

Nothing here depends on residue order or on residue adjacency, so linear,
branched and ring-closed (β-5, β-β) chains are treated identically.

Charges returned here are the *raw* fragment charges.  Per-residue
neutralisation is a separate step (:mod:`lignoforge.forcefield.charges`).
"""

from __future__ import annotations

from typing import Dict, Set, Tuple

__all__ = [
    "TypingError",
    "LINKAGE_ATOM_PAIRS",
    "EXTRA_DIHEDRALTYPES",
    "linkage_position",
    "type_residue",
]


class TypingError(ValueError):
    """Raised when a residue cannot be typed unambiguously."""


# (opls_type, charge)
AtomType = Tuple[str, float]

# ── Allowed inter-monomer atom pairs per linkage ──────────────────────────────
# Every inter-monomer bond of the graph must match one of these; anything else
# is a structure/typing mismatch and is reported instead of being guessed.
LINKAGE_ATOM_PAIRS: Dict[str, Set[frozenset]] = {
    "beta-O-4":  {frozenset(("CB", "O4H"))},
    "alpha-O-4": {frozenset(("CA", "O4H"))},
    "4-O-5":     {frozenset(("C5", "O4H"))},
    "5-5":       {frozenset(("C5",))},
    "beta-5":    {frozenset(("CB", "C5")), frozenset(("CA", "O4H"))},
    "beta-beta": {frozenset(("CB",)), frozenset(("CA", "OG"))},
    "beta-1":    {frozenset(("CB", "C1"))},
}

# Secondary (ring-closing) bond of ring linkages.
_RING_CLOSURE_PAIR = {
    "beta-5":    frozenset(("CA", "O4H")),
    "beta-beta": frozenset(("CA", "OG")),
}


def linkage_position(linkage: str, atom_name: str, partner_name: str) -> str:
    """
    Return the position label of *atom_name* inside *linkage*.

    The label is the atom name itself for the primary bond of a linkage and
    ``<name>_ringclose`` for the secondary bond of β-5 / β-β.

    Raises
    ------
    TypingError
        If the (atom, partner) pair is not valid for *linkage*.
    """
    allowed = LINKAGE_ATOM_PAIRS.get(linkage)
    if allowed is None:
        raise TypingError(f"Unsupported linkage type '{linkage}'")
    if frozenset((atom_name, partner_name)) not in allowed:
        raise TypingError(
            f"Atom pair {atom_name}–{partner_name} is not valid for "
            f"linkage '{linkage}'"
        )
    if frozenset((atom_name, partner_name)) == _RING_CLOSURE_PAIR.get(linkage):
        return f"{atom_name}_ringclose"
    return atom_name


# ── Fragment charges (OPLS-AA database values unless noted) ───────────────────
_ARO_C      = ("opls_145", -0.115)   # aromatic C-H
_ARO_H      = ("opls_146", +0.115)
_C_IPSO     = ("opls_221", -0.055)   # C1 of sp3 side chains
_C_IPSO_VIN = ("opls_221",  0.000)   # C1 of vinyl side chains (self-neutral)
_C_PHENOL   = ("opls_166", +0.150)
_O_PHENOL   = ("opls_167", -0.585)
_H_PHENOL   = ("opls_168", +0.435)
_C_ARYLETH  = ("opls_199", +0.085)   # ring C bearing OMe / ether O
_O_ARYLETH  = ("opls_179", -0.285)
_C_OME      = ("opls_181", +0.110)
_H_OME      = ("opls_185", +0.030)
_CM         = ("opls_142", -0.115)   # vinyl C
_HM         = ("opls_144", +0.115)
_CA_ALC     = ("opls_219", +0.260)   # Cα, benzyl alcohol
_HA_SP3     = ("opls_156", +0.060)
_OH_ALC     = ("opls_154", -0.683)
_HO_ALC     = ("opls_155", +0.418)
# Cβ/Cα of an ether: opls_183 (+0.170) with H opls_185 (+0.030).  The pair
# sums to +0.200, which exactly cancels the −0.200 of the aryl-ether C4–O
# pair, so every C–O–C bridge built from these types is neutral by itself.
_C_ETHER    = ("opls_183", +0.170)
_H_ETHER    = ("opls_185", +0.030)
_C_CH       = ("opls_137", -0.060)   # sp3 CH bonded to carbons only
_H_CH       = ("opls_140", +0.060)
_CG         = ("opls_157", +0.145)
_HG         = ("opls_156", +0.060)

_CT_CH2     = ("opls_136", -0.120)   # sp3 CH2 bonded to carbons only
_H_CH2      = ("opls_140", +0.060)

# Torsions absent from OPLS-AA ffbonded.itp but required by lignin side
# chains:
#  * CA-CA-CM-HC (aryl–vinyl, H on the vinyl carbon) is set to zero, exactly
#    as OPLS does for the analogous CA-CA-CT-HC (ethylbenzene).
#  * CA-CT-CT-CA (Ar–C–C–Ar, β-5 phenylcoumaran) takes the hydrocarbon
#    values, as OPLS does for CA-CT-CT-CT.
EXTRA_DIHEDRALTYPES = (
    "[ dihedraltypes ]\n"
    "; i    j    k    l   func    C0 .. C5 (Ryckaert-Bellemans)\n"
    "  CA   CA   CM   HC   3   0.0 0.0 0.0 0.0 0.0 0.0 "
    "; LignoForge: analogue of CA CA CT HC (ethylbenzene)\n"
    "  CA   CT   CT   CA   3   2.92880 -1.46440 0.20920 -1.67360 0.0 0.0 "
    "; LignoForge: analogue of CA CT CT CT\n"
)


def type_residue(
    names: Set[str],
    context: Dict[str, Tuple[str, str]],
) -> Dict[str, AtomType]:
    """
    Assign OPLS-AA types and raw charges to every atom of one residue.

    Parameters
    ----------
    names
        Atom names present in the residue.
    context
        ``{atom_name: (linkage, position)}`` for the atoms that take part in
        an inter-monomer bond (see :func:`linkage_position`).

    Returns
    -------
    dict
        ``{atom_name: (opls_type, charge)}`` covering *all* of ``names``.

    Raises
    ------
    TypingError
        If an atom name is not recognised.
    """
    n = names
    typed: Dict[str, AtomType] = {}

    def lk(atom: str) -> Tuple[str, str]:
        return context.get(atom, (None, None))

    # ── Side-chain topology flags (decide C1 and Cα/Cβ types) ─────────────────
    ca_ether = "CA" in n and lk("CA")[0] in ("alpha-O-4", "beta-5", "beta-beta")
    ca_alcohol = "CA" in n and "OA" in n and not ca_ether
    ca_sp3 = ca_ether or ca_alcohol
    if ca_ether and "OA" in n:
        raise TypingError(
            "Cα carries both an α-OH and an ether bond (two oxygens on one "
            "carbon); the structure is chemically invalid"
        )

    # ── Ring ──────────────────────────────────────────────────────────────────
    if "C1" in n:
        typed["C1"] = _C_IPSO if ca_sp3 else _C_IPSO_VIN
    for ring_c, ring_h in (("C2", "H2"), ("C6", "H6")):
        if ring_c in n:
            typed[ring_c] = _ARO_C
        if ring_h in n:
            typed[ring_h] = _ARO_H

    # C3 / C5: aromatic CH, OMe-bearing, or linkage-substituted
    for pos, ome_o in (("3", "OM3"), ("5", "OM5")):
        c, h = f"C{pos}", f"H{pos}"
        if c not in n:
            continue
        if ome_o in n:
            typed[c] = _C_ARYLETH
        elif h in n:
            typed[c] = _ARO_C
        elif pos == "5" and lk("C5")[0] == "4-O-5":
            typed[c] = _C_ARYLETH
        elif pos == "5" and lk("C5")[0] in ("5-5", "beta-5"):
            typed[c] = ("opls_145", 0.000)
        else:
            raise TypingError(f"{c} has neither H nor a known substituent")
        if h in n:
            typed[h] = _ARO_H
    for ome in ("3", "5"):
        if f"OM{ome}" in n:
            typed[f"OM{ome}"] = _O_ARYLETH
        if f"CM{ome}" in n:
            typed[f"CM{ome}"] = _C_OME
        for i in (1, 2, 3):
            if f"HM{ome}{i}" in n:
                typed[f"HM{ome}{i}"] = _H_OME

    # C4 / O4H / HO4
    if "C4" in n:
        typed["C4"] = _C_PHENOL if "HO4" in n else _C_ARYLETH
    if "O4H" in n:
        typed["O4H"] = _O_PHENOL if "HO4" in n else _O_ARYLETH
    if "HO4" in n:
        typed["HO4"] = _H_PHENOL

    # ── Cα ────────────────────────────────────────────────────────────────────
    if "CA" in n:
        if ca_ether:
            typed["CA"] = _C_ETHER
            if "HA" in n:
                typed["HA"] = _HA_SP3
        elif ca_alcohol:
            typed["CA"] = _CA_ALC
            if "HA" in n:
                typed["HA"] = _HA_SP3
        else:
            typed["CA"] = _CM
            if "HA" in n:
                typed["HA"] = _HM
    if "OA" in n:
        typed["OA"] = _OH_ALC
    if "HOA" in n:
        typed["HOA"] = _HO_ALC

    # ── Cβ ────────────────────────────────────────────────────────────────────
    if "CB" in n:
        link = lk("CB")[0]
        if link == "beta-O-4":
            typed["CB"] = _C_ETHER
            if "HB" in n:
                typed["HB"] = _H_ETHER
        elif link in ("beta-5", "beta-beta", "beta-1"):
            typed["CB"] = _C_CH
            if "HB" in n:
                typed["HB"] = _H_CH
        elif "HB1" in n and "HB2" in n:
            typed["CB"] = _CT_CH2          # e.g. Cβ left as CH2 after α-O-4
            typed["HB1"] = typed["HB2"] = _H_CH2
        else:
            typed["CB"] = _CM
            if "HB" in n:
                typed["HB"] = _HM

    # ── Cγ ────────────────────────────────────────────────────────────────────
    if "CG" in n:
        typed["CG"] = _CG
    for h in ("HG1", "HG2"):
        if h in n:
            typed[h] = _HG
    if "OG" in n:
        if lk("OG")[0] == "beta-beta":
            typed["OG"] = _O_ARYLETH       # tetrahydrofuran-type ether O
        else:
            typed["OG"] = _OH_ALC
            if "HOG" in n:
                typed["HOG"] = _HO_ALC

    unknown = n - typed.keys()
    if unknown:
        raise TypingError(f"Unrecognised atom name(s): {sorted(unknown)}")
    return typed
