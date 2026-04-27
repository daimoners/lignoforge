#!/usr/bin/env python3
"""
lignin_ff/tools/assign_chain_types.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Assign OPLS-AA atom types to a lignin polymer-chain PDB and write
per-chain custom topology RTP files.

Usage
-----
    python assign_chain_types.py chain.pdb [-o chain_custom.rtp]

The input PDB must use RTP-compatible atom names (as produced by the
lignoforge PDB writer with the updated generator.py).  The output RTP
can be used directly with GROMACS pdb2gmx together with the
``oplsaa.ff`` force field and ``residuetypes_lignin.dat``.

Design principles
-----------------
* **Deterministic** – given the same PDB the output is always identical.
* **Maximally stable types** – atom types are changed only at atoms
  directly involved in a linkage.
* **Per-residue net charge = 0.000 e** – any imbalance introduced by
  removed/added atoms at linkage sites is absorbed on C1 (ipso), which
  has no formal OPLS-AA charge constraint.
* Supports all 7 lignin linkage types:
  beta-O-4, 5-5, 4-O-5, beta-5, beta-beta, alpha-O-4, beta-1.
* Identical residue environments (same monomer type + same linkage
  fingerprint) are assigned the same canonical RTP name; unique
  environments receive a generated name (e.g. ``GY001``).
"""

from __future__ import annotations

import argparse
import collections
import os
import sys
from typing import Dict, List, Optional, Tuple

# ── Residue-code helpers ───────────────────────────────────────────────────────

# Map PDB residue names → H/G/S monomer base type
RES_TO_BASE = {
    "HPU": "H", "HNM": "H", "HPM": "H",
    "GYU": "G", "GNM": "G", "GYM": "G",
    "SYU": "S", "SNM": "S", "SYM": "S",
}
BASE_TO_NM   = {"H": "HNM", "G": "GNM", "S": "SNM"}

# ── OPLS-AA atom type records ──────────────────────────────────────────────────
# (oplsaa_code, gromacs_type, raw_charge)

_T = Dict[str, Tuple[str, str, float]]

ATYPE: _T = {
    "CA_ipso":    ("opls_145", "CA",  0.000),
    "CA_plain":   ("opls_145", "CA", -0.115),
    "HA_arom":    ("opls_146", "HA",  0.115),
    "CA_phenol":  ("opls_166", "CA",  0.150),
    "OH_phenol":  ("opls_167", "OH", -0.585),
    "HO_phenol":  ("opls_168", "HO",  0.435),
    "CA_ether":   ("opls_199", "CA",  0.085),
    "OS_ether":   ("opls_179", "OS", -0.285),
    "CT_methoxy": ("opls_181", "CT",  0.110),
    "HC_methoxy": ("opls_185", "HC",  0.030),
    "CT_alc":     ("opls_157", "CT",  0.145),   # raw; adjusted in neutral monomer defs
    "HC_alc":     ("opls_156", "HC",  0.040),
    "OH_alc":     ("opls_154", "OH", -0.683),
    "HO_alc":     ("opls_155", "HO",  0.418),
    "CT_CH2":     ("opls_136", "CT", -0.120),   # sp3 CH2, no heteroatom
    "HC_alkyl":   ("opls_140", "HC",  0.060),
    "CT_CH":      ("opls_137", "CT", -0.060),   # sp3 CH (tertiary/secondary)
    "CM_vinyl":   ("opls_142", "CM", -0.115),
    "HC_vinyl":   ("opls_144", "HC",  0.115),
}


def at(key: str) -> Tuple[str, str, float]:
    """Short-hand lookup for an ATYPE entry (opls_code, grmcs_type, charge)."""
    return ATYPE[key]


# ── Isolated neutral monomer atom definitions ──────────────────────────────────
# Layout per atom:
#   atom_name -> (opls_code, gromacs_type, charge, cgnr, comment)

def _ring_atoms(base: str) -> dict:
    """Aromatic ring atoms (C1–C6 + H), type-varied by H/G/S."""
    a: dict = {}
    a["C1"]  = (*at("CA_ipso"),  1, "Cipso")
    a["C2"]  = (*at("CA_plain"), 2, "aromatic CH")
    a["H2"]  = (*at("HA_arom"),  2, "aromatic H")

    if base == "H":
        a["C3"]  = (*at("CA_plain"), 3, "aromatic CH")
        a["H3"]  = (*at("HA_arom"),  3, "aromatic H")
    else:
        a["C3"]   = (*at("CA_ether"),   3, "C3 bearing OMe")
        a["OM3"]  = (*at("OS_ether"),   3, "OMe O at C3")
        a["CM3"]  = (*at("CT_methoxy"), 3, "OMe CH3 at C3")
        a["HM31"] = (*at("HC_methoxy"), 3, "OMe H")
        a["HM32"] = (*at("HC_methoxy"), 3, "OMe H")
        a["HM33"] = (*at("HC_methoxy"), 3, "OMe H")

    a["C4"]  = (*at("CA_phenol"), 4, "C4 bearing phenol OH")
    a["O4H"] = (*at("OH_phenol"), 4, "phenol O")
    a["HO4"] = (*at("HO_phenol"), 4, "phenol H")

    if base in ("H", "G"):
        a["C5"]  = (*at("CA_plain"), 5, "aromatic CH")
        a["H5"]  = (*at("HA_arom"),  5, "aromatic H")
    else:  # S
        a["C5"]   = (*at("CA_ether"),   5, "C5 bearing OMe")
        a["OM5"]  = (*at("OS_ether"),   5, "OMe O at C5")
        a["CM5"]  = (*at("CT_methoxy"), 5, "OMe CH3 at C5")
        a["HM51"] = (*at("HC_methoxy"), 5, "OMe H")
        a["HM52"] = (*at("HC_methoxy"), 5, "OMe H")
        a["HM53"] = (*at("HC_methoxy"), 5, "OMe H")

    a["C6"]  = (*at("CA_plain"), 6, "aromatic CH")
    a["H6"]  = (*at("HA_arom"),  6, "aromatic H")
    return a


def _sidechain_atoms() -> dict:
    """sp3 side-chain atoms for isolated neutral monomer (Cα-OH, Cβ-CH2, Cγ-OH)."""
    # CA charge adjusted to +0.205 so that cgnr7 = CA+HA+OA+HOA sums to 0.
    opls7, gt7, _ = at("CT_alc")
    _, _, hac  = at("HC_alc")
    _, _, oac  = at("OH_alc")
    _, _, hoac = at("HO_alc")
    # CA adjusted: +0.205 so 0.205+0.060-0.683+0.418 = 0.000
    a = {
        "CA":  (opls7,        gt7,   0.205,  7, "Cα sp3 secondary alcohol (adjusted)"),
        "HA":  (*at("HC_alc"),       7, "H on Cα"),
        "OA":  (*at("OH_alc"),       7, "α-OH O"),
        "HOA": (*at("HO_alc"),       7, "α-OH H"),
        "CB":  (*at("CT_CH2"),       8, "Cβ sp3 CH2"),
        "HB1": (*at("HC_alkyl"),     8, "H on Cβ"),
        "HB2": (*at("HC_alkyl"),     8, "H on Cβ"),
        "CG":  (*at("CT_alc"),       9, "Cγ CH2-OH"),
        "HG1": (*at("HC_alc"),       9, "H on Cγ"),
        "HG2": (*at("HC_alc"),       9, "H on Cγ"),
        "OG":  (*at("OH_alc"),       9, "γ-OH O"),
        "HOG": (*at("HO_alc"),       9, "γ-OH H"),
    }
    # Fix cgnr7 for HA: use +0.060 (1propanol value, already in HC_alc raw = 0.040,
    # but the RTP convention for lignin uses 0.060 throughout)
    _op, _gt, _ = at("HC_alc")
    a["HA"]  = (_op, _gt, 0.060, 7, "H on Cα")
    a["HG1"] = (_op, _gt, 0.060, 9, "H on Cγ")
    a["HG2"] = (_op, _gt, 0.060, 9, "H on Cγ")
    return a


# Build the full isolated neutral monomer dictionaries
ISOLATED: Dict[str, dict] = {}
for _base, _nm in [("H", "HNM"), ("G", "GNM"), ("S", "SNM")]:
    ISOLATED[_nm] = {**_ring_atoms(_base), **_sidechain_atoms()}


# ── Intra-residue bond lists for isolated neutral monomers ────────────────────

def _ring_bonds(base: str) -> list:
    b = [
        ("C1","C2"),("C2","C3"),("C3","C4"),("C4","C5"),("C5","C6"),("C6","C1"),
        ("C2","H2"),("C6","H6"),
    ]
    if base == "H":
        b += [("C3","H3"),("C5","H5")]
    else:
        b += [("C3","OM3"),("OM3","CM3"),("CM3","HM31"),("CM3","HM32"),("CM3","HM33")]
        if base == "G":
            b += [("C5","H5")]
        else:  # S
            b += [("C5","OM5"),("OM5","CM5"),("CM5","HM51"),("CM5","HM52"),("CM5","HM53")]
    b += [("C4","O4H"),("O4H","HO4")]
    return b


def _sidechain_bonds() -> list:
    return [
        ("C1","CA"),("CA","HA"),("CA","OA"),("OA","HOA"),
        ("CA","CB"),("CB","HB1"),("CB","HB2"),
        ("CB","CG"),("CG","HG1"),("CG","HG2"),("CG","OG"),("OG","HOG"),
    ]


def _ring_impropers(base: str) -> list:
    imp = [
        ("C2","C6","C1","CA"),
        ("C1","C3","C2","H2"),
    ]
    if base == "H":
        imp += [("C2","C4","C3","H3"), ("C4","C6","C5","H5")]
    else:
        imp += [("C2","C4","C3","OM3")]
        if base == "G":
            imp += [("C4","C6","C5","H5")]
        else:
            imp += [("C4","C6","C5","OM5")]
    imp += [("C3","C5","C4","O4H"), ("C5","C1","C6","H6")]
    return imp


ISOLATED_BONDS: Dict[str, list] = {}
ISOLATED_IMPROPERS: Dict[str, list] = {}
for _base, _nm in [("H", "HNM"), ("G", "GNM"), ("S", "SNM")]:
    ISOLATED_BONDS[_nm]     = _ring_bonds(_base) + _sidechain_bonds()
    ISOLATED_IMPROPERS[_nm] = _ring_impropers(_base)


# ── Linkage modification rules ────────────────────────────────────────────────
#
# LINKAGE_MODS[linkage_name][atom_name_in_inter_residue_bond] = {
#   "type_changes": {atom: (opls, gtype, q)},    # atoms that change type/charge
#   "remove":       [atom_names],                # atoms deleted from this residue
#   "rename":       {old: new},                  # atom renames (e.g. HB1→HB)
#   "add":          {atom: (opls, gtype, q, cgnr, desc)},  # new atoms
# }
#
# Charge imbalance after modifications is absorbed on C1 (ipso, charge=0.000).
# All rules are derived from the isolated neutral monomer reference (HNM/GNM/SNM).

LINKAGE_MODS: Dict[str, Dict[str, dict]] = {

    # ── β-O-4: most common lignin linkage (~50-65%) ────────────────────────────
    # Bond formed: O4H (phenol-O of this unit) — CB (β-C of other unit)
    # Donor (provides O4H): phenol C4→ether C4, phenol O→ether O, HO4 removed.
    # Acceptor (provides CB): sp3 CH2→sp3 ether CH (one HB removed, CB: sp3 C-O).
    "beta-O-4": {
        "O4H": {   # this residue is the ether-O DONOR (C4-O4H side)
            "type_changes": {
                "C4":  ("opls_199", "CA",  +0.085),  # phenol C → aryl ether C
                "O4H": ("opls_179", "OS",  -0.285),  # phenol O → aryl ether O
            },
            "remove":  ["HO4"],
            "rename":  {},
            "add":     {},
        },
        "CB": {    # this residue is the ether-C ACCEPTOR (Cβ side)
            "type_changes": {
                "CB":  ("opls_157", "CT",  +0.140),  # CH2 → ether CH (sp3 C-O)
                "HB1": ("opls_156", "HC",  +0.060),  # H on ether CH
            },
            "remove":  ["HB2"],   # CB: CH2 → CH (one H removed)
            "rename":  {"HB1": "HB"},
            "add":     {},
            # Note: OA (α-OH) is already present in the isolated neutral monomer.
        },
    },

    # ── 5-5 (biphenyl): C5-C5 direct aryl-aryl bond ────────────────────────────
    # Both residues contribute C5; H5 removed from each.
    "5-5": {
        "C5": {
            "type_changes": {
                "C5":  ("opls_145", "CA",  0.000),  # CH → substituted C (no H, q=0)
            },
            "remove":  ["H5"],
            "rename":  {},
            "add":     {},
        },
    },

    # ── 4-O-5 (diaryl ether): O bridge between C4 of one and C5 of another ─────
    "4-O-5": {
        "O4H": {   # C4-O4H side (same as beta-O-4 donor)
            "type_changes": {
                "C4":  ("opls_199", "CA",  +0.085),
                "O4H": ("opls_179", "OS",  -0.285),
            },
            "remove":  ["HO4"],
            "rename":  {},
            "add":     {},
        },
        "C5": {    # C5 side: receives the ether O from the other residue
            "type_changes": {
                "C5":  ("opls_199", "CA",  +0.085),  # CH → aryl ether C
            },
            "remove":  ["H5"],
            "rename":  {},
            "add":     {},
        },
    },

    # ── β-5 (phenylcoumaran): Cβ–C5 bond + ring closure O4H–Cα ─────────────────
    # Primary inter-residue bond: CB (unit A) — C5 (unit B).
    # Secondary ring-closure bond: O4H (unit B) — CA (unit A).
    # At C5: loses H5, C5 becomes sp2 substituted aromatic C.
    # At CB: loses one HB, CB becomes sp3 secondary C.
    # At O4H (ring-closure side): C4-O4H end of unit B bonds to CA of unit A
    #   → same phenol→ether change as beta-O-4 donor.
    # At CA (ring-closure side): CA bonds to O4H → COH → C-O-Ar; OA removed.
    "beta-5": {
        "CB": {  # CB of unit A bonds to C5 of unit B
            "type_changes": {
                "CB":  ("opls_137", "CT",  -0.060),  # CH2 → sp3 CH (tertiary-like)
                "HB1": ("opls_140", "HC",  +0.060),
            },
            "remove":  ["HB2"],
            "rename":  {"HB1": "HB"},
            "add":     {},
        },
        "C5": {  # C5 of unit B bonds to CB of unit A
            "type_changes": {
                "C5":  ("opls_145", "CA",  0.000),   # CH → substituted aromatic C
            },
            "remove":  ["H5"],
            "rename":  {},
            "add":     {},
        },
        # Ring closure: O4H of unit B bonds to CA of unit A.
        # This pair of atoms is also captured in the same CONECT scan.
        "O4H_ringclose": {  # same as beta-O-4 donor for C4-O4H
            "_alias": ("beta-O-4", "O4H"),
        },
        "CA_ringclose": {   # CA becomes ether C; OA (α-OH) removed, bond → O4H
            "type_changes": {
                "CA":  ("opls_157", "CT",  +0.140),  # Cα sp3 ether C (like CB in beta-O-4)
                "HA":  ("opls_156", "HC",  +0.060),
            },
            "remove":  ["OA", "HOA"],   # α-OH is replaced by the ether bond
            "rename":  {},
            "add":     {},
        },
    },

    # ── β-β (resinol): Cβ–Cβ bond + two furanoid O-ring closures ───────────────
    # Primary inter-residue bond: CB (unit A) — CB (unit B).
    # Ring closures: CA(A)—OG(B) and CA(B)—OG(A).
    "beta-beta": {
        "CB": {  # each unit contributes CB; one HB removed
            "type_changes": {
                "CB":  ("opls_137", "CT",  -0.060),  # CH2 → sp3 CH
                "HB1": ("opls_140", "HC",  +0.060),
            },
            "remove":  ["HB2"],
            "rename":  {"HB1": "HB"},
            "add":     {},
        },
        # Ring closure: CA bonds to OG of the OTHER unit (inter-residue).
        # CA side: CA changes from alcohol C (CA–OA) to a plain sp3 C with two
        # C-O bonds (the α-OH OA is retained; the new bond is to other-unit OG).
        "CA_ringclose": {
            "type_changes": {
                "CA":  ("opls_157", "CT",  +0.140),  # sp3 C bearing two oxygens
                "HA":  ("opls_156", "HC",  +0.060),
            },
            "remove":  [],
            "rename":  {},
            "add":     {},
        },
        # OG side: OG becomes a furanoid ether O (bridges two C units)
        "OG_ringclose": {
            "type_changes": {
                "OG":  ("opls_179", "OS",  -0.285),  # aliphatic OH → ether O
            },
            "remove":  ["HOG"],
            "rename":  {},
            "add":     {},
        },
    },

    # ── α-O-4: Cα–O–C4 ether (rarer) ───────────────────────────────────────────
    # Bond: O4H (unit A, C4-O end) — CA (unit B, Cα end).
    "alpha-O-4": {
        "O4H": {  # C4-O4H side: same as beta-O-4 donor
            "type_changes": {
                "C4":  ("opls_199", "CA",  +0.085),
                "O4H": ("opls_179", "OS",  -0.285),
            },
            "remove":  ["HO4"],
            "rename":  {},
            "add":     {},
        },
        "CA": {   # Cα becomes ether C: OA (α-OH) is freed/removed, CA→aryl-ether
            "type_changes": {
                "CA":  ("opls_157", "CT",  +0.140),
                "HA":  ("opls_156", "HC",  +0.060),
            },
            "remove":  ["OA", "HOA"],   # α-OH lost; CA now bonds to O4H of other unit
            "rename":  {},
            "add":     {},
        },
    },

    # ── β-1: Cβ–C1 (quinone methide rearrangement) ──────────────────────────────
    # Bond: CB (unit A) — C1 (unit B).
    # In this rearrangement the β-side loses one HB; the C1-side needs no
    # type change (C1 is already ipso at q=0.000).
    "beta-1": {
        "CB": {
            "type_changes": {
                "CB":  ("opls_137", "CT",  -0.060),
                "HB1": ("opls_140", "HC",  +0.060),
            },
            "remove":  ["HB2"],
            "rename":  {"HB1": "HB"},
            "add":     {},
        },
        "C1": {  # C1 receives the new bond; stays aromatic ipso, charge unchanged
            "type_changes": {},
            "remove":  [],
            "rename":  {},
            "add":     {},
        },
    },
}

# Alias resolution
for _lname, _positions in LINKAGE_MODS.items():
    for _pos, _rules in list(_positions.items()):
        if "_alias" in _rules:
            _tgt_ln, _tgt_pos = _rules["_alias"]
            LINKAGE_MODS[_lname][_pos] = LINKAGE_MODS[_tgt_ln][_tgt_pos]


# ── Normal intra-residue bond lookup ──────────────────────────────────────────
# After modifications, bonds involving removed atoms are dropped, and
# modifications to HB1→HB etc are reflected.

INTER_RESIDUE_BOND_PATTERN: Dict[Tuple[str,str], Tuple[str, str, str]] = {
    # (sorted_atom_a, sorted_atom_b): (linkage_name, pos_for_a, pos_for_b)
    # sorted means the pair is alphabetically ordered: a <= b
    ("CB",  "O4H"): ("beta-O-4",   "CB",  "O4H"),
    ("C5",  "C5"):  ("5-5",        "C5",  "C5"),
    ("C5",  "O4H"): ("4-O-5",      "C5",  "O4H"),  # O4H side is the donor
    ("C5",  "OA"):  ("4-O-5",      "C5",  "O4H"),  # fallback if O4H listed as OA
    ("C5",  "CB"):  ("beta-5",     "C5",  "CB"),
    ("CA",  "O4H"): ("alpha-O-4",  "CA",  "O4H"),
    ("CB",  "CB"):  ("beta-beta",  "CB",  "CB"),
    ("C1",  "CB"):  ("beta-1",     "C1",  "CB"),
    # Ring-closure bonds for beta-5 and beta-beta  (secondary bonds)
    ("CA",  "OG"):  ("beta-beta",  "CA_ringclose", "OG_ringclose"),
    ("C4", "OA"):   ("beta-5",     "O4H_ringclose","CA_ringclose"),  # C4 of "O4H"
}

def _sorted_pair(a: str, b: str) -> Tuple[str, str]:
    return (a, b) if a <= b else (b, a)


# ── PDB parsing ───────────────────────────────────────────────────────────────

def parse_pdb(path: str) -> Tuple[dict, list]:
    """
    Parse a lignoforge PDB.

    Returns
    -------
    residues : dict[res_seq -> {"name": str, "atoms": dict[atom_name -> serial]}]
    bonds    : list[(serial_a, serial_b)]
    serial_to_res : dict[serial -> res_seq]
    serial_to_name: dict[serial -> atom_name]
    """
    residues: dict = {}
    bonds:    list = []
    serial_to_res:  dict = {}
    serial_to_name: dict = {}

    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM"):
                try:
                    serial   = int(line[6:11])
                    name     = line[12:16].strip()
                    res_name = line[17:20].strip()
                    res_seq  = int(line[22:26])
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    continue
                if res_seq not in residues:
                    residues[res_seq] = {
                        "name": res_name,
                        "atoms": {},            # atom_name → serial
                        "coords": {},           # atom_name → (x,y,z)
                        "res_seq": res_seq,
                    }
                # Handle duplicate atom names within one residue (shouldn't happen
                # but be defensive): keep the first occurrence.
                if name not in residues[res_seq]["atoms"]:
                    residues[res_seq]["atoms"][name]  = serial
                    residues[res_seq]["coords"][name] = (x, y, z)
                serial_to_res[serial]  = res_seq
                serial_to_name[serial] = name

            elif rec == "CONECT":
                parts = line.split()
                src = int(parts[1])
                for tgt in parts[2:]:
                    try:
                        t = int(tgt)
                        pair = tuple(sorted((src, t)))
                        if pair not in bonds:
                            bonds.append(pair)
                    except ValueError:
                        pass

    return residues, bonds, serial_to_res, serial_to_name


def find_inter_residue_bonds(
    bonds: list,
    serial_to_res: dict,
    serial_to_name: dict,
) -> List[Tuple[int, str, int, int, str, int]]:
    """
    Return list of (serial_a, name_a, res_a, serial_b, name_b, res_b)
    for every bond that crosses a residue boundary.
    """
    inter = []
    for sa, sb in bonds:
        ra = serial_to_res.get(sa)
        rb = serial_to_res.get(sb)
        if ra is None or rb is None or ra == rb:
            continue
        na = serial_to_name[sa]
        nb = serial_to_name[sb]
        inter.append((sa, na, ra, sb, nb, rb))
    return inter


def classify_linkages(
    inter_bonds: list,
) -> List[Tuple[str, int, str, int, str]]:
    """
    Assign a (linkage_name, res_donor_seq, pos_donor, res_acceptor_seq, pos_acceptor)
    to each inter-residue bond.

    Bonds that cannot be classified are returned with linkage_name = "unknown".
    """
    result = []
    for sa, na, ra, sb, nb, rb in inter_bonds:
        pair = _sorted_pair(na, nb)
        rec  = INTER_RESIDUE_BOND_PATTERN.get(pair)
        if rec is None:
            result.append(("unknown", ra, na, rb, nb))
            continue
        lname, pos_a, pos_b = rec
        # Assign pos to the correct residue based on atom name order (a<=b).
        if na <= nb:  # na matched pair[0] → pos_a belongs to ra
            result.append((lname, ra, pos_a, rb, pos_b))
        else:
            result.append((lname, rb, pos_a, ra, pos_b))
    return result


# ── Type assignment ───────────────────────────────────────────────────────────

AtomRecord = Tuple[str, str, float, int, str]
# (opls_code, gromacs_type, charge, cgnr, comment)


def _apply_mod(
    atoms: dict,
    mod:   dict,
) -> dict:
    """Apply one linkage modification rule to a residue atom dict (in-place copy)."""
    atoms = dict(atoms)  # shallow copy

    # 1. type changes
    for aname, (opls, gtype, q) in mod.get("type_changes", {}).items():
        if aname in atoms:
            old = atoms[aname]
            atoms[aname] = (opls, gtype, q, old[3], old[4])

    # 2. remove atoms
    for aname in mod.get("remove", []):
        atoms.pop(aname, None)

    # 3. rename atoms
    for old_name, new_name in mod.get("rename", {}).items():
        if old_name in atoms:
            atoms[new_name] = atoms.pop(old_name)

    # 4. add atoms
    for aname, record in mod.get("add", {}).items():
        atoms[aname] = record

    return atoms


def _residue_nm_code(res_name: str) -> str:
    """Return the isolated neutral monomer code (HNM/GNM/SNM) for a residue."""
    base = RES_TO_BASE.get(res_name, "G")
    return BASE_TO_NM[base]


# Map (linkage_name, position_label) → atom name in THIS residue that forms
# the inter-residue bond.
_IB_TO_ATOM: Dict[Tuple[str, str], str] = {
    ("beta-O-4",  "O4H"):           "O4H",
    ("beta-O-4",  "CB"):            "CB",
    ("5-5",       "C5"):            "C5",
    ("4-O-5",     "O4H"):           "O4H",
    ("4-O-5",     "C5"):            "C5",
    ("beta-5",    "CB"):            "CB",
    ("beta-5",    "C5"):            "C5",
    ("beta-5",    "O4H_ringclose"): "O4H",
    ("beta-5",    "CA_ringclose"):  "CA",
    ("beta-beta", "CB"):            "CB",
    ("beta-beta", "CA_ringclose"):  "CA",
    ("beta-beta", "OG_ringclose"):  "OG",
    ("alpha-O-4", "O4H"):           "O4H",
    ("alpha-O-4", "CA"):            "CA",
    ("beta-1",    "CB"):            "CB",
    ("beta-1",    "C1"):            "C1",
}

# Standard intra-residue bond pairs (both atom names must be present to be kept)
_INTRA_BONDS = [
    # aromatic ring
    ("C1","C2"), ("C2","C3"), ("C3","C4"), ("C4","C5"), ("C5","C6"), ("C6","C1"),
    # ring H
    ("C2","H2"), ("C3","H3"), ("C5","H5"), ("C6","H6"),
    # phenol / ether
    ("C4","O4H"), ("O4H","HO4"),
    # OMe C3
    ("C3","OM3"), ("OM3","CM3"),
    ("CM3","HM31"), ("CM3","HM32"), ("CM3","HM33"),
    # OMe C5
    ("C5","OM5"), ("OM5","CM5"),
    ("CM5","HM51"), ("CM5","HM52"), ("CM5","HM53"),
    # C1→side-chain
    ("C1","CA"),
    # side chain
    ("CA","HA"), ("CA","OA"), ("OA","HOA"),
    ("CA","CB"), ("CB","HB"),
    ("CB","CG"), ("CG","HG1"), ("CG","HG2"), ("CG","OG"), ("OG","HOG"),
]

# Standard aromatic impropers (ring planarity)
_RING_IMPROPERS_H  = [("C2","C6","C1","CA"), ("C1","C3","C2","H2"),
                      ("C2","C4","C3","H3"), ("C3","C5","C4","O4H"),
                      ("C4","C6","C5","H5"), ("C5","C1","C6","H6")]
_RING_IMPROPERS_G  = [("C2","C6","C1","CA"), ("C1","C3","C2","H2"),
                      ("C2","C4","C3","OM3"), ("C3","C5","C4","O4H"),
                      ("C4","C6","C5","H5"), ("C5","C1","C6","H6")]
_RING_IMPROPERS_S  = [("C2","C6","C1","CA"), ("C1","C3","C2","H2"),
                      ("C2","C4","C3","OM3"), ("C3","C5","C4","O4H"),
                      ("C4","C6","C5","OM5"), ("C5","C1","C6","H6")]
_RING_IMPROPERS = {"H": _RING_IMPROPERS_H, "G": _RING_IMPROPERS_G, "S": _RING_IMPROPERS_S}


def _type_residue_atoms(
    rinfo: dict,
    linkage_context: Dict[str, Tuple[str, str]],
    base: str,
) -> dict:
    """
    Directly assign OPLS-AA types to every atom in one residue.

    The lignoforge PDB writer uses:
    - ``HB``  (singular) for the 1H on CB (vinyl or polymer sp3)
    - ``OA``, ``HOA`` only in polymer units where the α-OH is present
      (β-O-4 acceptor adds OA+HOA; other linkages do not)
    - ``HO4`` absent when C4-O4H has been converted to an aryl ether

    Types are assigned entirely from:
      1. Which atom names are present in the PDB (``rinfo["atoms"]``)
      2. The linkage context (``linkage_context``: atom → (lname, pos))

    Parameters
    ----------
    rinfo   : one entry from parse_pdb residues dict
    linkage_context : {atom_name: (linkage_name, position)} for this residue
    base    : "H", "G", or "S"
    """
    pdb   = set(rinfo["atoms"].keys())
    ic    = linkage_context  # {atom_name: (lname, pos)}
    typed: dict = {}

    # ── Aromatic ring ─────────────────────────────────────────────────────────
    if "C1" in pdb:
        typed["C1"] = ("opls_145", "CA",  0.000, 1, "Cipso (adjusted)")

    if "C2" in pdb:
        typed["C2"] = ("opls_145", "CA", -0.115, 2, "aromatic CH")
    if "H2" in pdb:
        typed["H2"] = ("opls_146", "HA",  0.115, 2, "aromatic H")

    if "C3" in pdb:
        if "OM3" in pdb:
            typed["C3"] = ("opls_199", "CA",  0.085, 3, "C3 aryl ether (OMe)")
        else:
            typed["C3"] = ("opls_145", "CA", -0.115, 3, "aromatic CH")
    if "H3" in pdb:
        typed["H3"] = ("opls_146", "HA",  0.115, 3, "aromatic H")
    if "OM3" in pdb:
        typed["OM3"]  = ("opls_179", "OS", -0.285, 3, "OMe O at C3")
    if "CM3" in pdb:
        typed["CM3"]  = ("opls_181", "CT",  0.110, 3, "OMe CH3 at C3")
    for _h in ("HM31", "HM32", "HM33"):
        if _h in pdb:
            typed[_h] = ("opls_185", "HC",  0.030, 3, "OMe H at C3")

    # C4 / O4H / HO4 – detect phenol (HO4 present) vs ether (HO4 absent)
    if "C4" in pdb:
        if "HO4" in pdb:
            typed["C4"]  = ("opls_166", "CA",  0.150, 4, "C4 phenol")
        else:
            typed["C4"]  = ("opls_199", "CA",  0.085, 4, "C4 aryl ether")
    if "O4H" in pdb:
        if "HO4" in pdb:
            typed["O4H"] = ("opls_167", "OH", -0.585, 4, "phenol O")
        else:
            typed["O4H"] = ("opls_179", "OS", -0.285, 4, "aryl ether O")
    if "HO4" in pdb:
        typed["HO4"] = ("opls_168", "HO",  0.435, 4, "phenol H")

    # C5 / H5 – aromatic CH, OMe-substituted, or linkage-substituted
    if "C5" in pdb:
        if "OM5" in pdb:
            typed["C5"] = ("opls_199", "CA",  0.085, 5, "C5 aryl ether (OMe)")
        elif "H5" in pdb:
            typed["C5"] = ("opls_145", "CA", -0.115, 5, "aromatic CH")
        else:
            lk5 = ic.get("C5", (None,))[0]
            if lk5 == "4-O-5":
                typed["C5"] = ("opls_199", "CA",  0.085, 5, "C5 aryl ether (4-O-5)")
            else:
                typed["C5"] = ("opls_145", "CA",  0.000, 5, "C5 substituted aromatic")
    if "H5" in pdb:
        typed["H5"] = ("opls_146", "HA",  0.115, 5, "aromatic H")
    if "OM5" in pdb:
        typed["OM5"]  = ("opls_179", "OS", -0.285, 5, "OMe O at C5")
    if "CM5" in pdb:
        typed["CM5"]  = ("opls_181", "CT",  0.110, 5, "OMe CH3 at C5")
    for _h in ("HM51", "HM52", "HM53"):
        if _h in pdb:
            typed[_h] = ("opls_185", "HC",  0.030, 5, "OMe H at C5")

    if "C6" in pdb:
        typed["C6"] = ("opls_145", "CA", -0.115, 6, "aromatic CH")
    if "H6" in pdb:
        typed["H6"] = ("opls_146", "HA",  0.115, 6, "aromatic H")

    # ── Side chain: CA / HA / OA / HOA ───────────────────────────────────────
    # Distinguishing contexts:
    #   vinyl isolated: OA absent, no CA linkage → CA/CB are sp2 (vinyl CM)
    #   beta-O-4 acceptor: OA present → CA sp3 alcohol, CB sp3 ether
    #   alpha-O-4 / beta-5 / beta-beta ring-closure at CA: OA absent, CA linkage
    if "CA" in pdb:
        lkCA, posCA = ic.get("CA", (None, None))
        if posCA in ("CA", "CA_ringclose"):
            # CA forms the inter-residue ether bond
            typed["CA"] = ("opls_157", "CT",  0.140, 7, "Cα sp3 ether C")
            if "HA" in pdb:
                typed["HA"] = ("opls_156", "HC",  0.060, 7, "H on Cα ether")
        elif "OA" in pdb:
            # sp3 secondary alcohol (β-O-4 acceptor side)
            typed["CA"]  = ("opls_157", "CT",  0.205, 7, "Cα sp3 alc (adjusted)")
            if "HA" in pdb:
                typed["HA"]  = ("opls_156", "HC",  0.060, 7, "H on Cα")
            typed["OA"]  = ("opls_154", "OH", -0.683, 7, "α-OH O")
            if "HOA" in pdb:
                typed["HOA"] = ("opls_155", "HO",  0.418, 7, "α-OH H")
        else:
            # vinyl sp2
            typed["CA"] = ("opls_142", "CM", -0.115, 7, "Cα vinyl sp2")
            if "HA" in pdb:
                typed["HA"] = ("opls_144", "HC",  0.115, 7, "Cα vinyl H")
    # OA/HOA that escaped the block above (shouldn't normally happen)
    if "OA" in pdb and "OA" not in typed:
        typed["OA"]  = ("opls_154", "OH", -0.683, 7, "α-OH O")
    if "HOA" in pdb and "HOA" not in typed:
        typed["HOA"] = ("opls_155", "HO",  0.418, 7, "α-OH H")

    # ── Side chain: CB / HB ────────────────────────────────────────────────────
    # Context hierarchy (in order of specificity):
    #   OA present     → β-O-4 acceptor: CB is sp3 ether C
    #   CB C-C linkage → β-5 / β-β / β-1: CB is sp3 CH
    #   otherwise      → vinyl sp2
    if "CB" in pdb:
        lkCB, posCB = ic.get("CB", (None, None))
        if "OA" in pdb:
            typed["CB"] = ("opls_157", "CT",  0.140, 8, "Cβ sp3 ether (β-O-4)")
            if "HB" in pdb:
                typed["HB"] = ("opls_156", "HC",  0.060, 8, "H on Cβ ether")
        elif lkCB in ("beta-5", "beta-beta", "beta-1"):
            typed["CB"] = ("opls_137", "CT", -0.060, 8, "Cβ sp3 CH (C-C bond)")
            if "HB" in pdb:
                typed["HB"] = ("opls_140", "HC",  0.060, 8, "H on Cβ")
        else:
            typed["CB"] = ("opls_142", "CM", -0.115, 8, "Cβ vinyl sp2")
            if "HB" in pdb:
                typed["HB"] = ("opls_144", "HC",  0.115, 8, "Cβ vinyl H")

    # ── Side chain: CG / HG1 / HG2 / OG / HOG ────────────────────────────────
    if "CG" in pdb:
        typed["CG"]  = ("opls_157", "CT",  0.145, 9, "Cγ CH2-OH")
    if "HG1" in pdb:
        typed["HG1"] = ("opls_156", "HC",  0.060, 9, "H on Cγ")
    if "HG2" in pdb:
        typed["HG2"] = ("opls_156", "HC",  0.060, 9, "H on Cγ")
    if "OG" in pdb:
        lkOG, posOG = ic.get("OG", (None, None))
        if lkOG == "beta-beta":
            typed["OG"]  = ("opls_179", "OS", -0.285, 9, "Oγ ether (β-β ring)")
        else:
            typed["OG"]  = ("opls_154", "OH", -0.683, 9, "γ-OH O")
            if "HOG" in pdb:
                typed["HOG"] = ("opls_155", "HO",  0.418, 9, "γ-OH H")

    return typed


def assign_types_to_chain(
    residues: dict,
    linkage_list: list,
) -> dict:
    """
    Build a per-residue dict of typed atoms using DIRECT typing from PDB atoms.

    Instead of starting from a sp3 reference and applying diffs, each atom is
    typed purely from:
      (a) its name and what other atoms are present in the same PDB residue, and
      (b) what inter-residue linkage events involve this residue.

    This is robust to differences between the theoretical sp3 reference and the
    actual vinyl/polymer PDB atoms emitted by lignoforge.

    Parameters
    ----------
    residues     : {res_seq → residue_info} from parse_pdb
    linkage_list : output of classify_linkages

    Returns
    -------
    typed : {res_seq → {base, nm, orig_name, res_seq, atoms, bonds, impropers,
                         linkages, inter_bonds}}
    """
    # Build per-residue linkage context: {rseq: {atom_name: (lname, pos)}}
    lk_ctx: Dict[int, Dict[str, Tuple[str, str]]] = collections.defaultdict(dict)
    for lname, r_a, pos_a, r_b, pos_b in linkage_list:
        if lname == "unknown":
            continue
        atom_a = _IB_TO_ATOM.get((lname, pos_a))
        atom_b = _IB_TO_ATOM.get((lname, pos_b))
        if atom_a:
            lk_ctx[r_a][atom_a] = (lname, pos_a)
        if atom_b:
            lk_ctx[r_b][atom_b] = (lname, pos_b)

    # Also collect linkage list per residue (for fingerprinting / naming)
    lk_list: Dict[int, list] = collections.defaultdict(list)
    for lname, r_a, pos_a, r_b, pos_b in linkage_list:
        if lname == "unknown":
            continue
        lk_list[r_a].append((lname, pos_a))
        lk_list[r_b].append((lname, pos_b))

    typed: dict = {}
    for rseq, rinfo in sorted(residues.items()):
        base = RES_TO_BASE.get(rinfo["name"], "G")
        nm   = BASE_TO_NM[base]
        ctx  = dict(lk_ctx[rseq])

        atom_types = _type_residue_atoms(rinfo, ctx, base)

        # Intra-residue bonds: include pairs where both atoms were typed
        present = set(atom_types.keys())
        bonds     = [b for b in _INTRA_BONDS     if set(b).issubset(present)]
        impropers = [b for b in _RING_IMPROPERS[base] if set(b).issubset(present)]

        typed[rseq] = {
            "base":       base,
            "nm":         nm,
            "orig_name":  rinfo["name"],
            "res_seq":    rseq,
            "atoms":      atom_types,
            "bonds":      bonds,
            "impropers":  impropers,
            "linkages":   lk_list[rseq],
            "inter_bonds": [],
        }

    # Add inter-residue bond entries for RTP ± notation
    for lname, r_a, pos_a, r_b, pos_b in linkage_list:
        if lname == "unknown":
            continue
        atom_a = _IB_TO_ATOM.get((lname, pos_a))
        atom_b = _IB_TO_ATOM.get((lname, pos_b))
        if atom_a and atom_b:
            if r_a in typed:
                typed[r_a]["inter_bonds"].append((atom_a, f"+{atom_b}", r_b))
            if r_b in typed:
                typed[r_b]["inter_bonds"].append((atom_b, f"-{atom_a}", r_a))

    return typed


# ── Per-residue charge balancing ──────────────────────────────────────────────

def balance_charges(typed: dict) -> dict:
    """
    Ensure each residue has net charge = 0.000 e by adjusting C1 charge.

    The ipso C1 has no formal OPLS-AA charge constraint (set to 0.000 by
    convention for substituted benzenes).  Any small imbalance introduced by
    linkage atom removal/type-changes is absorbed here.
    """
    for rseq, tr in typed.items():
        net = round(sum(rec[2] for rec in tr["atoms"].values()), 6)
        if abs(net) > 1e-5:
            if "C1" in tr["atoms"]:
                old = tr["atoms"]["C1"]
                new_q = round(old[2] - net, 6)
                tr["atoms"]["C1"] = (old[0], old[1], new_q, old[3], old[4])
            else:
                print(f"  [warn] res {rseq}: net charge {net:+.4f} e, C1 not found for adjustment",
                      file=sys.stderr)
    return typed


# ── Residue fingerprinting & canonical naming ─────────────────────────────────

def residue_fingerprint(tr: dict) -> str:
    """
    Produce a canonical string that uniquely identifies a residue's chemistry:
    monomer base type + sorted linkage set.

    For isolated residues (no inter-residue linkages) the side-chain form is
    also encoded:  ``|vinyl`` for the propenyl monolignol form (sp2 Cα=Cβ)
    and ``|sp3`` for the saturated propane-diol reference form.
    """
    base = tr["base"]
    lks  = sorted(f"{ln}@{pos}" for ln, pos in tr["linkages"])
    if not lks:
        cb_rec   = tr["atoms"].get("CB")
        sc_form  = "vinyl" if (cb_rec and cb_rec[0] == "opls_142") else "sp3"
        return f"{base}|{sc_form}"
    return f"{base}|{'|'.join(lks)}"


_CANONICAL_NAMES: Dict[str, str] = {
    # (fingerprint) → preferred 3-char RTP name
    #
    # Isolated forms: vinyl monolignol vs sp3 saturated reference
    "G|vinyl":  "GYM",
    "H|vinyl":  "HPM",
    "S|vinyl":  "SYM",
    "G|sp3":    "GNM",
    "H|sp3":    "HNM",
    "S|sp3":    "SNM",
    # Full β-O-4 internal units (both acceptor CB and donor O4H bonds present)
    "G|beta-O-4@CB|beta-O-4@O4H": "GYU",
    "H|beta-O-4@CB|beta-O-4@O4H": "HPU",
    "S|beta-O-4@CB|beta-O-4@O4H": "SYU",
}

# Single-character suffixes for generated names (base-36 minus 'U' which is
# already used by the canonical names GYU / HPU / SYU).
# 35 possible suffixes per prefix → more than enough for any real chain.
_SUFFIX_CHARS: str = "0123456789ABCDEFGHIJKLMNOPQRSTVWXYZ"


def assign_canonical_names(typed: dict) -> dict:
    """
    Assign each residue a canonical RTP name (max 3 characters for GROMACS
    PDB compatibility).  Identical fingerprints share a name; unique
    environments get an auto-generated name, e.g. ``GY0``, ``GY1``.
    """
    fp_to_name: dict = {}
    counters:   dict = collections.defaultdict(int)

    for rseq in sorted(typed.keys()):
        tr = typed[rseq]
        fp = residue_fingerprint(tr)
        if fp in _CANONICAL_NAMES:
            tr["rtp_name"] = _CANONICAL_NAMES[fp]
        elif fp in fp_to_name:
            tr["rtp_name"] = fp_to_name[fp]
        else:
            base   = tr["base"]
            prefix = {"H": "HP", "G": "GY", "S": "SY"}[base]
            idx    = counters[prefix]
            if idx >= len(_SUFFIX_CHARS):
                raise RuntimeError(
                    f"Too many unique residue types for prefix '{prefix}' "
                    f"(max {len(_SUFFIX_CHARS)})")
            name = prefix + _SUFFIX_CHARS[idx]   # e.g. "GY0", "GY1"
            counters[prefix] += 1
            fp_to_name[fp]   = name
            tr["rtp_name"]   = name

    return typed


# ── RTP writer ────────────────────────────────────────────────────────────────

def write_rtp(typed: dict, output_path: str) -> None:
    """Write a GROMACS .rtp file for all residues in the typed dict."""
    lines = [
        "; Custom OPLS-AA topology generated by assign_chain_types.py",
        f"; Source residues: {', '.join(tr['orig_name'] for tr in typed.values())}",
        ";",
        "[ bondedtypes ]",
        "; bonds  angles  dihedrals  impropers",
        "    1       1         3          1",
        "",
    ]

    # Track which rtp_names have already been written (write each unique block once)
    written_names: set = set()
    # But we need the atom data for each unique name.  Use first-seen residue.
    name_to_tr: dict = {}
    for rseq in sorted(typed.keys()):
        tr   = typed[rseq]
        rname = tr["rtp_name"]
        if rname not in name_to_tr:
            name_to_tr[rname] = tr

    for rname, tr in sorted(name_to_tr.items()):
        if rname in written_names:
            continue
        written_names.add(rname)

        fp  = residue_fingerprint(tr)
        lines.append(f"; ── {rname}  ({tr['base']}-type) ── fp: {fp}")
        lines.append(f"[ {rname} ]")
        lines.append(" [ atoms ]")
        lines.append(";  name   type          charge   cgnr")

        for aname, rec in tr["atoms"].items():
            opls, _gtype, q, cgnr, comment = rec
            lines.append(f"   {aname:<6s} {opls:<12s} {q:+8.4f}  {cgnr:3d}   ; {comment}")

        net = round(sum(rec[2] for rec in tr["atoms"].values()), 4)
        lines.append(f";  net charge: {net:+.4f} e")
        lines.append("")

        lines.append(" [ bonds ]")
        for b in tr["bonds"]:
            lines.append(f"   {b[0]:<6s} {b[1]}")
        for ib in tr.get("inter_bonds", []):
            lines.append(f"   {ib[0]:<6s} {ib[1]}   ; inter-residue (→ res {ib[2]})")
        lines.append("")

        if tr["impropers"]:
            lines.append(" [ impropers ]")
            for imp in tr["impropers"]:
                lines.append(f"   {'  '.join(f'{i:<6s}' for i in imp)} improper_Z_CA_X_Y")
            lines.append("")

        lines.append("")

    os.makedirs(os.path.dirname(os.path.abspath(output_path)) or ".", exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"[assign_chain_types] Written: {output_path}")


# ── PDB residue renamer ──────────────────────────────────────────────────────

def write_renamed_pdb(original_pdb: str, typed: dict, output_pdb: str) -> None:
    """
    Write a clean PDB with:
    - residue names replaced by the canonical RTP name (max 3 chars)
    - ATOM/HETATM records sorted by (resseq, original_serial) so that all
      atoms of each residue appear together in sequence
    - serials renumbered 1..N after sorting
    - CONECT records updated to the new serial numbers

    The resulting PDB can be passed directly to ``pdb2gmx`` together with
    the custom RTP and ``residuetypes_lignin.dat``.
    """
    res_to_name = {rseq: tr["rtp_name"] for rseq, tr in typed.items()}

    with open(original_pdb) as fh:
        raw_lines = fh.readlines()

    # ── Separate coordinate lines from everything else ────────────────────────
    coord_lines   = []   # (resseq, orig_serial, line)
    other_lines   = []   # header / REMARK / CONECT / END lines
    conect_lines  = []   # CONECT lines (kept separately for serial remapping)

    for line in raw_lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM") and len(line) >= 26:
            try:
                orig_serial = int(line[6:11])
                resseq      = int(line[22:26])
                rname       = res_to_name.get(resseq, line[17:20].strip())
                # Patch residue name (cols 17-20, 3 chars left-justified)
                patched = line[:17] + f"{rname:<3s}" + line[20:]
                coord_lines.append((resseq, orig_serial, patched))
            except ValueError:
                other_lines.append(line)
        elif rec == "CONECT":
            conect_lines.append(line)
        elif rec == "END":
            pass  # will be added at the end
        else:
            other_lines.append(line)

    # ── Sort coordinate records by (resseq, orig_serial) ─────────────────────
    coord_lines.sort(key=lambda t: (t[0], t[1]))

    # ── Renumber serials 1..N and build old→new mapping ───────────────────────
    old_to_new: dict = {}
    renum_coord: list = []
    for new_serial, (resseq, orig_serial, line) in enumerate(coord_lines, start=1):
        old_to_new[orig_serial] = new_serial
        # Replace serial field (cols 6-11, right-justified in 5 chars)
        patched = line[:6] + f"{new_serial:5d}" + line[11:]
        renum_coord.append(patched)

    # ── Rebuild CONECT lines with new serials ─────────────────────────────────
    new_conect: list = []
    for line in conect_lines:
        parts = line.split()
        try:
            new_parts = [parts[0]] + [str(old_to_new[int(s)]) for s in parts[1:]
                                       if int(s) in old_to_new]
            new_conect.append(f"{'CONECT':<6s}" +
                              "".join(f"{int(p):5d}" for p in new_parts[1:]) + "\n")
        except (ValueError, KeyError):
            new_conect.append(line)   # keep as-is if remapping fails

    # ── Assemble final output ─────────────────────────────────────────────────
    out = other_lines + renum_coord + new_conect + ["END\n"]

    os.makedirs(os.path.dirname(os.path.abspath(output_pdb)) or ".", exist_ok=True)
    with open(output_pdb, "w") as fh:
        fh.writelines(out)
    print(f"[assign_chain_types] Renamed PDB: {output_pdb}")


# ── High-level entry point ────────────────────────────────────────────────────

def assign_chain_types(
    pdb_path: str,
    output_rtp: Optional[str] = None,
    output_pdb: Optional[str] = None,
) -> dict:
    """
    Full pipeline: PDB → typed residues → custom RTP + renamed PDB.

    Both output files are auto-named from the input stem if not given:
    - ``<stem>_custom.rtp``
    - ``<stem>_renamed.pdb``

    Returns the typed residue dict (useful for testing / inspection).
    """
    stem = os.path.splitext(os.path.basename(pdb_path))[0]
    indir = os.path.dirname(pdb_path) or "."
    if output_rtp is None:
        output_rtp = os.path.join(indir, f"{stem}_custom.rtp")
    if output_pdb is None:
        output_pdb = os.path.join(indir, f"{stem}_renamed.pdb")

    print(f"[assign_chain_types] Reading {pdb_path}")
    residues, bonds, s2r, s2n = parse_pdb(pdb_path)
    print(f"  {len(residues)} residues, {len(bonds)} bonds")

    inter = find_inter_residue_bonds(bonds, s2r, s2n)
    print(f"  {len(inter)} inter-residue bonds")

    linkages = classify_linkages(inter)
    for ln, ra, pa, rb, pb in linkages:
        print(f"  linkage: {ln:12s}  res{ra}({pa}) — res{rb}({pb})")

    typed = assign_types_to_chain(residues, linkages)
    typed = balance_charges(typed)
    typed = assign_canonical_names(typed)

    write_rtp(typed, output_rtp)
    write_renamed_pdb(pdb_path, typed, output_pdb)
    return typed


# ── CLI ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Assign OPLS-AA atom types to a lignin chain PDB "
                    "and write a custom .rtp file + residue-renamed PDB."
    )
    parser.add_argument("pdb", help="Input PDB (lignoforge format, RTP-compatible atom names)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output .rtp path (default: <pdb_stem>_custom.rtp)")
    parser.add_argument("-p", "--output-pdb", default=None,
                        help="Output renamed PDB path (default: <pdb_stem>_renamed.pdb)")
    args = parser.parse_args(argv)
    assign_chain_types(args.pdb, args.output, args.output_pdb)


if __name__ == "__main__":
    main()
