"""
Chemical rules for lignin monomers and inter-monomer linkages.

Notation used throughout the package:
    bond     – covalent bond between two atoms
    linkage  – inter-monomer bond connecting two monomers
    C1       – first bonding carbon in a linkage (on the existing polymer end)
    C2       – second bonding carbon in a linkage (on the incoming monomer)

Carbon index mapping (1-based chemical numbering → 0-based graph index):
    aromatic ring: C1-C6  (indices 0-5)
    alpha:         C7     (index 6)
    beta:          C8     (index 7)
    gamma:         C9     (index 8)

Linkage types supported:
    C-O:  4-O-5, alpha-O-4, beta-O-4
    C-C:  5-5, beta-5, beta-beta, beta-1
"""

# ── Atom set ──────────────────────────────────────────────────────────────────
CHO = ["C", "H", "O"]

# Atomic weights (g/mol)
weight_CHO = {"C": 12.011, "H": 1.008, "O": 15.999}

# ── Monomer types ─────────────────────────────────────────────────────────────
monomer_types = ["H", "G", "S"]

# 3-letter residue codes for output (JSON, PDB, CSV)
# H = p-Hydroxyphenyl Unit, G = Guaiacyl Unit, S = Syringyl Unit
MONOMER_RESIDUE_CODE = {"H": "HPU", "G": "GYU", "S": "SYU"}

# Default node colours for visualisation
default_color = {"H": "lightcoral", "G": "lightgreen", "S": "lightblue"}

# Linkage types that form rings (aromatic-O or phenylcoumaran / resinol)
linkage_ring = ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta"]

# ── Include beta-1 linkage ────────────────────────────────────────────────────
include_beta_1 = True

if include_beta_1:
    linkage_names = [
        "4-O-5", "alpha-O-4", "beta-O-4",
        "5-5", "beta-5", "beta-beta", "beta-1",
    ]

    # Available bonding C1 positions per monomer type,
    # mapped to the allowed C2 positions on a new (incoming) monomer.
    monomer_select_C1_C2 = {
        "H": {1: [8],  4: [5, 7, 8], 5: [4, 5, 8], 7: [4], 8: [1, 4, 5, 8]},
        "G": {1: [8],  4: [5, 7, 8], 5: [4, 5, 8], 7: [4], 8: [1, 4, 5, 8]},
        "S": {1: [8],  4: [5, 7, 8],                7: [4], 8: [1, 4, 5, 8]},
    }

    # Allowed incoming monomer types for each (C1, C2) pair
    linkage_index_select_monomer = {
        (4, 5): ["H", "G"],
        (4, 7): ["H", "G", "S"],
        (4, 8): ["H", "G", "S"],
        (5, 4): ["H", "G", "S"],
        (5, 5): ["H", "G"],
        (5, 8): ["H", "G", "S"],
        (7, 4): ["H", "G", "S"],
        (8, 4): ["H", "G", "S"],
        (8, 5): ["H", "G"],
        (8, 8): ["H", "G", "S"],
        (1, 8): ["H", "G", "S"],
        (8, 1): ["H", "G", "S"],
    }

    linkage_name_select_C1_C2 = {
        "4-O-5":     {4: 5, 5: 4},
        "alpha-O-4": {4: 7, 7: 4},
        "beta-O-4":  {4: 8, 8: 4},
        "5-5":       {5: 5},
        "beta-5":    {8: 5, 5: 8},
        "beta-beta": {8: 8},
        "beta-1":    {1: 8, 8: 1},
    }

    linkage_index_to_name = {
        (4, 5): "4-O-5",
        (4, 7): "alpha-O-4",
        (4, 8): "beta-O-4",
        (5, 4): "4-O-5",
        (5, 5): "5-5",
        (5, 8): "beta-5",
        (7, 4): "alpha-O-4",
        (8, 4): "beta-O-4",
        (8, 5): "beta-5",
        (8, 8): "beta-beta",
        (1, 8): "beta-1",
        (8, 1): "beta-1",
    }

    # Linkages needing special atom-level handling
    # (oxygen bridging, ring closure, alpha-OH removal for beta-1)
    linkage_special_names = [
        "4-O-5", "alpha-O-4", "beta-O-4",
        "beta-5", "beta-beta", "beta-1",
    ]

    # Valid (monomer1_type, monomer2_type, linkage_name) combinations
    monomer1_select_monomer2_linkage_name = {
        "H": {
            "H": ["alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta", "beta-1"],
            "G": ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta", "beta-1"],
            "S": ["alpha-O-4", "beta-O-4", "beta-beta", "beta-1"],
        },
        "G": {
            "H": ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta", "beta-1"],
            "G": ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta", "beta-1"],
            "S": ["alpha-O-4", "beta-O-4", "beta-beta", "beta-1"],
        },
        "S": {
            "H": ["alpha-O-4", "beta-O-4", "beta-beta", "beta-1"],
            "G": ["alpha-O-4", "beta-O-4", "beta-beta", "beta-1"],
            "S": ["alpha-O-4", "beta-O-4", "beta-beta", "beta-1"],
        },
    }

else:
    # ── Without beta-1 ────────────────────────────────────────────────────────
    linkage_names = ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta"]

    monomer_select_C1_C2 = {
        "H": {4: [5, 7, 8], 5: [4, 5, 8], 7: [4], 8: [4, 5, 8]},
        "G": {4: [5, 7, 8], 5: [4, 5, 8], 7: [4], 8: [4, 5, 8]},
        "S": {4: [5, 7, 8],                7: [4], 8: [4, 5, 8]},
    }

    linkage_index_select_monomer = {
        (4, 5): ["H", "G"],
        (4, 7): ["H", "G", "S"],
        (4, 8): ["H", "G", "S"],
        (5, 4): ["H", "G", "S"],
        (5, 5): ["H", "G"],
        (5, 8): ["H", "G", "S"],
        (7, 4): ["H", "G", "S"],
        (8, 4): ["H", "G", "S"],
        (8, 5): ["H", "G"],
        (8, 8): ["H", "G", "S"],
    }

    linkage_name_select_C1_C2 = {
        "4-O-5":     {4: 5, 5: 4},
        "alpha-O-4": {4: 7, 7: 4},
        "beta-O-4":  {4: 8, 8: 4},
        "5-5":       {5: 5},
        "beta-5":    {8: 5, 5: 8},
        "beta-beta": {8: 8},
    }

    linkage_index_to_name = {
        (4, 5): "4-O-5",
        (4, 7): "alpha-O-4",
        (4, 8): "beta-O-4",
        (5, 4): "4-O-5",
        (5, 5): "5-5",
        (5, 8): "beta-5",
        (7, 4): "alpha-O-4",
        (8, 4): "beta-O-4",
        (8, 5): "beta-5",
        (8, 8): "beta-beta",
    }

    linkage_special_names = ["4-O-5", "alpha-O-4", "beta-O-4", "beta-5", "beta-beta"]

    monomer1_select_monomer2_linkage_name = {
        "H": {
            "H": ["alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta"],
            "G": ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta"],
            "S": ["alpha-O-4", "beta-O-4", "beta-beta"],
        },
        "G": {
            "H": ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta"],
            "G": ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta"],
            "S": ["alpha-O-4", "beta-O-4", "beta-beta"],
        },
        "S": {
            "H": ["alpha-O-4", "beta-O-4", "beta-beta"],
            "G": ["alpha-O-4", "beta-O-4", "beta-beta"],
            "S": ["alpha-O-4", "beta-O-4", "beta-beta"],
        },
    }
