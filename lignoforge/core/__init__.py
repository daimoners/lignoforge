"""
Core sub-package: graph-level molecular building blocks.

Provides monomer and polymer graph construction, chemical linkage rules,
structural characterisation, and graph utility functions.
"""

from lignoforge.core.rules import (
    monomer_types,
    MONOMER_RESIDUE_CODE,
    linkage_names,
    linkage_ring,
    default_color,
    CHO,
    weight_CHO,
    monomer_select_C1_C2,
    linkage_index_select_monomer,
    linkage_name_select_C1_C2,
    linkage_index_to_name,
    linkage_special_names,
)
from lignoforge.core.monomer import Monomer
from lignoforge.core.polymer import Polymer
from lignoforge.core.characterization import Characterize

__all__ = [
    "monomer_types",
    "MONOMER_RESIDUE_CODE",
    "linkage_names",
    "linkage_ring",
    "default_color",
    "Monomer",
    "Polymer",
    "Characterize",
]
