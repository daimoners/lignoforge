"""
Force-field sub-package: OPLS-AA typing, charge neutralisation and graph-based
GROMACS topology generation for lignin chains.
"""

from lignoforge.forcefield.opls import TypingError
from lignoforge.forcefield.topology import (
    Atom,
    ChainTopology,
    build_chain_topology,
    write_gromacs_system,
)

__all__ = [
    "Atom",
    "ChainTopology",
    "TypingError",
    "build_chain_topology",
    "write_gromacs_system",
]
