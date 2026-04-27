"""
Lignin monomer graph construction.

Each monomer is represented as a NetworkX graph where nodes are
heavy atoms (C, O) and edges are covalent bonds.  Node attributes:
    element  : str   – atomic symbol
    aromatic : bool  – part of the phenyl ring
    group    : str   – functional group label (None, 'OCH3', '4OH', '9OH')
    index    : int   – 1-based chemical carbon numbering
    mtype    : str   – monomer type ('H', 'G', 'S')
    bonding  : bool  – available for inter-monomer linkage
    color    : str   – visualisation colour
    mi       : int   – monomer index inside a polymer
"""

from __future__ import annotations

from typing import Optional
import networkx as nx

from lignoforge.core.rules import monomer_types, default_color
from lignoforge.core.utils import nxgraph, make_multi_available


# ── Low-level graph builders ──────────────────────────────────────────────────

def _add_OCH3(
    G:             nxgraph,
    monomer_type:  str,
    C_index:       int,
    color:         str,
    monomer_index: int = 0,
) -> nxgraph:
    """Attach a -OCH₃ group to carbon *C_index* (0-based graph index)."""
    if monomer_type not in monomer_types:
        raise ValueError("Monomer type must be 'H', 'G' or 'S'.")
    O_idx = len(G)
    C_idx = len(G) + 1
    G.add_node(O_idx, element="O", aromatic=False, group="OCH3",
               index=O_idx, mtype=monomer_type, bonding=False,
               color=color, mi=monomer_index)
    G.add_node(C_idx, element="C", aromatic=False, group="OCH3",
               index=C_idx, mtype=monomer_type, bonding=False,
               color=color, mi=monomer_index)
    G.add_edges_from(
        [(C_index, O_idx)], order=1, index=(C_index, O_idx),
        mtype=(monomer_type, monomer_type), btype=None,
    )
    G.add_edges_from(
        [(O_idx, C_idx)], order=1, index=(O_idx, C_idx),
        mtype=(monomer_type, monomer_type), btype=None,
    )
    return G


def _monomer_graph(
    monomer_type:  str,
    color:         str,
    monomer_index: int = 0,
) -> nxgraph:
    """Construct the base 9-carbon phenylpropane skeleton."""
    if monomer_type not in monomer_types:
        raise ValueError("Monomer type must be 'H', 'G' or 'S'.")
    G = nx.Graph()
    n_c = 9
    for i in range(1, n_c + 1):
        aromatic = i in range(1, 7)
        G.add_node(
            i - 1,
            element="C", aromatic=aromatic, group=None,
            index=i, mtype=monomer_type, bonding=False,
            color=color, mi=monomer_index,
        )
    # O at C4 (index 10 → 9 in 0-based) – the phenolic -OH
    G.add_node(9,  element="O", aromatic=False, group="4OH",
               index=10, mtype=monomer_type, bonding=False,
               color=color, mi=monomer_index)
    # O at C9 (index 11 → 10 in 0-based) – the γ-OH
    G.add_node(10, element="O", aromatic=False, group="9OH",
               index=11, mtype=monomer_type, bonding=False,
               color=color, mi=monomer_index)

    G = make_multi_available(G, monomer_type)

    # Single bonds (1-based indices converted to 0-based below)
    single = [(1,2),(2,3),(3,4),(4,5),(5,6),(6,1),(1,7),(7,8),(8,9),(10,4),(9,11)]
    for e in single:
        ei = (e[0]-1, e[1]-1)
        G.add_edges_from([ei], order=1, index=e,
                         mtype=(monomer_type, monomer_type), btype=None)
    # Double bond α-β  (C7-C8, i.e. 6-7 in 0-based)
    for e in [(7, 8)]:
        ei = (e[0]-1, e[1]-1)
        G.add_edges_from([ei], order=2, index=e,
                         mtype=(monomer_type, monomer_type), btype=None)
    return G


# ── Monomer type factories ────────────────────────────────────────────────────

def monomer_H(
    color: Optional[str] = default_color["H"],
    monomer_index: int = 0,
) -> nxgraph:
    """Build an H (p-hydroxyphenyl) monomer – no methoxy groups."""
    return _monomer_graph("H", color, monomer_index)


def monomer_G(
    color: Optional[str] = default_color["G"],
    monomer_index: int = 0,
) -> nxgraph:
    """Build a G (guaiacyl) monomer – one methoxy group at C3."""
    G = _monomer_graph("G", color, monomer_index)
    return _add_OCH3(G, "G", 2, color, monomer_index)  # C3 → 0-based index 2


def monomer_S(
    color: Optional[str] = default_color["S"],
    monomer_index: int = 0,
) -> nxgraph:
    """Build an S (syringyl) monomer – two methoxy groups at C3 and C5."""
    G = _monomer_graph("S", color, monomer_index)
    G = _add_OCH3(G, "S", 2, color, monomer_index)  # C3
    G = _add_OCH3(G, "S", 4, color, monomer_index)  # C5
    return G


# ── Monomer class ─────────────────────────────────────────────────────────────

class Monomer:
    """
    A single lignin monomer with a graph representation and a bigG node.

    Attributes
    ----------
    type  : str       – 'H', 'G' or 'S'
    G     : nxgraph   – atomistic graph (nodes are atoms)
    bigG  : nxgraph   – coarse-grained graph (one node per monomer)
    mi    : int       – monomer index inside a polymer
    """

    def __init__(self, monomer_type: str, monomer_index: int = 0):
        if monomer_type not in monomer_types:
            raise ValueError(f"Unknown monomer type '{monomer_type}'. Must be H, G or S.")
        self.type  = monomer_type
        self.mi    = monomer_index
        self.G     = None
        self.bigG  = None

    def create(self) -> nxgraph:
        """Build and store atomistic + coarse-grained graphs."""
        builders = {"H": monomer_H, "G": monomer_G, "S": monomer_S}
        self.G = builders[self.type](monomer_index=self.mi)

        big = nx.Graph()
        big.add_node(self.mi, mtype=self.type, color=default_color[self.type])
        self.bigG = big
        return self.G

    # Allow Monomer objects to be printed informatively
    def __repr__(self) -> str:
        return f"Monomer(type={self.type!r}, mi={self.mi})"
