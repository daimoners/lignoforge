"""
Lignin polymer graph construction.

PolymerGraph  – low-level graph operations (linkage formation, ring closure).
Polymer       – high-level interface for step-wise polymer growth.
"""

from __future__ import annotations

from copy import copy
from typing import Optional, Tuple
import warnings

import networkx as nx

from lignoforge.core.rules import (
    linkage_special_names,
    monomer_select_C1_C2,
    linkage_index_select_monomer,
    linkage_index_to_name,
    linkage_name_select_C1_C2,
    linkage_ring,
)
from lignoforge.core.monomer import Monomer
from lignoforge.core.utils import (
    nxgraph,
    make_unavailable,
    adjust_indices,
    generate_random_monomer,
    generate_random_linkage,
    generate_random_branching_state,
)
from lignoforge.core.rules import monomer_types, linkage_names, default_color


ENABLE_BETA_1_LINKAGE = False


class PolymerGraph:
    """Low-level graph operations for linking monomers."""

    def __init__(self, G: nxgraph, verbose: bool = True):
        self.G    = G.copy()
        self.verbose = verbose
        self.C1_indices_in_polymer = None

    # ── C1 / C2 discovery ─────────────────────────────────────────────────────

    def find_available_C1_in_polymer(
        self, branching_state: Optional[bool] = None
    ) -> list:
        """Return active bonding-carbon indices on the existing polymer."""
        C1s = [n for n, v in self.G.nodes(data=True) if v["bonding"]]
        if branching_state is not None:
            C1s = self._filter_by_branching(C1s, branching_state)
        return C1s

    def update_available_C1_in_polymer(self, used_C1: int) -> list:
        return [n for n in self.C1_indices_in_polymer if n != used_C1]

    def _filter_by_branching(
        self, C_indices: list, branching_state: bool
    ) -> list:
        """Keep only C1 indices that belong to terminal / interior monomers."""
        if len(self.bigG) < 3:
            return C_indices
        terminal_mi = {
            mi for mi, node in enumerate(self.bigG.nodes)
            if self.bigG.degree(node) == 1
        }
        result = []
        for ci in C_indices:
            mi = self.G.nodes[ci]["mi"]
            if branching_state and mi not in terminal_mi:
                result.append(ci)
            elif (not branching_state) and mi in terminal_mi:
                result.append(ci)
        return result

    def find_available_C2_in_monomer(
        self, C1_index: int, ring: bool = False
    ) -> list:
        node = self.G.nodes[C1_index]
        C2s  = monomer_select_C1_C2[node["mtype"]][node["index"]]
        if ring:
            C2s = [c for c in C2s if c != 1]
        return C2s

    def find_available_monomer_types(self, bond_index: tuple) -> list:
        return linkage_index_select_monomer[bond_index]

    def find_C2_index_in_polymer(self, C2_in_monomer: int) -> int:
        return C2_in_monomer + len(self.G) - 1

    def find_O_index_in_polymer(self, C_index: int) -> int:
        return next(
            i for i in self.G.neighbors(C_index)
            if self.G.nodes[i]["element"] == "O"
        )

    # ── Linkage formation ──────────────────────────────────────────────────────

    def connect_C1_C2(
        self,
        linkage_index:     Tuple[int, int],
        C1_index:          int,
        C2_index:          int,
    ) -> bool:
        """
        Form a new inter-monomer bond.  Returns True if the linkage was added.
        """
        new_flag = False
        C1_node  = self.G.nodes[C1_index]
        C2_node  = self.G.nodes[C2_index]
        lname    = linkage_index_to_name[linkage_index]
        ltypes   = (C1_node["mtype"], C2_node["mtype"])

        if lname == "beta-1" and not ENABLE_BETA_1_LINKAGE:
            return False

        bond_list: list[tuple] = []

        # Nodes that must be made unavailable after bonding
        extra_unavail = []

        # ── C-O bridged linkages (4-O-5, alpha-O-4, beta-O-4) ─────────────────
        if lname in linkage_special_names:

            if 4 in linkage_index:
                O_idx = (
                    self.find_O_index_in_polymer(C1_index)
                    if linkage_index[0] == 4
                    else self.find_O_index_in_polymer(C2_index)
                )
                bond_list = [(C1_index, O_idx), (O_idx, C2_index)]

                # β-O-4: Cβ (C8) bonds to the phenol-O → the vinyl double bond
                # Cα=Cβ is broken.  Cα needs an α-OH (secondary alcohol) to
                # restore sp3 tetrahedral valence.  beta_idx is always C8; Cα
                # is always the immediately preceding node (C8-1 in graph order).
                if lname == "beta-O-4":
                    beta_idx  = C1_index if linkage_index[0] == 8 else C2_index
                    alpha_idx = beta_idx - 1
                    alpha_node = self.G.nodes.get(alpha_idx, {})
                    # Verify this node really is Cα (1-based chemical index 7)
                    # and does not already carry an oxygen neighbour.
                    if (alpha_node.get("element") == "C"
                            and alpha_node.get("index") == 7
                            and not any(
                                self.G.nodes[nb].get("element") == "O"
                                for nb in self.G.neighbors(alpha_idx)
                            )):
                        O_alpha_idx = len(self.G)
                        self.G.add_node(
                            O_alpha_idx,
                            element="O", aromatic=False, group="alpha_OH",
                            index=O_alpha_idx, mtype=alpha_node["mtype"],
                            bonding=False, color=alpha_node["color"],
                            mi=alpha_node["mi"],
                        )
                        bond_list.append((alpha_idx, O_alpha_idx))

            # beta-5: β-5 bond + O-α ring closure
            elif lname == "beta-5":
                if linkage_index[0] == 8:
                    O_idx    = self.find_O_index_in_polymer(C2_index - 1)
                    alpha_i  = C1_index - 1
                    C4_i     = C2_index - 1
                else:
                    O_idx    = self.find_O_index_in_polymer(C1_index - 1)
                    alpha_i  = C2_index - 1
                    C4_i     = C1_index - 1
                if self.G.nodes[C4_i]["bonding"] and self.G.nodes[alpha_i]["bonding"]:
                    bond_list = [(C1_index, C2_index), (O_idx, alpha_i)]
                    extra_unavail = [alpha_i, C4_i]

            # beta-beta: β-β + two α-O ring closures
            elif lname == "beta-beta":
                O1_idx  = self.find_O_index_in_polymer(C1_index + 1)
                O2_idx  = self.find_O_index_in_polymer(C2_index + 1)
                alpha1  = C1_index - 1
                alpha2  = C2_index - 1
                if self.G.nodes[alpha1]["bonding"] and self.G.nodes[alpha2]["bonding"]:
                    bond_list = [
                        (C1_index, C2_index),
                        (alpha1, O2_idx),
                        (alpha2, O1_idx),
                    ]
                    extra_unavail = [alpha1, alpha2]

            # beta-1: rearrangement (remove α atoms from C2 side, add -OH)
            elif lname == "beta-1":
                O_new_idx = len(self.G)
                if linkage_index[0] == 8:
                    alpha1_i = C1_index - 1
                    alpha2_i = C2_index - 1
                else:
                    alpha1_i = C2_index - 1
                    alpha2_i = C1_index - 1

                alpha1_node = self.G.nodes[alpha1_i]
                # Remove atoms 7-9 and 11 relative to alpha2
                to_del = [alpha2_i + d for d in [7, 8, 9, 11]]
                self.G.remove_nodes_from(to_del)

                # Add a new -OH on alpha1
                self.G.add_node(
                    O_new_idx,
                    element="O", aromatic=False, group="7OH",
                    index=16, mtype=alpha1_node["mtype"],
                    bonding=False, color=alpha1_node["color"],
                    mi=alpha1_node["mi"],
                )
                # Reset double bond α-β → single
                alpha2_new = alpha1_i + 1
                if (alpha1_i, alpha2_new) in self.G.edges:
                    self.G.edges[alpha1_i, alpha2_new]["order"] = 1

                bond_list = [
                    (C1_index, C2_index),
                    (alpha1_i, O_new_idx),
                ]
                extra_unavail = [alpha1_i]

        # ── Simple C-C bond ────────────────────────────────────────────────────
        else:
            bond_list = [(C1_index, C2_index)]

        # ── Commit bonds ───────────────────────────────────────────────────────
        if bond_list:
            self.G.add_edges_from(
                bond_list, order=1, index=linkage_index,
                mtype=ltypes, btype=lname,
            )
            self.G = make_unavailable(self.G, C1_index)
            self.G = make_unavailable(self.G, C2_index)
            for ni in extra_unavail:
                self.G = make_unavailable(self.G, ni)

            if self.verbose:
                print(f"  {ltypes[0]}-{ltypes[1]} via {lname}")
            new_flag = True

        if lname == "beta-1" and new_flag:
            self.G = adjust_indices(self.G)

        return new_flag


# ── Polymer class ─────────────────────────────────────────────────────────────

class Polymer(PolymerGraph):
    """
    Growing lignin polymer.

    Can be initialised from a Monomer or an existing Polymer (deep copy).
    """

    def __init__(self, M_init: object, verbose: bool = True):
        if M_init.G is None and isinstance(M_init, Monomer):
            M_init.create()
        super().__init__(M_init.G, verbose)
        self.bigG = M_init.bigG.copy()
        self.mi   = 0 if isinstance(M_init, Monomer) else M_init.mi

    # ── Book-keeping ───────────────────────────────────────────────────────────

    def _add_monomer_to_bigG(
        self, new_mi: int, mtype: str, C1_mi: int, lname: str
    ) -> None:
        """Register a new monomer node and inter-monomer edge in the CG graph."""
        self.bigG.add_node(new_mi, mtype=mtype, color=default_color[mtype])
        self.bigG.add_edge(C1_mi, new_mi, btype=lname)

    # ── Monomer addition ───────────────────────────────────────────────────────

    def add_specific_monomer(
        self,
        monomer_type:    str,
        linkage_type:    str,
        branching_state: Optional[bool] = None,
    ) -> bool:
        """
        Attempt to add one new monomer via a specified linkage.

        Returns True if the monomer was successfully added.
        """
        C1_list = self.find_available_C1_in_polymer(branching_state)
        if not C1_list:
            return False

        # Get the (C1, C2) pair for the requested linkage
        linkage_C1_C2 = linkage_name_select_C1_C2.get(linkage_type, {})

        for C1_idx in C1_list:
            C1_node   = self.G.nodes[C1_idx]
            C1_pos    = C1_node["index"]  # 1-based chemical carbon number
            C1_mi_val = C1_node["mi"]

            if C1_pos not in linkage_C1_C2:
                continue

            C2_pos_in_monomer = linkage_C1_C2[C1_pos]
            if monomer_type not in linkage_index_select_monomer.get(
                (C1_pos, C2_pos_in_monomer), []
            ):
                continue

            # Build new monomer and join graphs
            new_mi = self.mi + 1
            M_new  = Monomer(monomer_type, new_mi)
            M_new.create()

            G_joined = _join_graphs(self.G, M_new.G)
            self_copy = Polymer.__new__(Polymer)
            self_copy.G        = G_joined
            self_copy.bigG     = self.bigG.copy()
            self_copy.verbose  = self.verbose
            self_copy.mi       = new_mi
            self_copy.C1_indices_in_polymer = None

            C2_idx_joined = C2_pos_in_monomer + len(self.G) - 1

            success = PolymerGraph.connect_C1_C2(
                self_copy,
                (C1_pos, C2_pos_in_monomer),
                C1_idx,
                C2_idx_joined,
            )
            if success:
                self.G    = self_copy.G
                self.bigG = self_copy.bigG
                self.mi   = new_mi
                self._add_monomer_to_bigG(new_mi, monomer_type, C1_mi_val, linkage_type)
                return True

        return False

    def add_random_monomer(
        self,
        monomer_distribution: list,
        linkage_distribution: list,
        branching_state:      Optional[bool]       = None,
        random_state:         Optional[object]     = None,
    ) -> bool:
        """Add a randomly sampled monomer–linkage pair."""
        ltype = generate_random_linkage(linkage_distribution, random_state)
        mtype = generate_random_monomer(monomer_distribution, random_state)
        return self.add_specific_monomer(mtype, ltype, branching_state)

    # ── Ring closure ───────────────────────────────────────────────────────────

    def add_specific_ring(self, linkage_type: str) -> bool:
        """
        Attempt to close a ring between two existing monomers in the polymer.
        """
        if linkage_type not in linkage_ring:
            return False

        linkage_C1_C2 = linkage_name_select_C1_C2.get(linkage_type, {})

        C1_list = [
            n for n, v in self.G.nodes(data=True)
            if v["bonding"] and v["index"] in linkage_C1_C2
        ]
        if len(C1_list) < 2:
            return False

        import random
        random.shuffle(C1_list)

        for i, C1_idx in enumerate(C1_list):
            C1_node  = self.G.nodes[C1_idx]
            C1_pos   = C1_node["index"]
            C2_pos   = linkage_C1_C2[C1_pos]

            for C2_idx in C1_list[i + 1:]:
                C2_node = self.G.nodes[C2_idx]
                if C2_node["index"] != C2_pos:
                    continue
                if C2_node["mi"] == C1_node["mi"]:
                    continue  # same monomer → skip

                allowed_mtypes = linkage_index_select_monomer.get(
                    (C1_pos, C2_pos), []
                )
                if C2_node["mtype"] not in allowed_mtypes:
                    continue

                success = PolymerGraph.connect_C1_C2(
                    self,
                    (C1_pos, C2_pos),
                    C1_idx,
                    C2_idx,
                )
                if success:
                    return True
        return False

    def __repr__(self) -> str:
        return f"Polymer(monomers={self.mi + 1}, nodes={len(self.G)})"


# ── Graph join helper ──────────────────────────────────────────────────────────

def _join_graphs(G1: nxgraph, G2: nxgraph) -> nxgraph:
    """Disjoint union of two graphs preserving node attributes."""
    return nx.disjoint_union(G1, G2)
