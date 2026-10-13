"""
Virtual chemical derivatisation of lignin chains.

Mimics the analytical/industrial modification of lignin hydroxyl groups
(methylation, acetylation, longer-chain esterification) directly on the polymer
graph, before 3-D generation, so that the modified molecules go through the same
topology, force-field and simulation pipeline.

Sites (free hydroxyls present in the LignoForge monomer model):

* ``phenolic``  – ring O4H (C4–OH) that is not already an ether
* ``benzylic``  – α-OH added to Cα by β-O-4 linkage formation
* ``gamma``     – γ-OH (Cγ–OH) that is not part of a β-β ring ether

Reagents:

* ``"methyl"``       – O–CH3 ether
* ``"acyl:<n>"``     – ester O–C(=O)–(CH2)n-1–CH3 with an acyl R group of *n*
  carbons (``acyl:1`` = acetyl, ``acyl:2`` = propionyl, ``acyl:3`` = butyryl, …)
* ``"acetyl"``       – alias of ``"acyl:1"``

The model has no carboxylic acids or carbonyl groups, so reactions of those
groups (e.g. COOH esterification, C=O reduction) are not available.

Atom naming of the added atoms (unique per residue, ≤ 4 characters; ``S`` is the
site letter P/A/G for phenolic/benzylic/gamma):

    methyl ether   C<S>M            H<S>M1..3
    ester          C<S>C  O<S>C     (carbonyl C and O)
    acyl chain     C<S>1 … C<S>n    H<S>k1..  (k = chain position)
"""

from __future__ import annotations

import copy
import re
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

__all__ = ["DerivatizationRule", "find_sites", "derivatize_polymer", "parse_reagent",
           "SITE_LETTER", "SITE_KINDS"]

SITE_LETTER = {"phenolic": "P", "benzylic": "A", "gamma": "G"}
SITE_KINDS = tuple(SITE_LETTER)
_GROUP_OF_SITE = {"phenolic": "4OH", "benzylic": "alpha_OH", "gamma": "9OH"}


class DerivatizationRule(dict):
    """``{"site": ..., "reagent": ..., "degree": 0..1}`` (a plain dict subclass for JSON)."""


def parse_reagent(reagent: str) -> Tuple[str, int]:
    """Return ``("methyl", 0)`` or ``("acyl", n_carbons_in_R)``."""
    r = reagent.strip().lower()
    if r == "methyl":
        return "methyl", 0
    if r == "acetyl":
        return "acyl", 1
    m = re.fullmatch(r"acyl:(\d+)", r)
    if m and 1 <= int(m.group(1)) <= 12:
        return "acyl", int(m.group(1))
    raise ValueError(f"Unknown reagent '{reagent}' (use 'methyl', 'acetyl' or 'acyl:<1-12>')")


def find_sites(G, kind: str) -> List[int]:
    """Graph nodes of free hydroxyl oxygens of the given site kind."""
    if kind not in SITE_KINDS:
        raise ValueError(f"Unknown site '{kind}' (use {SITE_KINDS})")
    group = _GROUP_OF_SITE[kind]
    return sorted(
        n for n, a in G.nodes(data=True)
        if a.get("element") == "O" and a.get("group") == group and G.degree(n) == 1
    )


def _new_node(G, element: str, parent: int, rtp_name: str) -> int:
    pa = G.nodes[parent]
    nid = max(G.nodes) + 1
    G.add_node(nid, element=element, aromatic=False, group="derivatization",
               index=0, mtype=pa["mtype"], bonding=False, color=pa.get("color"),
               mi=pa["mi"], rtp_name=rtp_name)
    return nid


def _bond(G, a: int, b: int, order: int = 1) -> None:
    G.add_edge(a, b, order=order, index=(a, b),
               mtype=(G.nodes[a]["mtype"], G.nodes[b]["mtype"]), btype=None)


def _attach(G, o_node: int, letter: str, reagent: str, n_chain: int) -> None:
    if reagent == "methyl":
        c = _new_node(G, "C", o_node, f"C{letter}M")
        _bond(G, o_node, c)
        return
    cc = _new_node(G, "C", o_node, f"C{letter}C")
    _bond(G, o_node, cc)
    oc = _new_node(G, "O", o_node, f"O{letter}C")
    _bond(G, cc, oc, order=2)
    prev = cc
    for k in range(1, n_chain + 1):
        ck = _new_node(G, "C", o_node, f"C{letter}{k}")
        _bond(G, prev, ck)
        prev = ck


def derivatize_polymer(
    polymer,
    rules: Iterable[dict],
    seed: int = 0,
) -> Tuple[object, dict]:
    """
    Return ``(new_polymer, report)`` with the hydroxyls of *polymer* modified.

    Parameters
    ----------
    rules
        Iterable of ``{"site", "reagent", "degree"}``.  ``degree`` is the
        fraction of the *available* free sites of that kind that is modified
        (the number reacted is ``round(degree × n_sites)``, chosen at random).
        Each site kind may appear in at most one rule.
    seed
        Random seed for the choice of reacted sites.

    The input polymer is not modified.  ``report`` lists, per rule, the number
    of available and reacted sites, the achieved degree of substitution and the
    graph nodes that reacted.
    """
    rules = [dict(r) for r in rules]
    kinds = [r["site"] for r in rules]
    if len(set(kinds)) != len(kinds):
        raise ValueError("Each site kind may appear in only one derivatisation rule")
    rng = np.random.RandomState(seed)
    new = copy.deepcopy(polymer)
    G = new.G
    report: Dict[str, dict] = {"rules": []}

    for r in rules:
        site, reagent_s = r["site"], r["reagent"]
        degree = float(r.get("degree", 1.0))
        if not 0.0 <= degree <= 1.0:
            raise ValueError("degree must be between 0 and 1")
        reagent, n_chain = parse_reagent(reagent_s)
        sites = find_sites(G, site)
        n_react = int(round(degree * len(sites)))
        chosen = sorted(rng.choice(sites, size=n_react, replace=False).tolist()) if n_react else []
        letter = SITE_LETTER[site]
        for o in chosen:
            _attach(G, o, letter, reagent, n_chain)
        report["rules"].append({
            "site": site, "reagent": reagent_s, "requested_degree": degree,
            "available_sites": len(sites), "reacted_sites": len(chosen),
            "degree_of_substitution": (len(chosen) / len(sites)) if sites else 0.0,
            "nodes": chosen,
        })
    return new, report
