"""
Polymer characterisation: monomer counts, linkage distributions,
branching coefficient, functional groups, and molecular weight.
"""

from __future__ import annotations

import os
import warnings
from collections import Counter
from typing import Optional

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", message="divide by zero")
warnings.filterwarnings("ignore", message="invalid value encountered")

from lignoforge.core.rules import linkage_names, monomer_types
from lignoforge.core.utils import (
    formula_to_MW,
    graph_to_smile,
    smiles_to_formula,
    nxgraph,
    nparray,
)
from lignoforge.core.polymer import Polymer


# ── Monomer-atom counts per type ───────────────────────────────────────────────
# H: 11 heavy atoms  |  G: 13  |  S: 15
_ATOMS_PER_MONOMER = {"H": 11, "G": 13, "S": 15}


# ── Convenience function ───────────────────────────────────────────────────────

def get_metrics_polymer(
    P: Polymer,
    additional: bool = False,
    cal_MW:     bool = False,
) -> tuple[nparray, int, Optional[float]]:
    """Return (metrics_array, monomer_count, MW) for a Polymer."""
    ch = Characterize(P)
    ch.get_metrics(cal_MW=cal_MW, additional=additional)
    MW = ch.MW if cal_MW else None
    return ch.metrics, ch.monomer_count, MW


def get_counts_polymer(
    P: Polymer,
    additional: bool = False,
    cal_MW:     bool = False,
) -> tuple[nparray, int, Optional[float]]:
    """Return (counts_array, monomer_count, MW) for a Polymer."""
    ch = Characterize(P)
    ch.get_counts(cal_MW=cal_MW, additional=additional)
    MW = ch.MW if cal_MW else None
    return ch.counts, ch.monomer_count, MW


# ── CharacterizeGraph: works on a raw NetworkX graph ─────────────────────────

class CharacterizeGraph:
    """Characterise a raw atomistic graph (used by CharacterizeGraph)."""

    def __init__(self, G: nxgraph):
        self.G = G
        self.mtype_count:   Optional[dict]  = None
        self.monomer_count: Optional[float] = None
        self.linkages_count:Optional[dict]  = None
        self.OCH3_count:    Optional[float] = None
        self.OH_count:      Optional[float] = None
        self.MW:            Optional[float] = None
        self.smiles:        Optional[str]   = None
        self.formula:       Optional[str]   = None
        self.metrics:       Optional[nparray] = None
        self.counts:        Optional[nparray] = None

    def count_types(self) -> dict:
        mtypes = [self.G.nodes[n]["mtype"] for n in self.G.nodes]
        self.mtype_count = {
            t: mtypes.count(t) / _ATOMS_PER_MONOMER[t]
            for t in monomer_types
        }
        return self.mtype_count

    def count_monomers(self) -> float:
        self.monomer_count = sum(self.mtype_count.values())
        return self.monomer_count

    def count_linkages(self) -> dict:
        bonds = [self.G.edges[e]["btype"] for e in self.G.edges]
        lc    = {name: float(bonds.count(name)) for name in linkage_names}
        # Correct for atoms that appear multiple times per linkage
        lc["4-O-5"]    /= 2
        lc["beta-O-4"] /= 2
        lc["beta-5"]   /= 2
        lc["beta-beta"]/= 3
        self.linkages_count = lc
        return lc

    def count_OCH3(self) -> float:
        groups = [self.G.nodes[n]["group"] for n in self.G.nodes]
        self.OCH3_count = float(groups.count("OCH3")) / 2
        return self.OCH3_count

    def count_OH(self) -> float:
        groups = list(self.G.nodes(data="group"))
        oh4 = [ni for ni, g in groups if g == "4OH" and self.G.degree(ni) == 1]
        oh9 = [ni for ni, g in groups if g == "9OH" and self.G.degree(ni) == 1]
        self.OH_count = float(len(oh4) + len(oh9))
        return self.OH_count

    def cal_MW(self) -> float:
        self.smiles  = graph_to_smile(self.G)
        self.formula = smiles_to_formula(self.smiles)
        self.MW      = formula_to_MW(self.formula)
        return self.MW

    def cal_all(self, cal_MW: bool = False, print_flag: bool = True):
        self.count_types()
        self.count_monomers()
        self.count_linkages()
        self.count_OCH3()
        self.count_OH()
        if cal_MW:
            self.cal_MW()

    def cal_metrics(self, cal_MW: bool = False) -> nparray:
        self.cal_all(cal_MW, print_flag=False)
        m = np.array(list(self.mtype_count.values()), dtype=float)
        l = np.array(list(self.linkages_count.values()), dtype=float)
        ms = m.sum(); ls = l.sum()
        m  = m / ms if ms > 0 else np.zeros_like(m)
        l  = l / ls if ls > 0 else np.zeros_like(l)
        if np.isnan(l).any():
            l = np.zeros_like(l)
        self.metrics = np.concatenate([m, l])
        return self.metrics


# ── Characterize: uses the polymer's bigG for monomer-level counts ─────────────

class Characterize(CharacterizeGraph):
    """Full polymer characterisation using both atomistic and CG graphs."""

    def __init__(self, P: Polymer):
        super().__init__(P.G)
        self.bigG = P.bigG
        self.connections_count: Optional[dict]  = None
        self.branching_coeff:   Optional[float] = None
        self.n_branched:        Optional[int]   = None

    # ── Override count_types / count_monomers using bigG ─────────────────────

    def count_types(self) -> dict:
        mtypes = [self.bigG.nodes[n]["mtype"] for n in self.bigG.nodes]
        self.mtype_count = {t: mtypes.count(t) for t in monomer_types}
        return self.mtype_count

    def count_monomers(self) -> int:
        self.monomer_count = len(self.bigG)
        return self.monomer_count

    def count_linkages(self) -> dict:
        bonds = [self.bigG.edges[e]["btype"] for e in self.bigG.edges]
        self.linkages_count = {name: bonds.count(name) for name in linkage_names}
        return self.linkages_count

    def count_connections(self) -> dict:
        self.connections_count = dict(
            Counter(self.bigG.degree(n) for n in self.bigG.nodes)
        )
        return self.connections_count

    def cal_branching(self) -> float:
        if self.connections_count is None:
            self.count_connections()
        if self.monomer_count is None:
            self.count_monomers()
        n_branched = sum(v for k, v in self.connections_count.items() if k >= 3)
        self.n_branched     = n_branched
        self.branching_coeff = n_branched / self.monomer_count if self.monomer_count else 0.0
        return self.branching_coeff

    def cal_all(self, cal_MW: bool = False, print_flag: bool = True):
        self.count_types()
        self.count_monomers()
        self.count_linkages()
        self.count_OCH3()
        self.count_OH()
        self.count_connections()
        self.cal_branching()
        if cal_MW:
            self.cal_MW()

    def get_metrics(self, additional: bool = False, cal_MW: bool = False) -> nparray:
        self.cal_all(cal_MW, print_flag=False)
        m = np.array(list(self.mtype_count.values()), dtype=float)
        l = np.array(list(self.linkages_count.values()), dtype=float)
        ms = m.sum(); ls = l.sum()
        m  = m / ms if ms > 0 else np.zeros_like(m)
        l  = l / ls if ls > 0 else np.zeros_like(l)
        if np.isnan(l).any():
            l = np.zeros_like(l)
        self.metrics = np.concatenate([m, l])
        if additional:
            self.metrics = np.concatenate([self.metrics, [self.branching_coeff]])
        return self.metrics

    def get_counts(self, additional: bool = False, cal_MW: bool = False) -> nparray:
        self.cal_all(cal_MW, print_flag=False)
        m = np.array(list(self.mtype_count.values()), dtype=float)
        l = np.array(list(self.linkages_count.values()), dtype=float)
        self.counts = np.concatenate([m, l])
        if additional:
            self.counts = np.concatenate([self.counts, [self.branching_coeff]])
        return self.counts

    def summary(self) -> dict:
        """Return a flat summary dictionary of all computed properties."""
        self.cal_all(cal_MW=True, print_flag=False)
        return {
            "monomer_count":   self.monomer_count,
            "H_count":         self.mtype_count["H"],
            "G_count":         self.mtype_count["G"],
            "S_count":         self.mtype_count["S"],
            "branching_coeff": self.branching_coeff,
            "OCH3_count":      self.OCH3_count,
            "OH_count":        self.OH_count,
            "MW":              self.MW,
            "smiles":          self.smiles or graph_to_smile(self.G),
            **{f"linkage_{k}": v for k, v in self.linkages_count.items()},
        }


# ── Population characterisation ────────────────────────────────────────────────

class Population:
    """Characterise an ensemble (population) of Polymer objects."""

    def __init__(
        self,
        P_population: list[Polymer],
        library_name: str,
        ResultsName:  str  = "results",
        TrialIndex:   Optional[str] = None,
        save_path:    str  = None,
    ):
        self.P_population = P_population
        self.library_name = library_name
        self.ResultsName  = ResultsName
        self.TrialIndex   = TrialIndex
        self.save_path    = save_path or os.getcwd()

        if TrialIndex is not None:
            self.ResultsPath = os.path.join(
                self.save_path, ResultsName, library_name, "i" + str(TrialIndex)
            )
        else:
            self.ResultsPath = os.path.join(self.save_path, ResultsName, library_name)

        self.summaries:   Optional[list]  = None
        self.df:          Optional[pd.DataFrame] = None

    def analyze(self) -> pd.DataFrame:
        """Characterise every polymer; store results in a DataFrame."""
        rows = []
        for pi, P in enumerate(self.P_population):
            ch = Characterize(P)
            row = ch.summary()
            row["polymer_id"] = pi
            rows.append(row)
        self.summaries = rows
        self.df = pd.DataFrame(rows)
        os.makedirs(self.ResultsPath, exist_ok=True)
        csv_path = os.path.join(self.ResultsPath, f"{self.library_name}_library.csv")
        self.df.to_csv(csv_path, index=False)
        return self.df

    def get_metrics_mean(self, additional: bool = False) -> dict:
        """Return the population-mean metrics dictionary."""
        if self.df is None:
            self.analyze()
        cols = (
            [f"linkage_{n}" for n in linkage_names]
            + [f"{t}_count" for t in monomer_types]
        )
        if additional:
            cols += ["branching_coeff"]
        subset = self.df[[c for c in cols if c in self.df.columns]]
        return subset.mean().to_dict()

    def population_stats(self) -> pd.DataFrame:
        """Return per-metric mean ± std for the population."""
        if self.df is None:
            self.analyze()
        stats = self.df.describe().T[["mean", "std", "min", "max"]]
        csv_path = os.path.join(
            self.ResultsPath, f"{self.library_name}_population_stats.csv"
        )
        stats.to_csv(csv_path)
        return stats
