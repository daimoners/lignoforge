"""
LigninPriorEstimator
====================

Translates high-level experimental / industrial JSON input into a complete
set of structural and statistical priors required by the kMC simulation engine.

Strategy (in order of precedence):
    1. If the input JSON already contains an experimental value → use it directly.
    2. If partial data are available → constrain the literature prior.
    3. If no data are available → sample from the literature prior (Gaussian).

Output ``priors`` dictionary keys
----------------------------------
    S_fraction, G_fraction, H_fraction   – monomer molar fractions (sum ~1)
    linkage_fractions                     – dict {linkage_name: fraction} (sum = 1)
    Mn, Mw, PDI                          – molecular weight statistics (g/mol)
    branching_index                       – fraction of monomers with ≥3 bonds
    condensation_degree                   – C-C / total inter-monomer bond ratio
    mean_DP, max_DP                      – degree of polymerisation (monomer count)
    avg_monomer_MW                        – composition-weighted monomer MW (g/mol)
"""

from __future__ import annotations

import warnings
from typing import Any, Optional

import numpy as np

from lignoforge.core.rules import linkage_names
from lignoforge.priors.literature_data import (
    SGH_BY_BIOMASS_PROCESS,
    LINKAGE_BY_BIOMASS_PROCESS,
    MW_BY_PROCESS,
    BRANCHING_BY_PROCESS,
    CONDENSATION_BY_PROCESS,
    MONOMER_MW,
)


class LigninPriorEstimator:
    """
    Estimate structural priors for a lignin sample from its experimental descriptor.

    Parameters
    ----------
    input_data : dict
        Validated experimental/industrial input dictionary
        (as defined by ``lignin_info_schema.json``).
    random_seed : int | None
        Seed for reproducibility. ``None`` → non-deterministic.
    """

    def __init__(self, input_data: dict, random_seed: Optional[int] = None):
        self.data   = input_data
        self.priors: dict[str, Any] = {}
        self._rng   = np.random.default_rng(random_seed)

    # ── Safe path-getter ───────────────────────────────────────────────────────

    def _get(self, path: str, default=None):
        """Retrieve a nested value via a dot-separated path."""
        keys = path.split(".")
        val  = self.data
        try:
            for k in keys:
                val = val[k]
            return val
        except (KeyError, TypeError):
            return default

    # ── Internal sampling helpers ─────────────────────────────────────────────

    def _sample(self, mean: float, std: float, lo: float = 0.0, hi: float = 1.0) -> float:
        """Draw a truncated-normal sample."""
        val = self._rng.normal(mean, std)
        return float(np.clip(val, lo, hi))

    def _lookup_process_key(self, table: dict, process: str, fallback: str = "kraft"):
        """Retrieve a process-keyed entry, falling back gracefully."""
        return table.get(process, table.get(fallback, next(iter(table.values()))))

    # ── 1. S/G/H fractions ────────────────────────────────────────────────────

    def infer_SGH(self) -> None:
        """Infer S, G, H molar fractions."""
        # Use directly provided values if available
        sgh = self._get("chemical_composition.S_G_H_ratio")
        if sgh and all(k in sgh for k in ("S", "G", "H")):
            S, G, H = float(sgh["S"]), float(sgh["G"]), float(sgh["H"])
            total   = S + G + H
            if total > 0:
                self.priors.update(
                    S_fraction=S / total,
                    G_fraction=G / total,
                    H_fraction=H / total,
                )
                return

        biomass = self._get("material_origin.biomass_type", "hardwood")
        process = self._get("extraction_process.process_type", "kraft")

        bio_table = SGH_BY_BIOMASS_PROCESS.get(biomass, SGH_BY_BIOMASS_PROCESS["hardwood"])
        Sm, Ss, Gm, Gs, Hm, Hs = self._lookup_process_key(bio_table, process)

        # Sample and clip to [0, 1]
        S = max(0.0, self._rng.normal(Sm, Ss))
        G = max(0.0, self._rng.normal(Gm, Gs))
        H = max(0.0, self._rng.normal(Hm, Hs))
        total = S + G + H or 1.0

        self.priors["S_fraction"] = float(S / total)
        self.priors["G_fraction"] = float(G / total)
        self.priors["H_fraction"] = float(H / total)

    # ── 2. Linkage distribution ───────────────────────────────────────────────

    def infer_linkage_distribution(self) -> None:
        """Infer the inter-monomer linkage distribution."""
        # Try to use explicitly provided bonding profile
        bp = self._get("chemical_composition.bonding_profile")
        provided: dict[str, float] = {}
        if bp:
            key_map = {
                "beta_O_4_fraction":  "beta-O-4",
                "beta_5_fraction":    "beta-5",
                "beta_beta_fraction": "beta-beta",
                "five_five_fraction": "5-5",
            }
            for json_key, lname in key_map.items():
                v = bp.get(json_key)
                if v is not None:
                    provided[lname] = float(v)

        biomass = self._get("material_origin.biomass_type", "hardwood")
        process = self._get("extraction_process.process_type", "kraft")

        bio_table   = LINKAGE_BY_BIOMASS_PROCESS.get(biomass, LINKAGE_BY_BIOMASS_PROCESS["hardwood"])
        lit_linkage = self._lookup_process_key(bio_table, process)

        raw: dict[str, float] = {}
        for lname in linkage_names:
            if lname in provided:
                raw[lname] = provided[lname]
            else:
                if lname in lit_linkage:
                    m, s = lit_linkage[lname]
                    raw[lname] = max(0.0, float(self._rng.normal(m, s)))
                else:
                    raw[lname] = 0.0

        # beta-1 rearrangement is currently unstable in the atomistic builder
        # and can produce disconnected graphs. Keep it disabled by default.
        allow_beta_1 = bool(self._get("simulation_config.allow_beta_1", False))
        if not allow_beta_1 and "beta-1" in raw:
            raw["beta-1"] = 0.0

        # Normalise to sum = 1
        total = sum(raw.values()) or 1.0
        self.priors["linkage_fractions"] = {k: v / total for k, v in raw.items()}

    # ── 3. Molecular weight ───────────────────────────────────────────────────

    def infer_MW(self) -> None:
        """Infer Mn, Mw, PDI from input or literature priors."""
        # Direct experimental values
        mw_input = self._get("batch_properties.Mw")
        if mw_input:
            Mn  = mw_input.get("Mn")
            Mw  = mw_input.get("Mw")
            PDI = mw_input.get("PDI")
            if Mn and Mw:
                self.priors["Mn"]  = float(Mn)
                self.priors["Mw"]  = float(Mw)
                self.priors["PDI"] = float(PDI) if PDI else float(Mw / Mn)
                return

        process   = self._get("extraction_process.process_type", "kraft")
        Mn_m, Mn_s, Mw_m, Mw_s = self._lookup_process_key(MW_BY_PROCESS, process)

        Mn  = max(200.0, float(self._rng.normal(Mn_m, Mn_s)))
        Mw  = max(Mn,    float(self._rng.normal(Mw_m, Mw_s)))
        PDI = Mw / Mn

        self.priors["Mn"]  = Mn
        self.priors["Mw"]  = Mw
        self.priors["PDI"] = PDI

    # ── 4. Branching index ─────────────────────────────────────────────────────

    def infer_branching(self) -> None:
        """Infer branching index."""
        # Direct experimental value
        bi = self._get("molecular_statistics.branching_index")
        if bi is not None:
            self.priors["branching_index"] = float(np.clip(bi, 0.0, 1.0))
            return

        process = self._get("extraction_process.process_type", "kraft")
        m, s    = self._lookup_process_key(BRANCHING_BY_PROCESS, process)
        self.priors["branching_index"] = self._sample(m, s, 0.0, 0.95)

    # ── 5. Condensation degree ─────────────────────────────────────────────────

    def infer_condensation(self) -> None:
        """Infer the degree of condensation (C-C / all inter-monomer bonds)."""
        cd = self._get("molecular_statistics.condensation_degree")
        if cd is not None:
            self.priors["condensation_degree"] = float(np.clip(cd, 0.0, 1.0))
            return

        process = self._get("extraction_process.process_type", "kraft")
        m, s    = self._lookup_process_key(CONDENSATION_BY_PROCESS, process)
        self.priors["condensation_degree"] = self._sample(m, s, 0.0, 1.0)

    # ── 6. Degree of polymerisation ───────────────────────────────────────────

    def infer_degree_of_polymerisation(self) -> None:
        """
        Compute mean_DP and max_DP from Mn / Mw and monomer composition.
        Must be called after infer_SGH() and infer_MW().
        """
        # Composition-weighted average monomer MW
        S = self.priors.get("S_fraction", 0.5)
        G = self.priors.get("G_fraction", 0.4)
        H = self.priors.get("H_fraction", 0.1)
        avg_monomer_MW = (
            H * MONOMER_MW["H"] + G * MONOMER_MW["G"] + S * MONOMER_MW["S"]
        )
        self.priors["avg_monomer_MW"] = avg_monomer_MW

        # Check molecular_statistics first
        chain_stats = self._get("molecular_statistics.chain_length_distribution")
        if chain_stats:
            mean_dp = chain_stats.get("mean_DP")
            var_dp  = chain_stats.get("variance_DP")
            if mean_dp:
                self.priors["mean_DP"] = float(mean_dp)
                self.priors["max_DP"]  = int(mean_dp * 2.5)
                return

        Mn = self.priors.get("Mn", 2000.0)
        Mw = self.priors.get("Mw", 5000.0)
        self.priors["mean_DP"] = max(2, int(round(Mn / avg_monomer_MW)))
        self.priors["max_DP"]  = max(4, int(round(Mw * 2.0 / avg_monomer_MW)))

    # ── Main entry point ───────────────────────────────────────────────────────

    def run(self) -> dict:
        """
        Execute the full prior estimation pipeline.

        Returns
        -------
        priors : dict
            Dictionary of all estimated structural priors.
        """
        self.infer_SGH()
        self.infer_linkage_distribution()
        self.infer_MW()
        self.infer_branching()
        self.infer_condensation()
        self.infer_degree_of_polymerisation()
        return self.priors

    def summary(self) -> str:
        """Return a human-readable summary string of the estimated priors."""
        lines = ["=== LignoForge Prior Estimates ==="]
        lines.append(
            f"Monomer fractions : H={self.priors.get('H_fraction',0):.3f}  "
            f"G={self.priors.get('G_fraction',0):.3f}  "
            f"S={self.priors.get('S_fraction',0):.3f}"
        )
        lines.append(
            f"Molecular weight  : Mₙ={self.priors.get('Mn',0):.0f}  "
            f"Mw={self.priors.get('Mw',0):.0f}  "
            f"PDI={self.priors.get('PDI',0):.2f}"
        )
        lines.append(
            f"Mean DP / max DP  : {self.priors.get('mean_DP','-')} / "
            f"{self.priors.get('max_DP','-')}"
        )
        lines.append(f"Branching index   : {self.priors.get('branching_index',0):.3f}")
        lines.append(f"Condensation deg. : {self.priors.get('condensation_degree',0):.3f}")
        if "linkage_fractions" in self.priors:
            lines.append("Linkage fractions :")
            for k, v in self.priors["linkage_fractions"].items():
                lines.append(f"    {k:<14s}: {v:.3f}")
        return "\n".join(lines)
