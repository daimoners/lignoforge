"""
Parameter translator from high-level priors to simulation-ready arguments.
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

from lignoforge.core.rules import linkage_names


class ParameterTranslator:
    """
    Convert prior estimates into the exact parameter format expected by
    `lignoforge.simulation.Simulation`.
    """

    def __init__(self, priors: dict, input_data: dict | None = None):
        self.priors = priors
        self.input_data = input_data or {}

    def _linkage_vector(self) -> List[float]:
        fractions = self.priors.get("linkage_fractions", {})
        v = np.array([float(fractions.get(name, 0.0)) for name in linkage_names], dtype=float)
        if np.sum(v) <= 0:
            v = np.ones(len(linkage_names), dtype=float)
        v /= np.sum(v)
        return v.tolist()

    def _monomer_vector(self) -> List[float]:
        H = float(self.priors.get("H_fraction", 0.2))
        G = float(self.priors.get("G_fraction", 0.5))
        S = float(self.priors.get("S_fraction", 0.3))
        v = np.array([H, G, S], dtype=float)
        if np.sum(v) <= 0:
            v = np.array([0.2, 0.5, 0.3], dtype=float)
        v /= np.sum(v)
        return v.tolist()

    def _sizes(self) -> tuple[int, int]:
        mean_dp = int(max(2, round(float(self.priors.get("mean_DP", 20)))))
        max_dp  = int(max(mean_dp + 1, round(float(self.priors.get("max_DP", mean_dp * 3)))))
        return mean_dp, max_dp

    def _distribution_scaling(self) -> float:
        pdi = float(self.priors.get("PDI", 2.0))
        return max(0.1, min(5.0, pdi - 1.0))

    def _temperatures(self) -> tuple[float, float]:
        sim_cfg = self.input_data.get("simulation_config", {})
        t_in  = sim_cfg.get("Tmetro")
        t_out = sim_cfg.get("Tmetro_out")
        Tmetro     = float(t_in)  if t_in  is not None else 298.15
        Tmetro_out = float(t_out) if t_out is not None else 500.0
        return Tmetro, Tmetro_out

    def _population_size(self) -> int:
        sim_cfg = self.input_data.get("simulation_config", {})
        n = sim_cfg.get("n_population")
        if n is not None:
            return max(1, int(n))
        return 50

    def to_simulation_kwargs(self) -> Dict:
        """
        Return kwargs dictionary ready to be expanded into `Simulation(...)`.
        """
        expected_size, max_size = self._sizes()
        Tmetro, Tmetro_out = self._temperatures()

        kwargs = {
            "linkage_distribution_input": self._linkage_vector(),
            "monomer_distribution_input": self._monomer_vector(),
            "expected_size": expected_size,
            "max_size": max_size,
            "distribution_scaling": self._distribution_scaling(),
            "Tmetro": Tmetro,
            "Tmetro_out": Tmetro_out,
            "branching_propensity": float(self.priors.get("branching_index", 0.1)),
            "n_population": self._population_size(),
            "population_metrics": [
                float(self.priors.get("Mn", 3000.0)),
                float(self.priors.get("Mw", 6000.0)),
            ],
        }
        return kwargs
