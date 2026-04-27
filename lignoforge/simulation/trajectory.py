"""
kMC / MCMC trajectory engine for growing a single lignin polymer chain.
"""

from __future__ import annotations

import time
from typing import Optional, Tuple

import numpy as np

from lignoforge.core.rules import monomer_types, linkage_names
from lignoforge.core.monomer import Monomer
from lignoforge.core.polymer import Polymer
from lignoforge.core.characterization import get_metrics_polymer, get_counts_polymer
from lignoforge.core.utils import (
    nparray,
    set_random_state,
    generate_random_monomer,
    generate_random_linkage,
    generate_random_branching_state,
    generate_random_size_from_distribution,
    cal_distance,
    counts_to_metrics,
)

# Boltzmann constant (eV K⁻¹)
_KB = 8.617333262145e-5


class Trajectory:
    """
    Single kMC/MCMC trajectory: grows one polymer chain from one initial
    monomer until the chain reaches the target size.

    The acceptance criterion uses a Metropolis–Hastings rule where the
    "energy" is the Euclidean distance between the current polymer metrics
    and the target metrics.
    """

    def __init__(
        self,
        linkage_distribution_input:  list,
        monomer_distribution_input:  list,
        Tmetro:                      float,
        expected_size:               float,
        max_size:                    float,
        distribution_scaling:        float = 1.0,
        additional_metrics:          Optional[list | float] = None,
        size_in_MW:                  bool  = False,
        branching_propensity:        Optional[float] = None,
        metrics_weights:             Optional[nparray] = None,
        verbose:                     bool  = True,
        log_file:                    Optional[object] = None,
    ):
        self.linkage_distribution = (
            np.array(linkage_distribution_input, dtype=float)
            / np.sum(linkage_distribution_input)
        )
        self.monomer_distribution = (
            np.array(monomer_distribution_input, dtype=float)
            / np.sum(monomer_distribution_input)
        )
        self.expected_size       = expected_size
        self.max_size            = max_size
        self.distribution_scaling = distribution_scaling
        self.size_in_MW          = size_in_MW
        self.Tmetro              = Tmetro
        self.branching_propensity = branching_propensity
        self.metrics_weights     = metrics_weights
        self.verbose             = verbose
        self.log_file            = log_file
        self.max_monomer_count   = 200  # hard safety cap

        # Build target metrics vector
        metrics_target = np.concatenate(
            [self.monomer_distribution, self.linkage_distribution]
        )
        self.additional = additional_metrics is not None
        if self.additional:
            if isinstance(additional_metrics, (int, float)):
                additional_metrics = [additional_metrics]
            additional_metrics = np.array(additional_metrics, dtype=float)
            metrics_target = np.concatenate([metrics_target, additional_metrics])

        self.metrics_target    = metrics_target
        self.additional_metrics = additional_metrics

    # ── Chain growth ──────────────────────────────────────────────────────────

    def run_MCMC(
        self,
        rseed:  int,
        i_max:  int = 500,
    ) -> Tuple[Polymer, list, int, int]:
        """
        Grow a single polymer chain via Metropolis MC.

        Returns
        -------
        polymer_final  : Polymer
        distance       : list of distance values at each accepted step
        monomer_count  : int    – final monomer count
        i_step         : int    – total number of MC steps attempted
        """
        rstate = set_random_state(rseed)

        # Draw a target size from the size distribution
        stop_size = generate_random_size_from_distribution(
            self.expected_size, self.max_size,
            self.distribution_scaling, self.size_in_MW,
        )

        branching_state = None

        # Initialise with a random monomer
        mtype   = generate_random_monomer(self.monomer_distribution, rstate)
        M_init  = Monomer(mtype)
        polymer = Polymer(M_init, verbose=False)

        metrics_P, monomer_count, MW_P = get_metrics_polymer(
            polymer, additional=self.additional, cal_MW=True
        )
        d = cal_distance(self.metrics_target, metrics_P, self.metrics_weights)

        # Form a dimer before entering the MC loop
        ltype = generate_random_linkage(self.linkage_distribution, rstate)
        mtype = generate_random_monomer(self.monomer_distribution, rstate)
        polymer.add_specific_monomer(mtype, ltype)
        metrics_P, monomer_count, MW_P = get_metrics_polymer(
            polymer, additional=self.additional, cal_MW=True
        )
        d = cal_distance(self.metrics_target, metrics_P, self.metrics_weights)
        current_size = MW_P if self.size_in_MW else monomer_count

        # ── Main MC loop ───────────────────────────────────────────────────────
        acceptance, distance, i_step = 0, [d], 0
        start = time.time()

        while current_size <= stop_size and i_step <= i_max:
            # Branching state for this step
            if self.branching_propensity is not None:
                if self.branching_propensity > 0.0:
                    branching_state = generate_random_branching_state(
                        self.branching_propensity
                    )
                else:
                    branching_state = False

            polymer_i = Polymer(polymer, verbose=False)
            ltype = generate_random_linkage(self.linkage_distribution, rstate)
            mtype = generate_random_monomer(self.monomer_distribution, rstate)
            added = polymer_i.add_specific_monomer(mtype, ltype, branching_state)

            if added:
                metrics_P, monomer_count, MW_P = get_metrics_polymer(
                    polymer_i, additional=self.additional, cal_MW=True
                )
                if MW_P is not None and MW_P < 100:
                    i_step += 1
                    continue

                d_new     = cal_distance(self.metrics_target, metrics_P, self.metrics_weights)
                delta_d   = d_new - d
                accepted  = False

                if delta_d <= 0 or self.Tmetro == np.inf:
                    accepted = True
                elif self.Tmetro > 0:
                    w = np.exp(-delta_d / _KB / self.Tmetro)
                    if np.random.rand() <= w:
                        accepted = True

                if accepted:
                    d        = d_new
                    polymer  = polymer_i
                    acceptance += 1
                    current_size = MW_P if self.size_in_MW else monomer_count

            i_step += 1
            distance.append(d)

        elapsed = time.time() - start
        if self.log_file:
            self.log_file.write(
                f"\t\tPolymerization: {elapsed/60:.2f} min, "
                f"{acceptance} accepted / {i_step} steps\n"
            )
        return polymer, distance, int(monomer_count), i_step

    # ── Ring closure ──────────────────────────────────────────────────────────

    def run_MCMC_ring(
        self,
        polymer_init: Polymer,
        rseed:        int,
        i_max_ring:   int = 100,
    ) -> Tuple[Polymer, list, int, int]:
        """
        Attempt ring closures on an existing polymer via Metropolis MC.

        Returns
        -------
        polymer_final  : Polymer
        distance       : list
        acceptance_count: int
        i_step         : int
        """
        rstate  = set_random_state(rseed)
        polymer = Polymer(polymer_init, verbose=False)

        metrics_P, monomer_count, MW_P = get_metrics_polymer(
            polymer, additional=self.additional, cal_MW=True
        )
        d = cal_distance(self.metrics_target, metrics_P, self.metrics_weights)

        acceptance, distance, i_step = 0, [d], 0
        start = time.time()

        while i_step <= i_max_ring:
            polymer_i = Polymer(polymer, verbose=False)
            ltype     = generate_random_linkage(self.linkage_distribution, rstate)
            added     = polymer_i.add_specific_ring(ltype)

            if added:
                metrics_P, monomer_count, MW_P = get_metrics_polymer(
                    polymer_i, additional=self.additional, cal_MW=True
                )
                if MW_P is not None and MW_P < 100:
                    i_step += 1
                    continue

                d_new   = cal_distance(self.metrics_target, metrics_P, self.metrics_weights)
                delta_d = d_new - d
                accepted = False

                if delta_d <= 0 or self.Tmetro == np.inf:
                    accepted = True
                elif self.Tmetro > 0:
                    w = np.exp(-delta_d / _KB / self.Tmetro)
                    if np.random.rand() <= w:
                        accepted = True

                if accepted:
                    d       = d_new
                    polymer = polymer_i
                    acceptance += 1

            distance.append(d)
            i_step += 1

        elapsed = time.time() - start
        if self.log_file:
            self.log_file.write(
                f"\t\tRing addition: {elapsed/60:.2f} min, "
                f"{acceptance} accepted / {i_step} steps\n"
            )
        return polymer, distance, acceptance, i_step
