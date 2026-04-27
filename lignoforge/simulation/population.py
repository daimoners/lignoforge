"""
Population-level MCMC simulation for generating an ensemble of lignin
polymer structures that collectively match a target metrics vector.
"""

from __future__ import annotations

import os
import shutil
import time
from typing import Optional, Tuple

import numpy as np

from lignoforge.core.rules import monomer_types, linkage_names
from lignoforge.core.characterization import get_counts_polymer, Population
from lignoforge.core.utils import (
    nparray,
    counts_to_metrics,
    MW_array_to_number_average,
    MW_array_to_weight_average,
    cal_distance,
    metrics_array_to_dict,
    plot_metrics,
    plot_distance_trajectory,
)
from lignoforge.simulation.trajectory import Trajectory

# Boltzmann constant (eV K⁻¹)
_KB = 8.617333262145e-5


class Simulation(Trajectory):
    """
    Full population MCMC: outer loop accepts/rejects entire polymer
    additions to an ensemble, inner loop (Trajectory) generates each chain.

    Parameters
    ----------
    linkage_distribution_input : list  – 7-element linkage probability vector
    monomer_distribution_input : list  – 3-element monomer probability vector [H, G, S]
    expected_size              : float – target mean polymer size
    max_size                   : float – maximum allowed polymer size
    distribution_scaling       : float – controls size distribution width
    Tmetro                     : float – inner-loop Metropolis temperature
    Tmetro_out                 : float – outer-loop Metropolis temperature
    seed_init                  : int   – global random seed
    library_name               : str   – name of this lignin library
    results_name               : str   – parent folder for results
    n_population               : int   – target ensemble size
    i_max                      : int   – max inner-loop MC steps per chain
    i_max_out                  : int   – max outer-loop iterations
    i_max_ring                 : int   – max ring-closure MC steps
    additional_metrics         : list  – extra target metrics (e.g. branching coeff)
    population_metrics         : list  – population-level targets [Mn, Mw] if desired
    size_in_MW                 : bool  – use MW instead of monomer count as size
    branching_propensity       : float – probability of branching per step
    metrics_weights            : array – per-metric weights for distance calculation
    verbose                    : bool
    show_plots                 : bool
    save_path                  : str   – root directory for output
    """

    def __init__(
        self,
        linkage_distribution_input:  list,
        monomer_distribution_input:  list,
        expected_size:               float,
        max_size:                    float,
        distribution_scaling:        float,
        Tmetro:                      float,
        Tmetro_out:                  float,
        seed_init:                   int   = 1,
        library_name:                str   = "lignin_x",
        results_name:                str   = "results",
        trial_index:                 Optional[int] = None,
        n_population:                int   = 100,
        i_max:                       int   = 1000,
        i_max_out:                   int   = 1000,
        i_max_ring:                  int   = 500,
        additional_metrics:          Optional[list] = None,
        population_metrics:          Optional[list] = None,
        size_in_MW:                  bool  = False,
        branching_propensity:        Optional[float] = None,
        metrics_weights:             Optional[nparray] = None,
        verbose:                     bool  = True,
        show_plots:                  bool  = True,
        save_path:                   str   = None,
    ):
        super().__init__(
            linkage_distribution_input,
            monomer_distribution_input,
            Tmetro,
            expected_size,
            max_size,
            distribution_scaling,
            additional_metrics,
            size_in_MW,
            branching_propensity,
            metrics_weights,
            verbose,
        )

        self.save_path = save_path or os.getcwd()

        # ── I/O paths ─────────────────────────────────────────────────────────
        results_parent = os.path.join(self.save_path, results_name, library_name)
        os.makedirs(results_parent, exist_ok=True)

        if trial_index is None:
            trial_index = sum(
                1 for e in os.scandir(results_parent) if e.is_dir()
            )
        results_path = os.path.join(results_parent, f"i{trial_index}")
        if os.path.exists(results_path):
            shutil.rmtree(results_path)
        os.makedirs(results_path, exist_ok=True)

        print(f"Starting trial No.{trial_index}  →  {results_path}")

        # ── Parameters ────────────────────────────────────────────────────────
        self.library_name        = library_name
        self.results_name        = results_name
        self.results_path        = results_path
        self.trial_index         = trial_index
        self.n_population        = n_population
        self.Tmetro_out          = Tmetro_out
        self.seed_init           = seed_init
        self.i_max               = i_max
        self.i_max_out           = i_max_out
        self.i_max_ring          = i_max_ring
        self.population_metrics  = population_metrics
        self.show_plots          = show_plots

        self.max_MW = self.max_size * 150 if not self.size_in_MW else self.max_size

        # Results attributes (populated after run())
        self.P_population:           Optional[list]    = None
        self.metrics_current_dict:   Optional[dict]    = None
        self.metrics_target_dict:    Optional[dict]    = None
        self.distance:               Optional[list]    = None
        self.distance_final:         Optional[float]   = None

    # ── Population generator ──────────────────────────────────────────────────

    def run(self) -> list:
        """
        Execute the full population MCMC and return the accepted polymer list.
        """
        np.random.seed(self.seed_init)

        log_path = os.path.join(self.results_path, self.library_name + ".out")
        with open(log_path, "w") as log:
            log.write(f"=== {self.library_name} | trial {self.trial_index} ===\n")

            # Per-chain weights (drop population-level metrics if present)
            w_individual = self.metrics_weights
            if self.metrics_weights is not None and self.population_metrics:
                w_individual = self.metrics_weights[:-2]

            traj = Trajectory(
                self.linkage_distribution,
                self.monomer_distribution,
                self.Tmetro,
                self.expected_size,
                self.max_size,
                self.distribution_scaling,
                self.additional_metrics,
                self.size_in_MW,
                self.branching_propensity,
                w_individual,
                verbose=self.verbose,
                log_file=log,
            )

            # Build target metrics vector (extended with population MW targets)
            metrics_names = monomer_types + linkage_names
            if self.additional:
                metrics_names += ["branching_coeff"]
            if self.population_metrics:
                metrics_names += ["Mn", "Mw"]

            metrics_target = traj.metrics_target.copy()
            metrics_target_orig = traj.metrics_target.copy()
            if self.population_metrics:
                pop_norm = [p / self.max_MW for p in self.population_metrics]
                metrics_target      = np.append(metrics_target, pop_norm)
                metrics_target_orig = np.append(
                    metrics_target_orig, self.population_metrics
                )
            self.metrics_target_dict = metrics_array_to_dict(
                metrics_target_orig, metrics_names
            )

            # ── Outer MC loop ──────────────────────────────────────────────────
            P_pop, counts_pop, MW_pop = [], [], []
            n_polymers = 0
            d_average  = np.inf
            rseed      = 0
            i_step     = 0
            start      = time.time()

            while n_polymers < self.n_population and i_step <= self.i_max_out:
                P_i, dist_i, mc_i, _ = traj.run_MCMC(rseed, self.i_max)
                counts_P, mc_P, MW_P = get_counts_polymer(
                    P_i, additional=self.additional, cal_MW=True
                )
                if MW_P is None or MW_P < 100:
                    rseed += 1; i_step += 1; continue

                metrics_P = counts_to_metrics(counts_P, additional=self.additional)

                # Tentatively add polymer
                counts_pop_try = counts_pop + [counts_P]
                MW_pop_try     = MW_pop + [MW_P]
                counts_sum     = np.sum(np.array(counts_pop_try), axis=0)
                m_avg          = counts_to_metrics(counts_sum, additional=self.additional)

                if self.population_metrics:
                    Mn_avg = MW_array_to_number_average(np.array(MW_pop_try))
                    Mw_avg = MW_array_to_weight_average(np.array(MW_pop_try))
                    m_avg  = np.append(m_avg, [Mn_avg/self.max_MW, Mw_avg/self.max_MW])

                d_new   = cal_distance(metrics_target, m_avg, self.metrics_weights)
                delta_d = d_new - d_average
                accept  = False

                if delta_d <= 0 or self.Tmetro_out == np.inf:
                    accept = True
                elif self.Tmetro_out > 0:
                    w = np.exp(-delta_d / _KB / self.Tmetro_out)
                    if np.random.rand() <= w:
                        accept = True

                if accept:
                    d_average = d_new
                    P_pop.append(P_i)
                    counts_pop = counts_pop_try
                    MW_pop     = MW_pop_try
                    n_polymers += 1
                    log.write(f"\tn_polymer {n_polymers-1} added on iter {i_step}\n")
                    if self.verbose:
                        print(f"\t  [{n_polymers}/{self.n_population}] polymer accepted")
                else:
                    log.write("\tPolymer addition rejected\n")

                rseed += 1; i_step += 1

            elapsed = time.time() - start
            log.write(f"All chains done: {elapsed/60:.2f} min\n")
            print(f"Chain generation: {elapsed/60:.2f} min")

            # ── Ring closure pass ──────────────────────────────────────────────
            P_pop = self._run_ring_pass(P_pop, [rseed + i for i in range(len(P_pop))],
                                        traj, metrics_target, counts_pop, MW_pop,
                                        d_average, log)

        # ── Characterise population ────────────────────────────────────────────
        pop_obj = Population(
            P_pop, self.library_name,
            ResultsName=self.results_name,
            TrialIndex=str(self.trial_index),
            save_path=self.save_path,
        )
        pop_obj.analyze()
        pop_obj.population_stats()
        self.metrics_current_dict = pop_obj.get_metrics_mean(additional=self.additional)

        self.P_population  = P_pop
        self.distance_final = d_average
        return P_pop

    # ── Ring-closure helper ───────────────────────────────────────────────────

    def _run_ring_pass(
        self, P_pop, rseeds, traj, metrics_target,
        counts_pop, MW_pop, d_current, log,
    ) -> list:
        """Attempt ring closures on all accepted polymers."""
        if self.i_max_ring <= 0:
            return P_pop
        if self.branching_propensity is not None and self.branching_propensity == 0.0:
            return P_pop

        P_ring, counts_ring, MW_ring = [], [], []
        start = time.time()

        for i, (Pi, ri) in enumerate(zip(P_pop, rseeds)):
            P_i, _, acc, _ = traj.run_MCMC_ring(Pi, ri, self.i_max_ring)
            counts_P, mc_P, MW_P = get_counts_polymer(
                P_i, additional=self.additional, cal_MW=True
            )
            if MW_P is None or MW_P < 100:
                P_ring.append(Pi); counts_ring.append(counts_pop[i]); MW_ring.append(MW_pop[i])
                continue
            log.write(f"\t{acc} ring(s) added to polymer {i}\n")
            if acc > 0 and self.verbose:
                print(f"\t  {acc} ring(s) closed on polymer {i}")
            P_ring.append(P_i)
            counts_ring.append(counts_P)
            MW_ring.append(MW_P)

        elapsed = time.time() - start
        print(f"Ring addition: {elapsed/60:.2f} min")

        # Accept the ring-closure pass only if it improves the distance
        counts_sum = np.sum(np.array(counts_ring), axis=0)
        m_ring     = counts_to_metrics(counts_sum, additional=self.additional)
        if self.population_metrics:
            Mn_avg = MW_array_to_number_average(np.array(MW_ring))
            Mw_avg = MW_array_to_weight_average(np.array(MW_ring))
            m_ring = np.append(m_ring, [Mn_avg/self.max_MW, Mw_avg/self.max_MW])
        d_ring = cal_distance(metrics_target, m_ring, self.metrics_weights)

        if d_ring <= d_current:
            return P_ring
        return P_pop
