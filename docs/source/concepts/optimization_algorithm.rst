=======================
Optimisation Algorithm
=======================

LignoForge employs a **hierarchical Metropolis Monte Carlo** (hMMC) scheme
to generate a population of polymer chains whose ensemble-averaged
structural metrics match experimental targets.  The algorithm proceeds in
three nested loops.

Overview
--------

::

    ┌──────────────────────────────────────────────────────────┐
    │  Outer loop  (population MCMC)                           │
    │  ─────────────────                                       │
    │  Accept / reject insertion of a new chain into the       │
    │  population by comparing the population-level distance.  │
    │                                                          │
    │  ┌─────────────────────────────────────────────────────┐ │
    │  │  Inner loop  (single-chain kMC / MCMC)              │ │
    │  │  ─────────────────────────────────────────          │ │
    │  │  Grow one polymer chain monomer by monomer.         │ │
    │  │  Accept each monomer addition by comparing the      │ │
    │  │  per-chain distance.                                │ │
    │  └─────────────────────────────────────────────────────┘ │
    │                                                          │
    │  Ring-closure loop  (intramolecular)                     │
    │  ─────────────────────────────────                       │
    │  Attempt ring motifs (β-5 phenylcoumaran, β-β resinol)   │
    │  on accepted chains.                                     │
    └──────────────────────────────────────────────────────────┘

Distance metric
---------------

The Euclidean distance between a simulated state and the experimental
target is:

.. math::

   d = \sqrt{\sum_{k} \, w_k \left( m_k^{\text{sim}} - m_k^{\text{exp}} \right)^2}

where:

* :math:`m_k` ranges over all active metrics: the 3 monomer fractions
  (H, G, S), the 7 linkage fractions, and optional additional metrics
  (branching coefficient, M\ :sub:`n`, M\ :sub:`w`).
* :math:`w_k` are optional per-metric weights (default: 1).

For a population the same formula is applied to the ensemble-averaged
metrics rather than per-chain values.

Metropolis acceptance criterion
----------------------------------

Each proposed change (monomer addition or chain insertion) is accepted
with probability:

.. math::

   P_{\text{accept}} = \min\!\left(1,\; e^{-\,\frac{\Delta d}{k_B T_{\text{metro}}}}\right)

where:

* :math:`\Delta d = d_{\text{new}} - d_{\text{old}}` (negative = improvement).
* :math:`k_B = 8.617 \times 10^{-5}\ \text{eV K}^{-1}` (Boltzmann constant).
* :math:`T_{\text{metro}}` is the Metropolis temperature — a hyperparameter
  that controls the acceptance rate of worsening moves.

Setting :math:`T_{\text{metro}} \to \infty` → all moves accepted (random
walk).  Lower temperatures → stricter acceptance, faster convergence but
risk of trapping.

Loop 1: Inner loop (single-chain growth)
-----------------------------------------

Implemented in :class:`~lignoforge.simulation.trajectory.Trajectory`.

1. Draw a target chain size from a truncated-normal distribution centred on
   ``expected_size`` with width controlled by ``distribution_scaling``.
2. Initialise with a random monomer sampled from
   ``monomer_distribution_input``.
3. Repeatedly propose monomer additions:

   a. Sample a linkage type from ``linkage_distribution_input``.
   b. Sample a monomer type from ``monomer_distribution_input``.
   c. Add the monomer via the rejection-free kMC step
      (:meth:`~lignoforge.core.polymer.Polymer.add_specific_monomer`).
   d. Compute the new per-chain distance.
   e. Accept / reject with the Metropolis criterion (``Tmetro``).

4. Stop when the chain reaches the target size or ``i_max`` steps are
   exhausted.

The inner loop is **rejection-free at the chemical level**: a chemically
valid (monomer, linkage) pair is always found using the O(1) lookup tables
in :mod:`lignoforge.core.rules`.  The Metropolis criterion controls
whether the accepted addition is kept in the growing chain, not whether
chemistry is possible.

Loop 2: Outer loop (population optimisation)
----------------------------------------------

Implemented in :class:`~lignoforge.simulation.population.Simulation`
(inherits from ``Trajectory``).

1. Start with an empty population.
2. For each outer iteration up to ``i_max_out``:

   a. Grow a candidate chain using the inner loop.
   b. Compute the new population-averaged metrics and population distance.
   c. Accept / reject the candidate with the Metropolis criterion
      (``Tmetro_out``).
   d. If accepted and population size < ``n_population``: insert the chain.
   e. If the population is full: optionally replace the worst-fitting chain.

3. Stop when ``n_population`` accepted chains are assembled or ``i_max_out``
   is exhausted.

Loop 3: Ring-closure loop
---------------------------

After the outer loop, an intramolecular ring-closure pass runs on every
accepted chain:

1. For each ring-closure attempt:

   a. Try to add a ring motif (β-5 or β-β) via
      :meth:`~lignoforge.core.polymer.Polymer.add_random_ring`.
   b. Accept / reject with the Metropolis criterion.

2. Stop after ``i_max_ring`` attempts.

Ring closures introduce the phenylcoumaran and resinol motifs that are
topologically impossible in a purely linear growth process.

Polymer size distribution
--------------------------

The target chain size for each inner-loop run is drawn from a
*truncated normal* distribution:

.. math::

   n_{\text{target}} \sim \mathcal{N}\!\left(\mu = n_{\text{expected}},\;
   \sigma = k \cdot n_{\text{expected}}\right) \;\;\text{truncated to } [2,\, n_{\text{max}}]

where :math:`k` = ``distribution_scaling``.  Larger :math:`k` → wider
size distribution → higher polydispersity index.

When ``size_in_MW = True``, the target is drawn in g/mol using the
composition-weighted average monomer molecular weight.

Convergence diagnostics
------------------------

The population MCMC writes a ``.out`` log file to ``results_path`` after
every accepted chain.  The log records per-chain acceptance rates and
elapsed times.

The exporter writes a ``population_statistics.json`` file containing the
mean and standard deviation of all structural metrics across the accepted
population.  Comparing these to the experimental targets gives a direct
measure of convergence.

Hyperparameters reference
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 28 14 55

   * - Parameter
     - Default
     - Effect
   * - ``Tmetro``
     - 298 K
     - Inner loop temperature.  Higher → more exploratory, slower
       convergence.
   * - ``Tmetro_out``
     - 500 K
     - Outer loop temperature.  Controls population diversity.
   * - ``expected_size``
     - from priors
     - Mean chain length (monomers or g/mol if ``size_in_MW``).
   * - ``max_size``
     - from priors
     - Hard upper cap on chain length.
   * - ``distribution_scaling``
     - from PDI
     - Width of the size distribution; higher → more polydisperse.
   * - ``n_population``
     - 50
     - Target ensemble size.
   * - ``i_max``
     - 1000
     - Max inner-loop steps per chain.
   * - ``i_max_out``
     - 1000
     - Max outer-loop iterations.
   * - ``i_max_ring``
     - 500
     - Max ring-closure attempts per population.
   * - ``branching_propensity``
     - from priors
     - Probability of attaching to an interior node (branching).
   * - ``metrics_weights``
     - None (all 1)
     - Per-metric weights vector for the distance function.
