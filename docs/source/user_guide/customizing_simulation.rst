==========================
Customising the Simulation
==========================

This guide explains every simulation hyperparameter and how to tune it
for specific use cases such as stiff or highly polydisperse lignins.

Performance and parallelism
-----------------------------

3-D coordinate generation (RDKit ETKDGv3 + MMFF geometry optimisation)
is the most time-consuming step.  LignoForge runs it in parallel across
all available CPU cores by default using :class:`~concurrent.futures.ProcessPoolExecutor`.

Control parallelism via the ``n_workers`` export option::

    results = pipeline.run(
        simulation_overrides={
            "n_population": 20,
            "n_workers":    4,    # cap at 4 processes
        }
    )

Or from the command line::

    python examples/demo_run.py --n-chains 20 --workers 4

Set ``n_workers=1`` to disable parallelism entirely (useful for debugging
or when profiling memory on a single chain).

To skip 3-D generation completely (topology and SMILES still produced)::

    results = pipeline.run(
        simulation_overrides={"generate_3d": False}
    )

Geometry optimisation tuning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The MMFF94 optimiser runs in chunks and reports progress every 25 % of
the iteration budget, so you always know how far along each chain is.
The default budget (``max_uff_iterations=300``) is a pragmatic balance
between structural quality and speed.  For publication-quality geometries
increase it; for rapid screening decrease it::

    simulation_overrides={
        "max_uff_iterations": 150,   # faster, coarser geometry
    }

The optimiser exits early if it converges before reaching the iteration
limit, so the actual number of iterations is often much lower than the
maximum.

Simulation parameter reference
--------------------------------

All parameters can be passed either via ``simulation_config`` in the JSON
input or via ``simulation_overrides`` in
:meth:`~lignoforge.pipeline.pipeline.LigninPipeline.run`.

Chain size
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 12 55

   * - Parameter
     - Type
     - Description
   * - ``expected_size``
     - float
     - Mean target chain size in monomer count (or g/mol when
       ``size_in_MW=True``).  Inferred from ``mean_DP`` unless overridden.
   * - ``max_size``
     - float
     - Hard upper cap.  Set to 3–4 × ``expected_size`` for realistic
       polydispersity.
   * - ``distribution_scaling``
     - float
     - Controls the standard deviation of the size distribution as
       :math:`\sigma = k \times \text{expected\_size}`.  Larger *k* →
       wider distribution → higher PDI.  Derived from PDI by default.
   * - ``size_in_MW``
     - bool
     - If ``True``, ``expected_size`` and ``max_size`` are in g/mol
       instead of monomer count.

Metropolis temperatures
~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 12 55

   * - Parameter
     - Typical range
     - Description
   * - ``Tmetro``
     - 1 – 300 K
     - Inner-loop temperature.  Lower → stricter per-chain acceptance
       → each chain converges faster but may be trapped in a local
       minimum.  Higher → more exploratory.
   * - ``Tmetro_out``
     - 1 – 1000 K
     - Outer-loop temperature.  Controls how much the population
       distance is allowed to worsen before rejecting a new chain.

**Rule of thumb**: start with ``Tmetro = Tmetro_out = 10`` for a good
balance of exploration and convergence.  Decrease both if the final
population statistics are far from the target; increase both if the
algorithm seems to get stuck.

Iteration limits
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 12 55

   * - Parameter
     - Default
     - Description
   * - ``i_max``
     - 1000
     - Maximum inner-loop MC steps per chain.  Increase for larger
       target chain sizes.
   * - ``i_max_out``
     - 1000
     - Maximum outer-loop iterations.  Increase if the population is
       not yet full after the default budget.
   * - ``i_max_ring``
     - 500
     - Maximum ring-closure attempts per population.  Set to 0 to
       disable ring formation entirely.

Population size
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 12 55

   * - Parameter
     - Default
     - Description
   * - ``n_population``
     - 50
     - Target number of accepted chains.  Larger populations give better
       statistics but take proportionally longer.

Branching
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 12 55

   * - Parameter
     - Range
     - Description
   * - ``branching_propensity``
     - [0, 1]
     - Probability of attaching the next monomer to an interior node
       (branching) rather than a terminal node.  Set to 0 for linear
       chains.  Inferred from ``branching_index`` by default.

Metric weights
~~~~~~~~~~~~~~

By default all metrics contribute equally to the distance function.  You
can upweight specific metrics by passing a ``metrics_weights`` array whose
length equals the number of active metrics (3 monomer + 7 linkage +
optional extras)::

    import numpy as np

    # Double the weight of beta-O-4 linkage fraction
    # metrics order: [H, G, S,  4-O-5, alpha-O-4, beta-O-4, 5-5, beta-5, beta-beta, beta-1]
    weights = np.ones(10)
    weights[5] = 2.0        # beta-O-4 is index 5

    results = pipeline.run(
        simulation_overrides={"metrics_weights": weights.tolist()}
    )

Practical recipes
------------------

High-PDI lignin (Mw/Mn > 3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    simulation_overrides = {
        "distribution_scaling": 2.5,
        "expected_size":        20,
        "max_size":             120,
    }

Softwood (G-rich, no branching)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    simulation_overrides = {
        "branching_propensity": 0.0,
        "i_max_ring":           0,
    }

Fast exploratory run
~~~~~~~~~~~~~~~~~~~~~

::

    simulation_overrides = {
        "n_population": 10,
        "i_max":        200,
        "i_max_out":    200,
        "i_max_ring":   50,
        "Tmetro":       50.0,
        "Tmetro_out":   100.0,
    }

Large high-quality library
~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    simulation_overrides = {
        "n_population": 200,
        "i_max":        5000,
        "i_max_out":    5000,
        "i_max_ring":   2000,
        "Tmetro":       5.0,
        "Tmetro_out":   10.0,
    }

Accessing the Simulation object directly
-----------------------------------------

For maximum control you can instantiate
:class:`~lignoforge.simulation.population.Simulation` directly::

    from lignoforge.simulation.population import Simulation

    sim = Simulation(
        linkage_distribution_input=[0.0, 0.05, 0.52, 0.12, 0.14, 0.10, 0.07],
        monomer_distribution_input=[0.01, 0.97, 0.02],
        expected_size=20,
        max_size=80,
        distribution_scaling=1.2,
        Tmetro=8.0,
        Tmetro_out=15.0,
        n_population=30,
        i_max=2000,
        i_max_out=2000,
        i_max_ring=500,
        seed_init=42,
        library_name="pine_library",
        results_name="output",
        branching_propensity=0.02,
        verbose=True,
    )
    sim.run()
    polymers = sim.P_population
