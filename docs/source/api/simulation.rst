==========
Simulation
==========

The ``lignoforge.simulation`` sub-package implements the hierarchical
Markov-Chain Monte Carlo (hMMC) algorithm that drives lignin chain
generation.

.. contents:: Contents
   :local:
   :depth: 2

Trajectory
----------

.. autoclass:: lignoforge.simulation.trajectory.Trajectory
   :members:
   :undoc-members:
   :show-inheritance:

:class:`~lignoforge.simulation.trajectory.Trajectory` runs the inner
Metropolis loop for a **single chain**.  Its primary method is
:meth:`~lignoforge.simulation.trajectory.Trajectory.run_MCMC`, which
accepts a random seed and a maximum step count and returns a tuple::

    (polymer, distance_list, monomer_count, i_step)

where ``distance_list`` records the evolution of the distance metric at
each accepted step — useful for diagnosing convergence.

Ring closures are handled by
:meth:`~lignoforge.simulation.trajectory.Trajectory.run_MCMC_ring`, which
attempts to close rings on an already-grown polymer::

    polymer_with_rings = traj.run_MCMC_ring(polymer, rseed=42)

Metropolis criterion
~~~~~~~~~~~~~~~~~~~~~

The acceptance probability for a proposed move with change in distance
:math:`\Delta d` is:

.. math::

   P_\text{accept} = \begin{cases}
       1 & \text{if } \Delta d \le 0 \\
       \exp\!\left(-\dfrac{\Delta d}{k_B \, T_{\text{metro}}}\right) & \text{if } \Delta d > 0
   \end{cases}

where :math:`k_B = 8.617 \times 10^{-5}` eV K⁻¹ and
:math:`T_{\text{metro}}` is the inner-loop temperature.

Simulation
----------

.. autoclass:: lignoforge.simulation.population.Simulation
   :members:
   :undoc-members:
   :show-inheritance:

:class:`~lignoforge.simulation.population.Simulation` extends
:class:`~lignoforge.simulation.trajectory.Trajectory` with the outer
population loop that fills a library of ``n_population`` accepted chains::

    from lignoforge.simulation.population import Simulation

    sim = Simulation(
        linkage_distribution_input=[0.0, 0.05, 0.52, 0.12, 0.14, 0.10, 0.07],
        monomer_distribution_input=[0.01, 0.97, 0.02],
        expected_size=18,
        max_size=70,
        n_population=50,
        i_max=2000,
        i_max_out=2000,
        branching_propensity=0.02,
        seed_init=0,
        library_name="pine_library",
        results_name="output",
    )
    sim.run()
    polymers = sim.P_population  # list[Polymer]

Outer-loop distance
~~~~~~~~~~~~~~~~~~~~

At the population level a second Metropolis criterion (governed by
``Tmetro_out``) decides whether a new chain should **replace** an existing
one already in the population.  The outer distance is computed on the
running mean composition of the population.

.. seealso::

   :doc:`/concepts/optimization_algorithm` for a full description of the
   hierarchical Monte Carlo algorithm, distance metric, and convergence
   criteria.
