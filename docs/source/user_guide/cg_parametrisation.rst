.. _cg_parametrisation:

=====================================
CG parametrisation of bulk lignin
=====================================

``lignoforge-cg-param`` (and the *CG parametrisation* page of the web
interface) derive a coarse-grained model of bulk lignin **from atomistic
simulations that the tool itself builds and runs** with GROMACS.  The package
is :mod:`lignoforge.cgfit`; it ships only the tools, all simulation data live
in a work directory of your choice.

Method
======

1. **Reference.**  An ensemble of chains (several distinct sequences, several
   copies each) is packed, minimised, compressed, equilibrated and simulated
   with OPLS-AA at the target temperature (default 500 K, i.e. a melt).
2. **Mapping.**  Every atom is assigned to a bead (centre of mass):

   * ``ring-side`` (default): two beads per monomer, the aromatic ring with its
     methoxyls and phenolic oxygen (``RH``/``RG``/``RS``) and the propyl side
     chain (``SC``);
   * ``monomer``: one bead per monomer (``MH``/``MG``/``MS``).

   The bead graph (bonds labelled by linkage, angles, 1-2/1-3 exclusions) is
   derived automatically from the atomistic topology.
3. **Bonded terms**: Boltzmann inversion of the mapped bond and angle
   distributions, harmonic wells, per linkage class with a bead-type fallback.
4. **Non-bonded terms**: 12-6 Lennard-Jones per bead-type pair (potential-shifted
   at the cut-off), optimised by *relative-entropy minimisation*: CG runs at the
   reference density are compared with the mapped reference through the
   non-excluded pair-distance histograms; a damped, trust-region Newton step
   updates ``sigma`` and ``epsilon``; a pressure step rescales all ``epsilon`` so
   that the CG pressure approaches the reference pressure.  ``mixing = "lb"``
   fits per-type parameters combined with Lorentz-Berthelot rules, which
   reduces the number of parameters and gives cross terms to other species.
5. **Validation**: a constant-pressure CG run is compared with the reference
   (density, :math:`g(r)`, chain radius of gyration) against the tolerances in
   ``config.json``.

The CG model is state-dependent: parameters are valid at the temperature (and
composition) they were fitted at.  It is intended for equilibration, chain
size and packing/contact properties, not for glass-transition or modulus
predictions.

Command line
============

.. code-block:: bash

   lignoforge-cg-param init   runs/cg1 [--size small|medium|large] [--temperature 500] [--gmx PATH]
   lignoforge-cg-param status runs/cg1 [--follow]      # progress and time remaining
   lignoforge-cg-param log    runs/cg1 [--md] [-f]
   lignoforge-cg-param run    runs/cg1   # resumable; --stage a,b to select stages
   lignoforge-cg-param slurm  runs/cg1 --gres gpu:1 --time 48:00:00   # writes submit.sh
   lignoforge-cg-param selftest runs/selftest   # optimiser check on a synthetic reference

``init`` writes a complete ``config.json`` for the machine: it finds ``gmx``,
detects GPU support (``--gpu auto|yes|no``), sets ``mdrun_args`` and threads
(``--threads``), and prints the estimated system (chains, atoms, box) with
warnings if the box is too small.  Other options: ``--scheme``, ``--mixing``,
``--monomers``, ``--chains``, ``--copies``, ``--H/--G/--S`` (monomer fractions),
``--compress-ns/--equil-ns/--prod-ns``, ``--iterations``, ``--iter-ns``,
``--cutoff``, ``--seed``, ``--force``.  The file can still be edited by hand
(linkage mix, chemical modification, tolerances).

``selftest`` generates a reference by CG MD with *known* parameters and checks
that the optimiser recovers them from a perturbed start (minutes, no atomistic
run).  ``init --quick`` followed by ``run`` is a minutes-long smoke test of the
whole chain of tools; its numbers are not physical.

Outputs (``<workdir>/result/``): ``cg_params.json`` (all fitted parameters),
``cg_bulk.top/.gro`` and ``cg_npt.mdp`` (a CG bulk system ready for ``grompp``),
``report.md``.  Intermediate data: ``reference/`` (atomistic run),
``cg/history.json`` (per-iteration RMSE of :math:`g(r)`, pressure, parameters),
``validation/validation.json``.

GPU arguments are written by ``init`` when the ``gmx`` found has GPU support; to set them by hand use ``atomistic.mdrun_args`` / ``cg.mdrun_args``, e.g. ``["-ntmpi", "1", "-ntomp", "8", "-nb", "gpu", "-pme", "gpu",
"-update", "gpu"]`` (atomistic) and ``["-ntmpi", "1", "-ntomp", "8", "-nb", "gpu"]`` (CG).

Web interface
=============

*CG parametrisation* in the project sidebar lists jobs, starts new ones (form
with a quick-test preset and an advanced JSON editor) and shows stage status,
convergence of :math:`g(r)` and pressure, validation checks, the fitted
parameters, the live log and downloads of the result files.  The same jobs are
available through ``Project.start_cgparam`` and the ``/cgparam`` REST endpoints.
