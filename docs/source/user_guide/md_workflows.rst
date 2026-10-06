.. _md_workflows:

=====================================
Simulation Workflows (atomistic & CG)
=====================================

``lignoforge-chain`` can write ready-to-run GROMACS inputs for every chain it
builds.  Three formats are available (combine them with commas):

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - Format
     - Output
   * - ``gromacs``
     - ``gromacs/<name>_<i>/``: OPLS-AA ``.top``, ``.gro``, ``em.mdp`` and a
       topology report (see :ref:`force_field`).
   * - ``md``
     - ``md/<name>_<i>/``: full atomistic workflow with ``run.sh``
       (box → solvation → EM → NVT → NPT → production) and ``workflow.json``.
   * - ``cg``
     - ``cg_md/``: one-bead-per-monomer CG system for *all* chains of the run,
       with ``run.sh`` (insert molecules → EM → stochastic dynamics).

.. code-block:: bash

   lignoforge-chain --n-monomers 10 --n-chains 5 --format md,cg \
       --md-temperature 300 --solvent tip3p --prod-ns 20 \
       --cg-copies 20 --cg-density 0.3 --cg-run-ns 200

Atomistic workflow
------------------

``run.sh`` uses only standard GROMACS tools and stops at the first error
(``GMX`` and ``NT`` environment variables select the binary and thread
count).  Timings are set with ``--prod-ns`` and, in Python,
``write_atomistic_workflow(..., nvt_ps=, npt_ps=)``.  ``--solvent none`` gives a
periodic vacuum box with EM → NVT → production (no NPT).  The ``.mdp``
templates (2 fs, h-bond constraints, PME, v-rescale / C-rescale /
Parrinello-Rahman) are plain files that can be edited before running.

Coarse-grained model
--------------------

Each monomer becomes one bead at its centre of mass (mass = monomer mass).
Bonded terms follow the linkage graph:

* harmonic bond per linkage (``r0``, ``k``); an optional type-specific entry
  ``"<linkage>:<A>-<B>"`` overrides it;
* harmonic angle for every pair of linkages sharing a bead, keyed by the
  sorted linkage pair;
* uncharged Lennard-Jones beads (H, G, S) with 1-2 and 1-3 exclusions;
* stochastic dynamics (implicit solvent).

.. warning::

   The shipped parameters are **provisional**.  Bond ``r0`` per linkage and a
   few angles come from an ensemble of MMFF-relaxed 6-mers; force constants,
   bead sizes and well depths are placeholders.  With attractive beads and no
   solvent a single chain collapses.  Derive parameters from atomistic MD
   (below) and tune ``epsilon`` to the solvent quality you want to model.

Deriving CG parameters from atomistic MD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Generate chains and run the atomistic workflows::

      lignoforge-chain --n-monomers 10 --n-chains 4 --format md,json-atomistic
      (cd chain_output/md/chain_0 && ./run.sh)        # repeat for each chain

2. Extract the whole-molecule chain trajectory (the solute atoms come first
   in the system)::

      n=$(awk 'NR==2{print $1}' chain_output/md/chain_0/chain_0.gro)
      echo "a 1-$n
      q" | gmx make_ndx -f chain_output/md/chain_0/chain_0.gro -o chain.ndx
      echo "a_1-$n" | gmx trjconv -s prod.tpr -f prod.xtc -n chain.ndx \
          -pbc mol -o frames0.pdb

3. Boltzmann-invert, pooling all runs::

      lignoforge-cg-fit chain_output/chain_atomistic_topology.json \
          --run 0:frames0.pdb --run 1:frames1.pdb --temperature 300

4. Use the result::

      lignoforge-chain --n-monomers 10 --format cg --cg-params cg_parameters.json

Bonds use :math:`U(r) = -k_BT \ln[P(r)/r^2]` and angles
:math:`U(\theta) = -k_BT \ln[P(\theta)/\sin\theta]`, each fitted with a
harmonic well.  Entries with fewer than ``--min-samples`` samples keep the
defaults and are listed in the output.

Python API
----------

.. code-block:: python

   from lignoforge.md import write_atomistic_workflow
   from lignoforge.cg import write_cg_system, BoltzmannFitter, CGParameters

   write_atomistic_workflow(chain_topology, "run_md", name="lig", solvent="spce")
   write_cg_system([chain_a, chain_b], "run_cg", copies=10,
                   parameters=CGParameters.load("cg_parameters.json"))
