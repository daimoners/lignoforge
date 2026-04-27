=======================
Command-Line Tools
=======================

LignoForge provides the ``lignoforge-chain`` command-line tool for building
individual lignin chains directly from a terminal, without writing any Python
code.  It is registered as a script entry point during installation, so it is
available anywhere in the environment after ``pip install lignoforge`` (or
``pip install -e .`` for a development install).

.. note::

   For generating a *population* of chains with full statistical analysis, use
   :doc:`pipeline_usage` or ``examples/demo_run.py`` instead.
   ``lignoforge-chain`` is designed for single-chain or small-batch runs, but
   with ``--format all --population-stats --export-priors --export-sim-params``
   it can produce the same complete artefact set as the full pipeline demo.


``lignoforge-chain``
====================

Synopsis
--------

.. code-block:: text

   lignoforge-chain [INPUT_JSON] [options]

   # or, without installation:
   python -m lignoforge.cli.build_chain [INPUT_JSON] [options]

``INPUT_JSON`` is the optional path to a :doc:`LignoForge JSON input file
<input_format>` that supplies the biomass type and extraction process context.
When omitted, built-in hardwood-kraft defaults are used automatically.

Quick Examples
--------------

Minimal run — 10-monomer chain, SMILES + PDB output::

   lignoforge-chain examples/lignin_input_example.json --n-monomers 10

No JSON required — uses built-in defaults::

   lignoforge-chain --n-monomers 10

Target a molecular weight (g/mol) instead of a monomer count::

   lignoforge-chain input.json --mw-target 3000

Pure guaiacyl softwood chain, β-O-4-rich, 3 independent chains with all
output formats, using 4 parallel workers::

   lignoforge-chain input.json \
       --monomer-type G --beta-O-4 0.80 \
       --n-chains 3 --format all --n-workers 4 --seed 7

Branched chain, no 3-D generation (fast, SMILES only)::

   lignoforge-chain input.json \
       --n-monomers 20 --branching 0.05 \
       --format smiles --no-3d

Demo-equivalent command (same artefacts as ``python examples/demo_run.py``)::

   lignoforge-chain examples/lignin_input_example.json \
       --n-chains 2 --seed 42 \
       --format all \
       --population-stats --export-priors --export-sim-params \
       --n-workers 4 --output demo_cli/

Option Reference
----------------

Chain size
~~~~~~~~~~

These two options are mutually exclusive.

.. list-table::
   :header-rows: 1
   :widths: 28 65

   * - Option
     - Description
   * - ``--n-monomers N``
     - Target degree of polymerisation (number of monomeric units).  Uses a
       direct growth loop that produces **exactly** N monomers, bypassing the
       stochastic Trajectory engine.
   * - ``--mw-target MW``
     - Approximate target molecular weight in g/mol.  Converted to a monomer
       count using the composition-weighted average monomer MW minus the
       condensation loss (~18 g/mol per bond), then grown exactly.

Monomer composition
~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 65

   * - Option
     - Description
   * - ``--monomer-type T``
     - Force a single monomer type: ``H``, ``G``, or ``S``.  Use ``mix``
       (default) to keep the composition inferred from the prior.
   * - ``--G-fraction F``
     - Override the guaiacyl fraction (0–1).  Re-normalised with the other
       fractions.
   * - ``--S-fraction F``
     - Override the syringyl fraction (0–1).
   * - ``--H-fraction F``
     - Override the p-hydroxyphenyl fraction (0–1).

Linkage distribution
~~~~~~~~~~~~~~~~~~~~

Individual linkage fractions can be set independently; unspecified fractions
are taken from the estimated prior, and the whole vector is re-normalised
before use.

.. list-table::
   :header-rows: 1
   :widths: 20 73

   * - Option
     - Linkage
   * - ``--beta-O-4 F``
     - β-O-4 ether
   * - ``--alpha-O-4 F``
     - α-O-4 ether
   * - ``--4-O-5 F``
     - 4-O-5 diaryl ether
   * - ``--beta-5 F``
     - β-5 phenylcoumaran
   * - ``--5-5 F``
     - 5-5 biphenyl
   * - ``--beta-beta F``
     - β-β resinol
   * - ``--beta-1 F``
     - β-1 spirodienone

Branching & topology
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 65

   * - Option
     - Description
   * - ``--branching F``
     - Propensity of a branch point per MC step (0–1).  ``0`` (default)
       gives a linear chain; values around 0.02–0.05 produce lightly
       branched chains typical of technical lignins.

Simulation
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 65

   * - Option
     - Description
   * - ``--seed N``
     - Integer random seed for full reproducibility (default: 42).
   * - ``--n-chains N``
     - Number of independent chains to generate (default: 1).
       Each chain receives seed ``N + i`` so results are reproducible
       even when generating multiple chains in one call.
   * - ``--Tmetro T``
     - Inner Metropolis temperature for the chain-growth MCMC
       (default: from priors / input JSON).
   * - ``--max-steps N``
     - Maximum MC steps per chain (default: 2000).  Set higher for
       very large chains or highly constrained linkage distributions.
   * - ``--n-workers N``
     - Number of parallel worker processes for 3-D coordinate embedding
       (default: all logical CPU cores).  Set to ``1`` to force sequential
       execution, which is useful for debugging or very small runs.

Output
~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 63

   * - Option
     - Description
   * - ``--output DIR``
     - Output directory (default: ``chain_output/``).  Created
       automatically if it does not exist.
   * - ``--name NAME``
     - Base name for all output files (default: ``chain``).
   * - ``--format FMT``
     - Comma-separated list of output formats (default: ``smiles,pdb``).
       See :ref:`cli-formats` for the full list.
   * - ``--no-3d``
     - Skip 3-D coordinate generation entirely.  All formats that require
       coordinates (``pdb``, ``pdb-cg``, ``sdf``, ``json-atomistic``,
       ``json-cg``, ``html``) are silently skipped.
   * - ``--no-optimize``
     - Embed 3-D coordinates with RDKit ETKDGv3 but skip MMFF94/UFF
       geometry optimisation.  Faster; geometry is less accurate.
   * - ``--no-H``
     - Strip explicit hydrogens from 3-D structures (PDB, SDF).
   * - ``--max-iter N``
     - Maximum MMFF94/UFF geometry optimisation iterations per chain
       (default: 300).
   * - ``--export-priors``
     - Write ``<name>_estimated_priors.json`` containing the structural
       priors estimated from the literature database for the given input.
   * - ``--export-sim-params``
     - Write ``<name>_simulation_parameters.json`` with the translated
       simulation kwargs passed to the chain-growth engine.
   * - ``--population-stats``
     - Write ``<name>_population_stats.json`` with ensemble-level
       statistics (mean, std, min, max) for every numeric property.
   * - ``--verbose``
     - Print per-step chain-growth progress to stdout.

.. _cli-formats:

Available Formats
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 24 69

   * - Format token
     - Output
   * - ``smiles``
     - One ``.smi`` file per chain containing the canonical SMILES string.
   * - ``population-smiles``
     - A single ``<name>_population.smi`` with all chains (SMILES +
       identifier per line), compatible with RDKit batch processing.
   * - ``sdf``
     - One SDF structure file per chain (3-D coordinates if enabled).
   * - ``pdb``
     - Atomistic PDB files in ``pdb/`` subfolder.  Requires 3-D.
   * - ``pdb-cg``
     - Coarse-grained pseudo-atom PDB files (one bead per monomer) in
       ``pdb_cg/`` subfolder.  Coordinates are monomer centres-of-mass
       derived from the atomistic topology.  Requires 3-D.
   * - ``json-atomistic``
     - Population-level atomistic topology JSON
       (``<name>_atomistic_topology.json``), containing all atom
       coordinates, bonds, and metadata for the entire run.  Requires 3-D.
   * - ``json-cg``
     - Population-level coarse-grained topology JSON
       (``<name>_cg_topology.json``).  One bead per monomer, with
       3-D coordinates and inter-monomer bond list.  Requires 3-D.
   * - ``html``
     - Self-contained interactive Plotly 3-D viewer
       (``<name>_cg_viewer.html``) showing all CG beads coloured by
       monomer type and bonds coloured by linkage.  Open in any browser.
       Requires ``plotly`` and 3-D.
   * - ``all``
     - All of the above.

Output Files
------------

For a run with ``--name chain``, ``--n-chains 3``, and ``--format all`` the
output directory will contain::

   chain_output/
   ├── chain.smi                        # chain 0 — SMILES
   ├── chain_001.smi                    # chain 1
   ├── chain_002.smi                    # chain 2
   ├── chain_population.smi             # all chains combined
   ├── chain.sdf                        # chain 0 — SDF
   ├── chain_001.sdf
   ├── chain_002.sdf
   ├── pdb/
   │   ├── chain_0.pdb                  # atomistic PDB, chain 0
   │   ├── chain_1.pdb
   │   └── chain_2.pdb
   ├── pdb_cg/
   │   ├── chain_cg_0.pdb               # coarse-grained PDB
   │   ├── chain_cg_1.pdb
   │   └── chain_cg_2.pdb
   ├── chain_atomistic_topology.json    # population-level atomistic JSON
   ├── chain_cg_topology.json           # population-level CG JSON
   ├── chain_cg_viewer.html             # interactive 3-D viewer
   ├── chain_chain_stats.json           # per-chain characterisation
   ├── chain_population_stats.json      # ensemble statistics (--population-stats)
   ├── chain_estimated_priors.json      # prior estimates  (--export-priors)
   ├── chain_simulation_parameters.json # sim kwargs       (--export-sim-params)
   └── chain_manifest.json              # full run record

``<name>_chain_stats.json``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A JSON array where each element corresponds to one generated chain.
Fields are identical to those returned by
:meth:`~lignoforge.core.characterization.Characterize.summary`:

.. code-block:: json

   [
     {
       "monomer_count": 10,
       "H_count": 0, "G_count": 10, "S_count": 0,
       "branching_coeff": 0.1,
       "OCH3_count": 10.0,
       "OH_count": 11.0,
       "MW": 1839.9,
       "smiles": "...",
       "linkage_beta-O-4": 5,
       "linkage_5-5": 1,
       "chain_index": 0,
       "seed": 42
     }
   ]

``<name>_manifest.json``
~~~~~~~~~~~~~~~~~~~~~~~~

Records all run parameters exactly so that the same chain can be fully
reproduced at any time:

.. code-block:: json

   {
     "input_file": "/path/to/input.json",
     "target_dp": 10,
     "monomer_dist": {"H": 0.0, "G": 1.0, "S": 0.0},
     "linkage_dist": {"beta-O-4": 0.80},
     "branching": 0.0,
     "seed": 42,
     "n_chains": 3,
     "Tmetro": 298.15,
     "generate_3d": true,
     "optimize_3d": true,
     "n_workers": null,
     "formats": ["pdb", "smiles"]
   }

Programmatic API
----------------

``lignoforge-chain`` uses two internal code paths depending on whether an
explicit size is given:

- **Explicit size** (``--n-monomers`` / ``--mw-target``): a direct growth loop
  via :class:`~lignoforge.core.polymer.Polymer` that produces **exactly** N
  monomers regardless of MC stochasticity.
- **No size flag**: the stochastic :class:`~lignoforge.simulation.trajectory.Trajectory`
  engine, which draws a target DP from a truncated normal prior and runs an
  MCMC loop, giving a realistic MW distribution across chains.

Both paths are accessible directly from Python::

   # --- Exact N monomers via Polymer API ---
   from lignoforge.core.monomer import Monomer
   from lignoforge.core.polymer import Polymer
   from lignoforge.core.utils import (
       set_random_state, generate_random_monomer, generate_random_linkage
   )
   import numpy as np

   rstate = set_random_state(42)
   m_dist = np.array([0.0, 1.0, 0.0])   # pure G
   l_dist = np.array([0.0, 0.0, 0.80, 0.07, 0.06, 0.07, 0.0])   # beta-O-4 rich

   polymer = Polymer(Monomer(generate_random_monomer(m_dist, rstate)))
   for _ in range(9):   # adds 9 more monomers → total 10
       polymer.add_random_monomer(m_dist, l_dist, random_state=rstate)

   # --- Stochastic growth via Trajectory ---
   from lignoforge.simulation.trajectory import Trajectory

   traj = Trajectory(
       linkage_distribution_input=l_dist.tolist(),
       monomer_distribution_input=m_dist.tolist(),
       expected_size=10.0,
       max_size=15.0,
       distribution_scaling=1.0,
       Tmetro=298.15,
   )
   polymer, distances, n_monomers, n_steps = traj.run_MCMC(rseed=42, i_max=2000)

The returned ``polymer`` is a :class:`~lignoforge.core.polymer.Polymer`
object that can be passed to
:meth:`~lignoforge.structure.generator.MolecularStructureGenerator.polymer_atomistic_topology`
or :meth:`~lignoforge.core.characterization.Characterize.summary` for
further analysis.
