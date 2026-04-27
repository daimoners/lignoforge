===========
Quick Start
===========

This page shows how to get meaningful results in under five minutes,
assuming LignoForge is already installed (see :doc:`installation`).

Minimal pipeline run
---------------------

The following three lines run the full pipeline — input validation, prior
estimation, population simulation, and export — using the bundled example
JSON::

    from lignoforge.pipeline import LigninPipeline

    pipeline = LigninPipeline.from_json("examples/lignin_input_example.json")
    results  = pipeline.run(library_name="quickstart")

``results`` is a :class:`~lignoforge.pipeline.pipeline.LigninResults`
dataclass.  Print the produced files::

    for name, path in results.artifacts.items():
        print(f"  {name:30s} {path}")

Typical output::

  input_data                     quickstart_results/input_high_level.json
  priors                         quickstart_results/estimated_priors.json
  simulation_parameters          quickstart_results/simulation_parameters.json
  chain_statistics               quickstart_results/chain_statistics.json
  population_statistics          quickstart_results/population_statistics.json
  atomistic_topology             quickstart_results/atomistic_topology.json
  coarse_grained_topology        quickstart_results/coarse_grained_topology.json
  coarse_grained_viewer          quickstart_results/coarse_grained_topology_viewer.html
  coarse_grained_pdb             quickstart_results/pdb_cg/
  pdb                            quickstart_results/pdb/
  smiles                         quickstart_results/population.smi

Command-line demo
------------------

Run the full population pipeline via the demo script::

    python examples/demo_run.py --n-chains 5 --seed 42

Optional flags:

.. list-table::
   :header-rows: 1
   :widths: 30 55

   * - Flag
     - Description
   * - ``--n-chains INT``
     - Number of chains to generate (default: 2)
   * - ``--seed INT``
     - Random seed for full reproducibility (default: 42)
   * - ``--workers INT``
     - Number of parallel worker processes for 3-D embedding
       (default: all CPU cores)
   * - ``--output PATH``
     - Output directory (default: ``demo_output/``)

Or build a single chain with full morphology control using
:doc:`/user_guide/cli_tools`::

    # 10-monomer pure-G chain, β-O-4-rich
    lignoforge-chain examples/lignin_input_example.json \
        --n-monomers 10 --monomer-type G --beta-O-4 0.80

    # Target MW, all output formats
    lignoforge-chain examples/lignin_input_example.json \
        --mw-target 3000 --format all

Inspecting the results
-----------------------

Load chain statistics::

    import json

    with open("demo_output/chain_statistics.json") as f:
        chains = json.load(f)

    for c in chains[:3]:
        print(c["monomer_count"], c["MW"], c["branching_coeff"])

Load population statistics::

    with open("demo_output/population_statistics.json") as f:
        stats = json.load(f)

    print(stats["n_polymers"])
    print(stats["numeric"]["MW"]["mean"])

Open the interactive chain viewer by opening
``demo_output/coarse_grained_topology_viewer.html`` in any modern browser.

What to read next
-----------------

* :doc:`first_pipeline` — detailed walkthrough with a realistic input file
* :doc:`../concepts/index` — how the graph model and algorithm work
* :doc:`../user_guide/input_format` — complete JSON input specification
