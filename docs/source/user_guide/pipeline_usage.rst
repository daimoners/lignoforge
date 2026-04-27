==============
Pipeline Usage
==============

This guide covers all ways to drive the LignoForge pipeline, from the
simplest one-liner to fine-grained step-by-step control.

Creating a pipeline
--------------------

From a JSON file (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    from lignoforge.pipeline import LigninPipeline

    pipeline = LigninPipeline.from_json(
        "my_lignin.json",
        output_dir="results/my_run",
    )

From a Python dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~

::

    data = {
        "material_origin": {"lignin_type": "kraft_lignin",
                            "biomass_type": "hardwood"},
        "extraction_process": {"process_type": "kraft"},
    }

    pipeline = LigninPipeline.from_dict(data, output_dir="results/")

Pass ``validate=False`` to skip schema validation (e.g. in tests)::

    pipeline = LigninPipeline.from_dict(data, validate=False)

Running the full population
----------------------------

::

    results = pipeline.run(
        random_seed=42,
        library_name="poplar_library",
    )

``results`` is a :class:`~lignoforge.pipeline.pipeline.LigninResults`
dataclass with these attributes:

.. list-table::
   :header-rows: 1
   :widths: 28 60

   * - Attribute
     - Contents
   * - ``polymers``
     - ``list[Polymer]`` — all accepted chains
   * - ``priors``
     - Estimated prior dictionary
   * - ``simulation_kwargs``
     - Translated simulation parameters
   * - ``input_data``
     - Original validated input dict
   * - ``output_dir``
     - Absolute path to the output directory
   * - ``artifacts``
     - ``dict {name: path}`` of all written files

Generating a single chain (fast)
----------------------------------

Use :meth:`~lignoforge.pipeline.pipeline.LigninPipeline.run_single_chain`
when you only need one chain (e.g. for debugging or a quick test)::

    results = pipeline.run_single_chain(
        random_seed=7,
        simulation_overrides={"Tmetro": 50.0},
    )

This forces ``n_population=1``, ``i_max=40``, ``i_max_out=40``,
``i_max_ring=0`` unless overridden.

Overriding simulation parameters
----------------------------------

Pass ``simulation_overrides`` to ``run()`` to change translated parameters
without modifying the priors::

    results = pipeline.run(
        random_seed=42,
        simulation_overrides={
            "n_population":       100,
            "expected_size":      25,
            "branching_propensity": 0.25,
            "Tmetro":             5.0,
            "Tmetro_out":         10.0,
        },
    )

Step-by-step control
---------------------

Calling each pipeline stage explicitly gives you access to the
intermediate results::

    # 1. Estimate priors
    priors = pipeline.estimate_priors(random_seed=42)
    print(f"Estimated Mn = {priors['Mn']:.0f} g/mol")

    # 2. (Optional) adjust priors manually
    priors["branching_index"] = 0.20

    # 3. Translate to simulation kwargs
    pipeline.priors = priors
    sim_kwargs = pipeline.translate_parameters()

    # 4. Run simulation
    results = pipeline.run(random_seed=42)

Exporting PDB files after the run
-----------------------------------

PDB export is optional (expensive for large populations) and can be
triggered after ``run()``::

    # Export PDB for every chain
    pdb_paths = pipeline.export_population_pdbs(
        folder="pdb",
        explicit_hydrogens=True,
        optimize_3d=True,
        max_uff_iterations=200,
    )

    # Export PDB for a single chain only
    pdb_path = pipeline.export_single_pdb(
        filename="best_chain.pdb",
        optimize_3d=True,
    )

Enabling 3-D export options
-----------------------------

Export options (3-D generation, hydrogen inclusion, geometry optimisation,
and parallelism) are passed inside ``simulation_overrides``.  The pipeline
extracts them automatically before forwarding the remaining keys to the
simulator::

    results = pipeline.run(
        random_seed=42,
        simulation_overrides={
            "n_population":       20,
            "generate_3d":        True,
            "include_hydrogens":  True,
            "optimize_3d":        True,
            "max_uff_iterations": 300,
            "n_workers":          None,   # None = all CPU cores
        },
    )

To disable 3-D generation entirely (topology and SMILES are still
exported)::

    results = pipeline.run(
        simulation_overrides={"generate_3d": False}
    )

See :ref:`export-options-reference` for the full list of options.

Reproducibility
---------------

Every stochastic operation in LignoForge is seeded through the
``random_seed`` argument.  Using identical seeds on identical input
always produces identical results::

    r1 = pipeline.run(random_seed=99)
    r2 = pipeline.run(random_seed=99)

    import json
    assert (json.load(open(r1.artifacts["chain_statistics"])) ==
            json.load(open(r2.artifacts["chain_statistics"])))
