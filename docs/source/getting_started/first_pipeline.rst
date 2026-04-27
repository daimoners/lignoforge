============================
Your First Pipeline Run
============================

This walkthrough takes you through a complete end-to-end run step by step,
explaining every stage of the pipeline and matching each code block to the
corresponding module in the package.

Scenario
--------

We generate a structural library for a **kraft-process hardwood lignin**
(poplar) with NMR-derived experimental data as input.


Step 1 — Prepare the input file
---------------------------------

LignoForge accepts a JSON file that describes the sample.  Create
``poplar_kraft.json``::

    {
      "material_origin": {
        "lignin_type": "kraft_lignin",
        "biomass_type": "hardwood"
      },
      "extraction_process": {
        "process_type": "kraft"
      },
      "structural_characterization": {
        "monomer_composition": {
          "S_fraction": 0.55,
          "G_fraction": 0.43,
          "H_fraction": 0.02
        },
        "linkage_distribution": {
          "beta_O_4":   0.48,
          "alpha_O_4":  0.05,
          "4_O_5":      0.05,
          "5_5":        0.10,
          "beta_5":     0.12,
          "beta_beta":  0.11,
          "beta_1":     0.09
        },
        "molecular_weight": {
          "Mn_g_per_mol": 3200,
          "Mw_g_per_mol": 7500
        }
      },
      "simulation_config": {
        "n_population": 30,
        "Tmetro":       10.0,
        "Tmetro_out":   20.0
      }
    }

Keys not provided by you are filled automatically from the literature
database in :mod:`lignoforge.priors.literature_data`.
See :doc:`../user_guide/input_format` for the complete schema reference.

Step 2 — Validate and create the pipeline
------------------------------------------

::

    from lignoforge.pipeline import LigninPipeline

    pipeline = LigninPipeline.from_json(
        "poplar_kraft.json",
        output_dir="output/poplar_kraft",
    )

:meth:`~lignoforge.pipeline.pipeline.LigninPipeline.from_json` calls
:class:`~lignoforge.io.schema.InputSchemaValidator` internally.  If
validation fails you will receive a ``jsonschema.ValidationError`` with a
clear message identifying the offending key.

Step 3 — Estimate structural priors
-------------------------------------

::

    priors = pipeline.estimate_priors(random_seed=42)
    print(priors)

:class:`~lignoforge.priors.estimator.LigninPriorEstimator` produces a flat
dictionary with these keys:

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Key
     - Description
   * - ``S_fraction``, ``G_fraction``, ``H_fraction``
     - Monomer molar fractions (sum ≈ 1)
   * - ``linkage_fractions``
     - Dict ``{linkage_name: fraction}`` (sum = 1)
   * - ``Mn``, ``Mw``, ``PDI``
     - Number-average MW, weight-average MW, dispersity
   * - ``mean_DP``, ``max_DP``
     - Degree of polymerisation (monomer count)
   * - ``branching_index``
     - Fraction of monomers with ≥ 3 inter-monomer bonds
   * - ``condensation_degree``
     - C–C bond fraction among all inter-monomer bonds
   * - ``avg_monomer_MW``
     - Composition-weighted monomer molecular weight

Step 4 — Translate priors to simulation parameters
----------------------------------------------------

::

    sim_kwargs = pipeline.translate_parameters()

:class:`~lignoforge.pipeline.translator.ParameterTranslator` converts
the priors dictionary into the exact keyword arguments expected by
:class:`~lignoforge.simulation.population.Simulation`:

* ``linkage_distribution_input`` — normalised 7-element probability vector
* ``monomer_distribution_input`` — normalised [H, G, S] probability vector
* ``expected_size`` / ``max_size`` — drawn from ``mean_DP`` / ``max_DP``
* ``distribution_scaling`` — derived from PDI
* ``Tmetro`` / ``Tmetro_out`` — Metropolis temperatures
* ``branching_propensity`` — from ``branching_index``

Step 5 — Run the population simulation
----------------------------------------

::

    results = pipeline.run(
        random_seed=42,
        library_name="poplar_kraft",
    )

The call executes the three MCMC loops (see :doc:`../concepts/optimization_algorithm`)
and populates ``results.polymers`` — a list of
:class:`~lignoforge.core.polymer.Polymer` objects.

Progress is printed to stdout::

    Starting trial No.0  →  output/poplar_kraft/...
    Chain 1 / 30  |  16 monomers  |  d = 0.0041
    Chain 2 / 30  |  18 monomers  |  d = 0.0039
    ...

Step 6 — Inspect results
--------------------------

::

    from lignoforge.core import Characterize

    for i, polymer in enumerate(results.polymers[:3]):
        ch = Characterize(polymer)
        s  = ch.summary()
        print(f"Chain {i}: {s['monomer_count']} monomers, "
              f"MW = {s['MW']:.0f} g/mol, "
              f"β-O-4 = {s['linkage_fracs']['beta-O-4']:.2f}")

Step 7 — Export formats
------------------------

The default pipeline bundle is written automatically.  For additional
formats call :class:`~lignoforge.io.exporters.LigninExporter` directly::

    from lignoforge.io.exporters import LigninExporter

    exp = LigninExporter("output/poplar_kraft")

    # Generate PDB files with 3-D coordinates
    pdb_paths = exp.export_population_pdbs(
        results.polymers,
        folder="pdb",
        explicit_hydrogens=True,
        optimize_3d=True,
    )

    # Generate SDF files for cheminformatics tools
    sdf_path = exp.export_sdf(results.polymers)

See :doc:`../user_guide/output_formats` for a complete description of every
output file.

Step 8 — Reproducibility
--------------------------

All stochastic choices are seeded through a single ``random_seed`` integer.
Setting the same seed on identical input always produces byte-identical output::

    r1 = pipeline.run(random_seed=7)
    r2 = pipeline.run(random_seed=7)

    assert r1.artifacts == r2.artifacts   # same files produced
