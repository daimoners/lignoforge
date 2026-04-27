==================
Input File Format
==================

LignoForge accepts a single JSON file as its primary input.  This page
describes every supported key, their types, constraints, and default
behaviour when omitted.

The JSON schema is defined in ``lignoforge/io/lignin_info_schema.json``
(bundled with the package) and is enforced by
:class:`~lignoforge.io.schema.InputSchemaValidator` at pipeline
initialisation.

Top-level structure
--------------------

The input object has four optional top-level sections:

.. code-block:: text

    {
      "material_origin"            :  { ... },   // what the lignin is
      "extraction_process"         :  { ... },   // how it was isolated
      "structural_characterization":  { ... },   // NMR / GPC data
      "simulation_config"          :  { ... }    // simulation hyperparameters
    }

All four sections are optional; the minimal valid input is an empty object
``{}``.

section: ``material_origin``
-----------------------------

.. code-block:: json

    "material_origin": {
        "lignin_type": "kraft_lignin",
        "biomass_type": "hardwood"
    }

.. list-table::
   :header-rows: 1
   :widths: 25 15 55

   * - Key
     - Type
     - Allowed values
   * - ``lignin_type``
     - string
     - ``"kraft_lignin"``, ``"organosolv_lignin"``, ``"soda_lignin"``,
       ``"sulfite_lignin"``, ``"native_lignin"``, ``"des_lignin"``,
       ``"steam_explosion_lignin"``, ``"technical_lignin"``
   * - ``biomass_type``
     - string
     - ``"hardwood"``, ``"softwood"``, ``"agricultural_residue"``,
       ``"mixed"``

Not providing these keys causes the estimator to use the global
(process-independent) literature averages.

section: ``extraction_process``
--------------------------------

.. code-block:: json

    "extraction_process": {
        "process_type": "kraft"
    }

.. list-table::
   :header-rows: 1
   :widths: 25 15 55

   * - Key
     - Type
     - Allowed values
   * - ``process_type``
     - string
     - ``"kraft"``, ``"organosolv"``, ``"soda"``, ``"sulfite"``,
       ``"native"``, ``"des"``, ``"steam_explosion"``

Used jointly with ``biomass_type`` to look up the two-way prior table.

section: ``structural_characterization``
-----------------------------------------

This section holds the experimental NMR / GPC data.  All sub-keys are
optional; any combination may be supplied.

``monomer_composition``
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

    "monomer_composition": {
        "S_fraction": 0.55,
        "G_fraction": 0.43,
        "H_fraction": 0.02
    }

* All values are molar fractions in [0, 1].
* Fractions are normalised to sum to 1 if they do not already.
* Providing only two of the three is allowed; the third is inferred as
  ``1 – (sum of provided)``.

``linkage_distribution``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

    "linkage_distribution": {
        "beta_O_4":   0.48,
        "alpha_O_4":  0.05,
        "4_O_5":      0.05,
        "5_5":        0.10,
        "beta_5":     0.12,
        "beta_beta":  0.11,
        "beta_1":     0.09
    }

* JSON keys use underscores (``beta_O_4``); they map internally to
  hyphenated linkage names (``"beta-O-4"``).
* All values are molar fractions; normalised automatically.
* Any subset of the 7 linkage types may be provided; missing ones are
  filled from literature.

``molecular_weight``
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: json

    "molecular_weight": {
        "Mn_g_per_mol": 3200,
        "Mw_g_per_mol": 7500
    }

* Both ``Mn_g_per_mol`` and ``Mw_g_per_mol`` are in g/mol.
* PDI is derived as Mw / Mn.
* Providing only one is allowed.

section: ``simulation_config``
-------------------------------

Direct overrides for simulation hyperparameters.  Values supplied here
take precedence over the translated priors.

.. code-block:: json

    "simulation_config": {
        "n_population":  30,
        "Tmetro":        10.0,
        "Tmetro_out":    20.0
    }

.. list-table::
   :header-rows: 1
   :widths: 25 12 55

   * - Key
     - Type
     - Description
   * - ``n_population``
     - int
     - Number of chains in the simulated library
   * - ``Tmetro``
     - float
     - Inner-loop Metropolis temperature (K)
   * - ``Tmetro_out``
     - float
     - Outer-loop Metropolis temperature (K)
   * - ``i_max``
     - int
     - Max inner-loop MC steps per chain
   * - ``i_max_out``
     - int
     - Max outer-loop iterations
   * - ``i_max_ring``
     - int
     - Max ring-closure attempts

Minimal valid examples
----------------------

**Only biomass / process type (priors fully inferred)**::

    {
      "material_origin": {"lignin_type": "kraft_lignin",
                          "biomass_type": "hardwood"},
      "extraction_process": {"process_type": "kraft"}
    }

**Full NMR + GPC experimental data**::

    {
      "material_origin": {"lignin_type": "organosolv_lignin",
                          "biomass_type": "softwood"},
      "extraction_process": {"process_type": "organosolv"},
      "structural_characterization": {
        "monomer_composition": {"S_fraction": 0.02,
                                "G_fraction": 0.97,
                                "H_fraction": 0.01},
        "linkage_distribution": {"beta_O_4": 0.52,
                                 "5_5": 0.12,
                                 "beta_5": 0.14,
                                 "beta_beta": 0.10,
                                 "alpha_O_4": 0.05,
                                 "4_O_5": 0.04,
                                 "beta_1": 0.03},
        "molecular_weight": {"Mn_g_per_mol": 4500,
                             "Mw_g_per_mol": 10000}
      },
      "simulation_config": {"n_population": 50}
    }

Programmatic construction and validation
-----------------------------------------

You can construct the input dict in Python and validate it before
creating the pipeline::

    from lignoforge.io.schema import InputSchemaValidator

    data = {
        "material_origin": {"lignin_type": "kraft_lignin",
                            "biomass_type": "hardwood"},
    }

    validator = InputSchemaValidator()
    validator.validate_dict(data)   # raises ValidationError if invalid

    from lignoforge.pipeline import LigninPipeline
    pipeline = LigninPipeline.from_dict(data)
