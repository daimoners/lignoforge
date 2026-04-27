==========================
Priors and Literature Data
==========================

LignoForge estimates structural priors automatically from a curated
literature database, so users do not need to provide every parameter
explicitly.  This page describes the prior estimation strategy, the
reference data tables, and how experimental input overrides defaults.

What are priors?
----------------

"Priors" in the LignoForge context are the set of parameters required by
the kMC simulation engine:

* Monomer molar fractions (S, G, H)
* Linkage molar fractions (7 types)
* Number-average molecular weight M\ :sub:`n`
* Weight-average molecular weight M\ :sub:`w`
* Polydispersity index PDI
* Mean and maximum degree of polymerisation (``mean_DP``, ``max_DP``)
* Branching index
* Condensation degree (C–C fraction among inter-monomer bonds)
* Composition-weighted average monomer molecular weight

The user may provide any subset of these; values not specified are
inferred from the literature database.

Prior estimation strategy
--------------------------

:class:`~lignoforge.priors.estimator.LigninPriorEstimator` follows a
three-level precedence rule for each parameter:

1. **Direct input** — if the user supplies the value in the JSON input
   (e.g. ``"S_fraction": 0.55``), it is used as-is (clipped to valid
   range if necessary).
2. **Constrained prior** — if partial information is available (e.g.
   biomass type and extraction process), the corresponding literature
   distribution is looked up and a truncated-normal sample is drawn.
3. **Unconditional prior** — if no relevant information is available, a
   sample is drawn from the global (process-independent) distribution.

All sampling uses a seeded :class:`numpy.random.Generator` so results are
fully reproducible.

Literature reference tables
-----------------------------

The reference data live in :mod:`lignoforge.priors.literature_data`.  The
module defines five dictionaries:

``SGH_BY_BIOMASS_PROCESS``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Monomer molar fractions indexed by ``(biomass_type, process_type)``.
Each entry is a tuple ``(S_mean, S_std, G_mean, G_std, H_mean, H_std)``.

Supported biomass types: ``hardwood``, ``softwood``,
``agricultural_residue``, ``mixed``.

Supported process types: ``native``, ``kraft``, ``organosolv``,
``soda``, ``sulfite``, ``des``, ``steam_explosion``.

Example entries:

.. code-block:: python

    "hardwood": {
        "native": (0.55, 0.08, 0.43, 0.08, 0.02, 0.01),
        "kraft":  (0.48, 0.08, 0.50, 0.08, 0.02, 0.01),
        ...
    }

``LINKAGE_BY_BIOMASS_PROCESS``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Linkage molar fractions indexed by ``(biomass_type, process_type,
linkage_name)``.  Each entry is a ``(mean, std)`` tuple.

The 7 linkage types are: ``beta-O-4``, ``4-O-5``, ``alpha-O-4``,
``5-5``, ``beta-5``, ``beta-beta``, ``beta-1``.

The table reflects the well-known correlation between extraction severity
and C–C condensation: kraft processing reduces β-O-4 from ~60 % to ~40 %
in hardwood and increases 5-5 and β-5 accordingly.

``MW_BY_PROCESS``
~~~~~~~~~~~~~~~~~

Number-average M\ :sub:`n` and weight-average M\ :sub:`w` (g/mol) as
``(mean, std)`` tuples indexed by process type.  PDI is derived as
M\ :sub:`w` / M\ :sub:`n`.

``BRANCHING_BY_PROCESS``
~~~~~~~~~~~~~~~~~~~~~~~~~

Branching index (fraction of monomers with degree ≥ 3 in ``bigG``)
indexed by process type.

``CONDENSATION_BY_PROCESS``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Condensation degree (C–C linkage fraction) indexed by process type.

Monomer molecular weights
~~~~~~~~~~~~~~~~~~~~~~~~~~

``MONOMER_MW`` gives the average C9 unit molecular weight for each type::

    MONOMER_MW = {"H": 150.17, "G": 180.20, "S": 210.23}

These are used to compute ``avg_monomer_MW`` from the composition and to
convert between M\ :sub:`n` / M\ :sub:`w` and mean/max DP.

Overriding priors manually
---------------------------

You can inspect or override the estimated priors before running the
simulation::

    pipeline = LigninPipeline.from_json("input.json")
    priors   = pipeline.estimate_priors(random_seed=42)

    # Override a single prior
    priors["branching_index"] = 0.30

    # Re-translate with the modified priors
    pipeline.priors = priors
    pipeline.translate_parameters()

    results = pipeline.run()

Alternatively, supply ``simulation_overrides`` to ``run()`` to override
the translated simulation kwargs directly without touching the priors::

    results = pipeline.run(
        simulation_overrides={
            "branching_propensity": 0.30,
            "n_population":         80,
        }
    )

Data sources
-------------

The reference tables in :mod:`lignoforge.priors.literature_data` are
compiled from the following primary literature (selected):

* Rinaldi *et al.*, *Angew. Chem.* **2016** — linkage distribution overview
* Ragauskas *et al.*, *Science* **2014** — biomass composition
* Chio *et al.*, *Bioresour. Technol.* **2019** — kraft / organosolv contrast
* Constant *et al.*, *Green Chem.* **2016** — technical lignin MW
* Lupoi *et al.*, *Front. Bioeng. Biotechnol.* **2015** — S/G ratios
