======
Priors
======

The ``lignoforge.priors`` sub-package estimates structural priors from
the literature and from user-supplied constraints.

.. contents:: Contents
   :local:
   :depth: 2

LigninPriorEstimator
---------------------

.. autoclass:: lignoforge.priors.estimator.LigninPriorEstimator
   :members:
   :undoc-members:
   :show-inheritance:

The estimator applies a three-level precedence rule:

1. **Direct input** — any value explicitly provided in the input JSON is
   used as-is.
2. **Constrained prior** — if the biomass type and process are specified,
   the corresponding literature mean is used as a prior for any missing
   quantities.
3. **Unconditional prior** — if neither biomass type nor process is
   provided, a broad average across all entries in the literature
   database is returned.

Usage::

    from lignoforge.priors.estimator import LigninPriorEstimator

    est = LigninPriorEstimator(
        biomass_type="poplar",
        process="kraft",
        S_fraction=0.52,          # override one value
    )
    priors = est.estimate()

Output keys are documented in :doc:`/user_guide/output_formats`.

Literature data tables
-----------------------

.. automodule:: lignoforge.priors.literature_data
   :members:
   :undoc-members:

The module exposes the following lookup tables:

.. list-table::
   :header-rows: 1
   :widths: 35 60

   * - Table
     - Description
   * - ``SGH_BY_BIOMASS_PROCESS``
     - Mean ± std of S, G, H fractions indexed by biomass type × process
   * - ``LINKAGE_BY_BIOMASS_PROCESS``
     - Mean ± std of all seven linkage fractions indexed by process
   * - ``MW_BY_PROCESS``
     - Mn, Mw, PDI statistics per extraction process
   * - ``BRANCHING_BY_PROCESS``
     - Branching index per extraction process
   * - ``CONDENSATION_BY_PROCESS``
     - Condensation degree per extraction process
   * - ``MONOMER_MW``
     - Molecular weight (g/mol) of each monomer type

All values are compiled from:

* Rinaldi *et al.*, *Angew. Chem. Int. Ed.* **55**, 8164 (2016)
* Ragauskas *et al.*, *Science* **344**, 1246843 (2014)
* Chio *et al.*, *Renew. Sustain. Energy Rev.* **107**, 232 (2019)
* Constant *et al.*, *Green Chem.* **18**, 2651 (2016)
* Lupoi *et al.*, *Front. Bioeng. Biotechnol.* **3**, 50 (2015)

.. seealso::

   :doc:`/concepts/priors_and_literature` for a conceptual explanation
   of how priors are estimated and their uncertainty propagated.
