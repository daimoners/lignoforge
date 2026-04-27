========
Pipeline
========

The ``lignoforge.pipeline`` sub-package provides the high-level entry
point for end-to-end lignin generation.

.. contents:: Contents
   :local:
   :depth: 2

LigninPipeline
--------------

.. autoclass:: lignoforge.pipeline.pipeline.LigninPipeline
   :members:
   :undoc-members:
   :show-inheritance:

The most common usage pattern is to load a JSON file and call
:meth:`~lignoforge.pipeline.pipeline.LigninPipeline.run`::

    from lignoforge.pipeline.pipeline import LigninPipeline

    pipeline = LigninPipeline.from_json("poplar_kraft.json")
    pipeline.estimate_priors()
    pipeline.translate_parameters()
    results = pipeline.run()

Alternatively the full workflow can be triggered in a single call::

    results = LigninPipeline.from_json("poplar_kraft.json").run()

LigninResults
-------------

.. autoclass:: lignoforge.pipeline.pipeline.LigninResults
   :members:
   :undoc-members:
   :show-inheritance:

A :class:`~lignoforge.pipeline.pipeline.LigninResults` object is returned
by every call to
:meth:`~lignoforge.pipeline.pipeline.LigninPipeline.run`.  It exposes:

* ``input_data`` — validated input dictionary
* ``priors`` — estimated prior values (same structure as
  ``estimated_priors.json``)
* ``simulation_kwargs`` — keyword arguments passed to
  :class:`~lignoforge.simulation.population.Simulation`
* ``polymers`` — list of accepted :class:`~lignoforge.core.polymer.Polymer`
  instances
* ``output_dir`` — path to the output directory
* ``artifacts`` — list of all written file paths

ParameterTranslator
--------------------

.. autoclass:: lignoforge.pipeline.translator.ParameterTranslator
   :members:
   :undoc-members:
   :show-inheritance:

Translates the flat prior dictionary (e.g., ``S_fraction``,
``linkage_fractions``, ``mean_DP``) into the keyword-argument dictionary
expected by :class:`~lignoforge.simulation.population.Simulation`.  The
main entry point is
:meth:`~lignoforge.pipeline.translator.ParameterTranslator.to_simulation_kwargs`::

    from lignoforge.pipeline.translator import ParameterTranslator

    priors = {
        "S_fraction":  0.45,
        "G_fraction":  0.50,
        "H_fraction":  0.05,
        "linkage_fractions": { ... },
        "mean_DP":     18,
        "max_DP":      70,
        "PDI":         2.2,
        "branching_index": 0.03,
    }
    kwargs = ParameterTranslator.to_simulation_kwargs(priors)
