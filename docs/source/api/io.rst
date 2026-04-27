==
IO
==

The ``lignoforge.io`` sub-package handles input validation and output
export.

.. contents:: Contents
   :local:
   :depth: 2

InputSchemaValidator
---------------------

.. autoclass:: lignoforge.io.schema.InputSchemaValidator
   :members:
   :undoc-members:
   :show-inheritance:

The validator checks the user-supplied dictionary against the
``lignin_info_schema.json`` JSON Schema bundled inside the package at
``lignoforge/io/lignin_info_schema.json``.  It is
called automatically by
:meth:`~lignoforge.pipeline.pipeline.LigninPipeline.from_json` and
:meth:`~lignoforge.pipeline.pipeline.LigninPipeline.from_dict`::

    from lignoforge.io.schema import InputSchemaValidator

    validator = InputSchemaValidator(data)
    validator.validate()   # raises jsonschema.ValidationError on failure

LigninExporter
---------------

.. autoclass:: lignoforge.io.exporters.LigninExporter
   :members:
   :undoc-members:
   :show-inheritance:

Convenience bundle
~~~~~~~~~~~~~~~~~~~

Most users should call
:meth:`~lignoforge.io.exporters.LigninExporter.export_full_pipeline_bundle`
which writes every artefact in a single call::

    from lignoforge.io.exporters import LigninExporter

    exporter = LigninExporter(
        polymers=results.polymers,
        output_dir="output",
        library_name="pine_library",
        export_options={"generate_3d": True, "include_hydrogens": True},
    )
    exporter.export_full_pipeline_bundle(results)

Individual exporters
~~~~~~~~~~~~~~~~~~~~~

Each output file can also be generated independently:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Method
     - Output
   * - ``export_smiles``
     - ``population.smi``
   * - ``export_sdf``
     - ``population.sdf``
   * - ``export_summary``
     - (returns dict or writes JSON)
   * - ``export_chain_statistics``
     - ``chain_statistics.json``
   * - ``export_population_statistics``
     - ``population_statistics.json``
   * - ``export_atomistic_topology_json``
     - ``atomistic_topology.json``
   * - ``export_coarse_grained_topology_json``
     - ``coarse_grained_topology.json``
   * - ``export_coarse_grained_viewer_html``
     - ``coarse_grained_topology_viewer.html``

.. seealso::

   :doc:`/user_guide/output_formats` for detailed schema documentation
   for every output file.
