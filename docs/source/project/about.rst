==================
About LignoForge
==================

LignoForge is a Python framework for stochastic generation of lignin
polymer models from experimental constraints.  It is developed by the
**DAIMON** research group.

Design goals
------------

LignoForge was designed around four principles:

1. **Reproducibility** — every run is seeded and the full parameter set
   is serialised alongside outputs.
2. **Scientific accessibility** — domain scientists can drive the entire
   workflow through a JSON file without writing Python code.
3. **Modularity** — each sub-package (``core``, ``simulation``,
   ``pipeline``, ``priors``, ``io``, ``structure``) has well-defined
   interfaces and can be used independently.
4. **Compatibility** — outputs follow standard cheminformatics conventions
   (PDB residue codes, SMILES, SDF) for immediate use in downstream
   molecular-dynamics or cheminformatics workflows.

Citation
--------

If you use LignoForge in academic work please cite:

.. code-block:: text

    LignoForge v0.2.0 — DAIMON research group, 2026.
    https://github.com/daimoners/lignoforge

Licence
-------

LignoForge is released under the MIT Licence.  See the ``LICENSE`` file
in the repository root for the complete text.

Acknowledgements
-----------------

The prior estimation literature tables are compiled from contributions
by Rinaldi *et al.*, Ragauskas *et al.*, Chio *et al.*, Constant *et al.*,
and Lupoi *et al.*  Full references are listed in
:doc:`/concepts/priors_and_literature`.

