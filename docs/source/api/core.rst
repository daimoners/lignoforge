====
Core
====

The ``lignoforge.core`` sub-package contains the fundamental building
blocks: chemistry rules, monomer objects, polymer graph data structure,
characterisation metrics, and utility functions.

.. contents:: Contents
   :local:
   :depth: 2

Rules and constants
-------------------

.. automodule:: lignoforge.core.rules
   :members:
   :undoc-members:

Key constants at a glance:

.. list-table::
   :header-rows: 1
   :widths: 35 60

   * - Constant
     - Description
   * - ``monomer_types``
     - List of recognised monomer codes: ``["H", "G", "S"]``
   * - ``MONOMER_RESIDUE_CODE``
     - Mapping to three-letter PDB codes: ``{"H": "HPU", "G": "GYU", "S": "SYU"}``
   * - ``linkage_names``
     - List of seven recognised inter-unit linkage types
   * - ``CHO``
     - Elemental formula arrays *[C, H, O]* per monomer type
   * - ``weight_CHO``
     - Atomic masses of C, H, O in g/mol

Monomer
--------

.. autoclass:: lignoforge.core.monomer.Monomer
   :members:
   :undoc-members:
   :show-inheritance:

Factory functions
~~~~~~~~~~~~~~~~~

.. autofunction:: lignoforge.core.monomer.monomer_H
.. autofunction:: lignoforge.core.monomer.monomer_G
.. autofunction:: lignoforge.core.monomer.monomer_S

Quick reference: monomer atom counts

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 50

   * - Monomer
     - Carbon
     - Oxygen
     - Notes
   * - H
     - 9
     - 2
     - *p*-hydroxyphenyl; no methoxy group
   * - G
     - 10
     - 3
     - guaiacyl; one ring methoxy at C3
   * - S
     - 11
     - 4
     - syringyl; two ring methoxys at C3 and C5

PolymerGraph
------------

.. autoclass:: lignoforge.core.polymer.PolymerGraph
   :members:
   :undoc-members:
   :show-inheritance:

:class:`~lignoforge.core.polymer.PolymerGraph` is the low-level graph
container.  It wraps a :class:`networkx.Graph` (``G``) and a higher-level
monomer graph (``bigG``).  The key node / edge attributes are documented
in :doc:`/concepts/graph_representation`.

Polymer
--------

.. autoclass:: lignoforge.core.polymer.Polymer
   :members:
   :undoc-members:
   :show-inheritance:

:class:`~lignoforge.core.polymer.Polymer` extends
:class:`~lignoforge.core.polymer.PolymerGraph` with a growth API::

    from lignoforge.core.polymer import Polymer
    from lignoforge.core.monomer import monomer_G

    p = Polymer(monomer_G())
    p.add_specific_monomer(monomer_G(), "beta-O-4")
    p.add_random_monomer()     # picks type & linkage stochastically
    p.add_specific_ring(0, 2)  # ring closure between monomers 0 and 2

Characterize
------------

.. autoclass:: lignoforge.core.characterization.Characterize
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: lignoforge.core.characterization.CharacterizeGraph
   :members:
   :undoc-members:
   :show-inheritance:

The :meth:`~lignoforge.core.characterization.Characterize.summary` method
returns a flat dictionary — see :doc:`/user_guide/output_formats` for the
full key listing.

Population
-----------

.. autoclass:: lignoforge.core.characterization.Population
   :members:
   :undoc-members:
   :show-inheritance:

Utility functions
-----------------

.. automodule:: lignoforge.core.utils
   :members:
   :undoc-members:

Highlights:

* :func:`~lignoforge.core.utils.graph_to_smile` — convert the internal
  graph to a SMILES string via *pysmiles*.
* :func:`~lignoforge.core.utils.graph_to_mol` — convert to an RDKit
  :class:`Mol` object with optional 3-D coordinate generation.
* :func:`~lignoforge.core.utils.draw_graph` — 2-D node-level drawing.
* :func:`~lignoforge.core.utils.draw_big_graph` — coarse-grained
  monomer-level drawing.
* :func:`~lignoforge.core.utils.draw_atomic_graph` — full atomic
  representation with element labels.
