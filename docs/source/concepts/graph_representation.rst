====================
Graph Representation
====================

LignoForge uses **NetworkX undirected graphs** to represent lignin
structures at two resolutions: atomistic and coarse-grained.  This section
describes both representations, their attributes, and the 3-letter residue
code convention used in all output formats.

Why graphs?
-----------

Graphs are a natural choice for covalently bonded molecules:

* Nodes map to atoms; edges map to bonds.
* NetworkX provides O(1) neighbour lookup, fast iteration, and
  out-of-the-box drawing utilities.
* The graph structure stores all chemical information as node / edge
  attributes — no additional data structures are needed.
* Both topological metrics (linkage counts, branching) and SMILES strings
  can be derived directly from the graph.

Atomistic graph (``G``)
-----------------------

``polymer.G`` is a :class:`networkx.Graph` where:

* Every node is a **heavy atom** (C or O).
* Every edge is a **covalent bond**.

Node attributes
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 12 60

   * - Attribute
     - Type
     - Description
   * - ``element``
     - str
     - Atomic symbol: ``"C"`` or ``"O"``
   * - ``aromatic``
     - bool
     - ``True`` if the atom belongs to the phenyl ring (indices 1–6)
   * - ``group``
     - str | None
     - Functional group: ``None``, ``"OCH3"``, ``"4OH"``, or ``"9OH"``
   * - ``index``
     - int
     - 1-based chemical carbon numbering (see :doc:`lignin_chemistry`)
   * - ``mtype``
     - str
     - Monomer type: ``"H"``, ``"G"``, or ``"S"``
   * - ``bonding``
     - bool
     - ``True`` if the atom is still available for a new inter-monomer bond
   * - ``mi``
     - int
     - Monomer index within the polymer chain (0-based)
   * - ``color``
     - str
     - CSS colour for visualisation (set per monomer type)

Edge attributes
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 12 60

   * - Attribute
     - Type
     - Description
   * - ``order``
     - int
     - Bond order: ``1`` = single, ``2`` = double
   * - ``btype``
     - str | None
     - Inter-monomer linkage name (e.g. ``"beta-O-4"``); ``None`` for
       intra-monomer bonds
   * - ``index``
     - tuple
     - (C1 index, C2 index) of the bonding atom pair
   * - ``mtype``
     - tuple
     - (monomer_type_atom1, monomer_type_atom2)

Atom counts per monomer type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

=====  ============  ===========================================================
Type   Heavy atoms   Notes
=====  ============  ===========================================================
H      11            9 C + 2 O (no OCH₃)
G      13            9 C + 2 O + 1 OCH₃ = 11 C + 2 O... wait: 9C + 2O + 1O + 1C = 11C + 3O → 13 total
S      15            9 C + 2 O + 2 OCH₃ = 13 C + 4 O — 2 shared O... = 15 total
=====  ============  ===========================================================

These counts are used by ``CharacterizeGraph.count_types()`` to infer the
monomer composition from the raw atom counts.

Coarse-grained graph (``bigG``)
--------------------------------

``polymer.bigG`` is derived automatically alongside ``G``.  It collapses
every monomer into a single node and every inter-monomer linkage into a
single edge.

Node attributes
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 12 60

   * - Attribute
     - Type
     - Description
   * - ``mtype``
     - str
     - Monomer type (``"H"``, ``"G"``, or ``"S"``)
   * - ``color``
     - str
     - Visualisation colour

Edge attributes
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 12 60

   * - Attribute
     - Type
     - Description
   * - ``btype``
     - str
     - Inter-monomer linkage name (e.g. ``"5-5"``)

The ``bigG`` size is ``len(polymer.bigG) = monomer_count``, while
``len(polymer.G) ≈ monomer_count × (11–15)`` heavy atoms.

3-letter residue codes
-----------------------

All exported topologies (JSON, PDB, SMILES statistics) use 3-letter
residue codes following the PDB/mmCIF convention:

.. list-table::
   :header-rows: 1
   :widths: 15 50 25

   * - Code
     - Full name
     - Internal type
   * - ``HPU``
     - p-Hydroxyphenyl Unit
     - H
   * - ``GYU``
     - Guaiacyl Unit
     - G
   * - ``SYU``
     - Syringyl Unit
     - S

The mapping is defined in :data:`lignoforge.core.rules.MONOMER_RESIDUE_CODE`
and applied in :class:`~lignoforge.structure.generator.MolecularStructureGenerator`
when building the topology dicts.

Multi-scale hierarchy
----------------------

::

    Population  (list[Polymer])
        │
        ├── Polymer[0]
        │       ├── G        – atomistic graph (~200–600 nodes)
        │       └── bigG     – CG graph         (~15–50 nodes)
        │
        ├── Polymer[1]
        │       ├── G
        │       └── bigG
        │
        └── ...

The hierarchy is reflected in the JSON topology outputs:

* ``atomistic_topology.json`` — one record per chain, listing monomers and
  atoms with 3-letter codes and (optionally) 3-D coordinates.
* ``coarse_grained_topology.json`` — one record per chain, listing beads
  (monomers) and links, with coarse-grained spring-layout coordinates.

Linkage encoding
-----------------

A linkage between two monomers is recorded as:

1. One or more **edges** in ``G`` with ``btype`` set to the linkage name
   (e.g. ``"beta-O-4"``).  Some linkages span multiple edges (ring motifs).
2. One **edge** in ``bigG`` with ``btype`` set to the same linkage name.
3. The bonding atoms have ``bonding = False`` after the linkage is formed,
   marking them as unavailable for future bonds.

The ``btype`` attribute on edges in ``G`` is set during
:meth:`~lignoforge.core.polymer.PolymerGraph.connect_C1_C2` by looking
up the ``(C1, C2)`` index pair in
:data:`~lignoforge.core.rules.linkage_index_to_name`.

Visualisation
--------------

Three convenience drawing functions are available in
:mod:`lignoforge.core.utils`:

.. code-block:: python

    from lignoforge.core.utils import (
        draw_big_graph,
        draw_atomic_graph,
        draw_graph,
    )
    import matplotlib.pyplot as plt

    draw_big_graph(polymer.bigG)    # hexagonal nodes, monomer types as labels
    plt.show()

    draw_atomic_graph(polymer.G)    # circular nodes, element symbols as labels
    plt.show()

The interactive HTML coarse-grained viewer (``coarse_grained_topology_viewer.html``)
provides a more polished rendering; see :doc:`../user_guide/output_formats`.
