====================
Manual Construction
====================

Advanced users can bypass the pipeline entirely and construct polymer
chains atom by atom using the low-level :class:`~lignoforge.core.monomer.Monomer`
and :class:`~lignoforge.core.polymer.Polymer` API.

Creating individual monomers
------------------------------

::

    from lignoforge.core import Monomer

    g = Monomer("G")
    g.create()                          # builds the atomistic graph

    print(g.G.number_of_nodes())        # 13 heavy atoms
    print(g.bigG.number_of_nodes())     # 1 coarse-grained bead

Low-level factory functions are also available::

    from lignoforge.core.monomer import monomer_H, monomer_G, monomer_S

    g_graph = monomer_G()               # returns nx.Graph directly

Initialising a polymer
-----------------------

::

    from lignoforge.core import Monomer, Polymer

    seed = Monomer("G")
    seed.create()
    polymer = Polymer(seed)

A ``Polymer`` inherits all attributes from its seed monomer and adds
step-wise growth methods.

Adding monomers
----------------

Specific monomer and specific linkage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    polymer.add_specific_linkage(linkage_type="beta-O-4", monomer_type="S")

Both arguments are required.  Raises ``ValueError`` if the combination is
chemically invalid.

Specific linkage, random compatible monomer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    # Adds a random H-or-G monomer via a 5-5 bond
    polymer.add_specific_linkage(linkage_type="5-5")

Specific monomer, random compatible linkage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    polymer.add_specific_monomer(monomer_type="G")

Random monomer and random linkage
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    polymer.add_random_monomer()

Adding ring closures
----------------------

::

    # Specific ring type
    polymer.add_specific_ring(linkage_type="beta-5")

    # Random ring (beta-5 or beta-beta)
    polymer.add_random_ring()

Ring-closure methods require at least two available bonding atoms that
satisfy the ring geometry; they return ``False`` silently if no valid
ring is possible.

Branching
----------

To attach a new monomer to an **interior** node (branching) rather than a
terminal node, pass ``branching_state=True`` to ``add_specific_monomer``::

    polymer.add_specific_monomer(monomer_type="G", branching_state=True)

Pre-defined chain recipes
--------------------------

::

    from lignoforge.core import Monomer, Polymer

    seed = Monomer("G")
    seed.create()
    p = Polymer(seed)

    # Build a short linear G-S-G-S tetramer via beta-O-4
    for mtype in ["S", "G", "S"]:
        p.add_specific_linkage(linkage_type="beta-O-4", monomer_type=mtype)

    # Close a beta-5 ring
    p.add_specific_ring(linkage_type="beta-5")

    print(len(p.bigG), "monomers")
    print(len(p.G), "heavy atoms")

Inspecting the graphs
----------------------

NetworkX operations apply directly::

    # Iterate monomer nodes in bigG
    for node, attrs in p.bigG.nodes(data=True):
        print(node, attrs["mtype"])

    # Iterate inter-monomer bonds
    for u, v, attrs in p.bigG.edges(data=True):
        print(f"  {u} — {v}  ({attrs['btype']})")

    # Count available bonding carbons
    available_C1 = [n for n, d in p.G.nodes(data=True) if d["bonding"]]
    print("Available C1 sites:", len(available_C1))

Characterising a manually built chain
---------------------------------------

::

    from lignoforge.core import Characterize

    ch = Characterize(p)
    summary = ch.summary()

    print(f"Monomers   : {summary['monomer_count']}")
    print(f"β-O-4 frac : {summary['linkage_fracs']['beta-O-4']:.2f}")
    print(f"MW         : {summary['MW']:.0f} g/mol")
    print(f"Branching  : {summary['branching_coeff']:.3f}")
    print(f"SMILES     : {summary['smiles']}")

Converting to RDKit
--------------------

::

    from lignoforge.core.utils import graph_to_mol

    mol = graph_to_mol(p.G)           # heavy atoms only, no 3D

    # With 3-D coordinates and explicit H
    mol_3d = graph_to_mol(
        p.G,
        generate_3d=True,
        add_hs_for_3d=True,
        optimize_3d=True,
        max_uff_iterations=200,
    )

    from rdkit import Chem
    Chem.MolToMolFile(mol_3d, "chain.mol")

Generating SMILES
~~~~~~~~~~~~~~~~~

::

    from lignoforge.core.utils import graph_to_smile

    smiles = graph_to_smile(p.G)
    print(smiles)
