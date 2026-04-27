=========
Structure
=========

The ``lignoforge.structure`` sub-package converts the internal graph
representation into molecular structure files (PDB, SDF, 3-D mol).

.. contents:: Contents
   :local:
   :depth: 2

MolecularStructureGenerator
----------------------------

.. autoclass:: lignoforge.structure.pdb.MolecularStructureGenerator
   :members:
   :undoc-members:
   :show-inheritance:

This class wraps
:func:`~lignoforge.core.utils.graph_to_mol` for pipeline-level use.  It
adds chain-ID assignment and optional 3-D coordinate export::

    from lignoforge.structure.generator import MolecularStructureGenerator

    gen = MolecularStructureGenerator(
        polymer,
        generate_3d=True,
        include_hydrogens=True,
        optimize_3d=True,
    )
    mol = gen.to_mol()
    gen.write_sdf("chain_A.sdf")

3-D embedding strategy
~~~~~~~~~~~~~~~~~~~~~~~

1. The heavy-atom graph is embedded with RDKit ETKDG (``EmbedMolecule``).
2. Hydrogens are added *with coordinate inference* (``AddHs(addCoords=True)``).
3. If ``optimize_3d=True``, MMFF94 geometry optimisation is attempted
   first; UFF is used as a fallback for monomers that are outside the
   MMFF94 parameter space.

PDBStructureWriter
------------------

.. autoclass:: lignoforge.structure.pdb.PDBStructureWriter
   :members:
   :undoc-members:
   :show-inheritance:

Writes ``HETATM``-formatted PDB files with:

* Residue names from :data:`~lignoforge.core.rules.MONOMER_RESIDUE_CODE`
  (``GYU`` / ``SYU`` / ``HPU``)
* Sequential atom serial numbers across the full chain
* One file per chain (``chain_A.pdb``, ``chain_B.pdb``, …)

Usage::

    from lignoforge.structure.generator import PDBStructureWriter

    writer = PDBStructureWriter(polymer, output_dir="pdb/")
    writer.write("chain_A")

.. seealso::

   :doc:`/user_guide/output_formats` — full description of the PDB
   output format and residue naming conventions.
