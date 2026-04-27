.. _force_field:

=======================================
Force-Field Setup (OPLS-AA / GROMACS)
=======================================

LignoForge ships OPLS-AA force-field parameter files in ``lignin_ff/`` that
enable all-atom molecular dynamics simulations of lignin polymer chains with
**GROMACS**.  This page describes the residue definitions, explains how to
assign per-atom types to any LignoForge PDB, and shows the full workflow from
structure generation to GROMACS topology.

.. contents:: On this page
   :depth: 2
   :local:

----

Overview
--------

The ``lignin_ff/`` directory contains:

.. code-block:: text

   lignin_ff/
   ├── lignin.rtp                   GROMACS residue topology (9 residues)
   ├── residuetypes_lignin.dat      pdb2gmx residuetype entries
   ├── notes/
   │   └── atom_type_assignment.md  Per-atom type justification
   └── tools/
       ├── assign_chain_types.py    Automatic type-assignment script
       └── test_assign.py           Validation tests

All nine residues carry **net charge = 0.000 e** (verified per OPLS-AA
charge group; see :ref:`charge-groups`).

----

Residue Definitions
-------------------

Three categories of residues are defined in ``lignin.rtp``:

β-O-4 internal chain units
~~~~~~~~~~~~~~~~~~~~~~~~~~~

These represent interior monomers in a β-O-4 linked polymer chain.  The
side chain is saturated (sp³).  Each unit has **two inter-residue bonds** in
``[ bonds ]`` that encode the β-O-4 ether linkage:

.. code-block:: text

   O4H  -CB     ; aryl ether O (this unit) → Cβ of previous residue
   CB   +O4H    ; Cβ (this unit) → aryl ether O of next residue

.. list-table::
   :header-rows: 1
   :widths: 12 20 68

   * - Name
     - Monomer type
     - Description
   * - ``GYU``
     - G (guaiacyl)
     - –OCH₃ at C3; sp³ α-OH at Cα; β-O-4 ether at Cβ
   * - ``HPU``
     - H (*p*-hydroxyphenyl)
     - No OMe substituents; sp³ α-OH at Cα; β-O-4 ether at Cβ
   * - ``SYU``
     - S (syringyl)
     - –OCH₃ at C3 **and** C5; sp³ α-OH at Cα; β-O-4 ether at Cβ

Isolated free monolignols
~~~~~~~~~~~~~~~~~~~~~~~~~~

Free monolignol model compounds with a vinyl (sp² E-configured) side chain
and a free phenol –OH at C4.  No inter-residue bonds.

.. list-table::
   :header-rows: 1
   :widths: 12 35 20 15

   * - Name
     - Compound
     - Formula
     - MW (g/mol)
   * - ``GYM``
     - coniferyl alcohol
     - C₁₀H₁₂O₃
     - 180.20
   * - ``HPM``
     - *p*-coumaryl alcohol
     - C₉H₁₀O₂
     - 150.17
   * - ``SYM``
     - sinapyl alcohol
     - C₁₁H₁₄O₄
     - 210.23

Neutral saturated references
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Dihydro- analogues used during parametrisation verification and for MD of
non-polymerised units.

.. list-table::
   :header-rows: 1
   :widths: 12 45

   * - Name
     - Compound
   * - ``GNM``
     - dihydroconiferyl alcohol (G)
   * - ``HNM``
     - dihydro-*p*-coumaryl alcohol (H)
   * - ``SNM``
     - dihydrosinapyl alcohol (S)

----

Atom-Name Conventions
---------------------

LignoForge generates PDB files that use the atom names expected by
``lignin.rtp``.  The table below lists all expected names:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Position
     - Atom names
   * - Aromatic ring
     - ``C1 C2 H2 C3 C4 C5 C6`` (plus heteroatoms at substituted positions)
   * - Free phenol (monolignols only)
     - ``O4H HO4``
   * - OMe at C3 (G, S types)
     - ``OM3 CM3 HM31 HM32 HM33``
   * - OMe at C5 (S type only)
     - ``OM5 CM5 HM51 HM52 HM53``
   * - Vinyl side chain (monolignols)
     - ``CA HA CB HB``
   * - sp³ side chain (chain units)
     - ``CA HA OA HOA CB HB``
   * - γ-Alcohol
     - ``CG HG1 HG2 OG HOG``

----

Automatic Type Assignment
--------------------------

A raw LignoForge PDB contains residue names but no atom-type information.
The script ``lignin_ff/tools/assign_chain_types.py`` reads the PDB, detects
all inter-residue bonds, identifies each linkage type from the bonded atom
names, and produces:

- A **custom per-chain RTP** with correct OPLS-AA types at all linkage sites.
- An optional **renamed PDB** with ≤3-character GROMACS residue names.

The script supports all seven standard lignin linkage types:
β-O-4, 5-5, 4-O-5, β-5, β-β, α-O-4, and β-1.

.. code-block:: bash

   # Minimal usage — outputs chain_custom.rtp and chain_renamed.pdb
   python lignin_ff/tools/assign_chain_types.py chain.pdb

   # Explicit output paths
   python lignin_ff/tools/assign_chain_types.py chain.pdb \
       -o chain_custom.rtp \
       --renamed-pdb chain_renamed.pdb

Using the API directly:

.. code-block:: python

   import sys
   sys.path.insert(0, "lignin_ff/tools")
   from assign_chain_types import assign_chain_types

   typed = assign_chain_types(
       "chain.pdb",
       rtp_out="chain_custom.rtp",
       renamed_pdb_out="chain_renamed.pdb",
   )

   # typed[residue_seq]["rtp_name"]  → e.g. "GYU", "HPM", …
   # typed[residue_seq]["atoms"]     → {atom_name: (opls_type, charge, cgnr), …}
   # typed[residue_seq]["linkages"]  → [(linkage_name, position), …]

----

Full Workflow: LignoForge → GROMACS
-------------------------------------

**Step 1 — Generate an atomistic PDB with LignoForge**

.. code-block:: python

   from lignoforge.core.monomer import Monomer
   from lignoforge.core.polymer import Polymer
   from lignoforge.structure.pdb import PDBStructureWriter

   m = Monomer("G", monomer_index=0)
   m.create()
   p = Polymer(m, verbose=False)
   p.add_specific_monomer("G", "beta-O-4")
   p.add_specific_monomer("S", "beta-O-4")

   PDBStructureWriter().write_polymer_pdb(p, "chain.pdb", optimize_3d=True)

Or via the pipeline (see :doc:`pipeline_usage`):

.. code-block:: python

   from lignoforge.pipeline import LigninPipeline

   results = LigninPipeline.from_json("examples/lignin_input_example.json").run()
   # PDB files are written to results.output_dir/pdb/

**Step 2 — Assign OPLS-AA types**

.. code-block:: bash

   python lignin_ff/tools/assign_chain_types.py chain.pdb
   # → chain_custom.rtp
   # → chain_renamed.pdb

**Step 3 — Install the force-field files**

For a project-local force field (recommended for reproducibility):

.. code-block:: bash

   cp -r /usr/share/gromacs/top/oplsaa.ff ./oplsaa.ff
   cp lignin_ff/lignin.rtp oplsaa.ff/
   cat lignin_ff/residuetypes_lignin.dat >> oplsaa.ff/residuetypes.dat

**Step 4 — Run pdb2gmx**

.. code-block:: bash

   gmx pdb2gmx \
       -f chain_renamed.pdb \
       -ff ./oplsaa \
       -water tip3p \
       -o processed.gro \
       -p topol.top \
       -rtpres chain_custom.rtp

``pdb2gmx`` will look up each residue name in ``chain_custom.rtp``, assign
bonds, angles, dihedrals, and impropers automatically, and produce a complete
GROMACS topology.

----

Running the Validation Tests
------------------------------

.. code-block:: bash

   python lignin_ff/tools/test_assign.py
   # or
   python -m pytest lignin_ff/tools/test_assign.py -v

The test suite (13 tests) checks:

- Isolated H, G, S monomers → correct residue name and free-phenol types.
- GGG, HH, SSS β-O-4 chains → correct aryl-ether and sp³-ether types.
- G–G 5-5, G–H 4-O-5, G–G β-1 dimers → linkage-specific type modifications.
- Net charge = 0.000 e for **every** residue in every test.

See :ref:`atom-type-table` for the scientific background behind the type
assignments.
