================
Lignin Chemistry
================

This page summarises the chemical background needed to understand how
LignoForge constructs polymer models.

Lignin in context
-----------------

Lignin is a complex aromatic heteropolymer that constitutes 15–35 % by
mass of lignocellulosic biomass.  Together with cellulose and hemicellulose
it forms the structural scaffold of plant cell walls.  Its high carbon
content and aromatic backbone make it a valuable feedstock for renewable
aromatic chemicals and carbon fibres, but its structural heterogeneity
complicates characterisation and valorisation.

Lignin is biosynthesised by radical end-wise coupling of three
phenylpropanoid monolignols:

* **p-coumaryl alcohol** → H (p-hydroxyphenyl) unit
* **coniferyl alcohol** → G (guaiacyl) unit
* **sinapyl alcohol** → S (syringyl) unit

The relative abundance of H, G, and S monomers depends on the plant
species and tissue type:

.. list-table::
   :header-rows: 1
   :widths: 25 20 20 20

   * - Biomass type
     - H fraction
     - G fraction
     - S fraction
   * - Softwood (e.g. pine)
     - 0–3 %
     - 90–99 %
     - 0–3 %
   * - Hardwood (e.g. poplar)
     - 0–5 %
     - 40–55 %
     - 45–58 %
   * - Grass / agricultural residue
     - 10–30 %
     - 45–65 %
     - 15–35 %

Monomer structure
-----------------

Every monomer shares the same **phenylpropane carbon skeleton** (C9 unit):
an aromatic ring (C1–C6), an α-carbon (C7), a β-carbon (C8), and a
γ-carbon (C9).  Substituents differentiate the three types:

.. list-table::
   :header-rows: 1
   :widths: 10 25 20 20 15

   * - Type
     - Full name
     - C3-OCH₃
     - C5-OCH₃
     - Residue code
   * - H
     - p-Hydroxyphenyl unit
     - —
     - —
     - **HPU**
   * - G
     - Guaiacyl unit
     - ✓
     - —
     - **GYU**
   * - S
     - Syringyl unit
     - ✓
     - ✓
     - **SYU**

In LignoForge all output formats (JSON, PDB, statistics) use the 3-letter
residue codes above to follow the PDB/mmCIF convention for biomolecular
modelling.

Atom indexing
~~~~~~~~~~~~~

Within a monomer the atoms are numbered using the standard 1-based
chemical convention:

.. list-table::
   :header-rows: 1
   :widths: 20 15 60

   * - Index
     - Symbol
     - Position
   * - 1–6
     - C
     - Aromatic ring (C1–C6)
   * - 7
     - C
     - α-carbon
   * - 8
     - C
     - β-carbon
   * - 9
     - C
     - γ-carbon
   * - 10
     - O
     - C4-phenolic OH
   * - 11
     - O
     - C9-hydroxyl (γ-OH)
   * - +2 / +4
     - O, C
     - OCH3 group(s) for G and S

The integer index is stored in each graph node under the attribute
``"index"`` and drives the O(1) linkage lookup rules.

Inter-monomer linkages
-----------------------

Monomers are covalently connected through six primary **inter-monomer
linkage** types.  The linkage is identified by the pair of bonding carbon
indices (C1, C2):

.. list-table::
   :header-rows: 1
   :widths: 15 10 10 50

   * - Linkage name
     - C1 index
     - C2 index
     - Description
   * - β-O-4
     - 8
     - 4 (via O)
     - Most abundant (~45–60 % in native hardwood); β-aryl ether.  Cleaved selectively during Kraft and organosolv pulping.
   * - α-O-4
     - 4 (via O)
     - 7
     - α-aryl ether; less abundant than β-O-4.
   * - 4-O-5
     - 4 (via O)
     - 5
     - Diaryl ether; two aromatic rings bridged by oxygen.
   * - 5-5
     - 5
     - 5
     - Biphenyl; direct C–C bond between two aromatic C5.
   * - β-5
     - 8
     - 5
     - Phenylcoumaran ring system; β–C5 bond + O–α ring closure.
   * - β-β
     - 8
     - 8
     - Resinol; Cβ–Cβ bond + two O–α ring closures.
   * - β-1
     - 8
     - 1
     - Spirodienone motif; optional (disabled by default).

Effect of extraction process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The extraction (delignification) process profoundly alters the linkage
distribution:

* **Kraft** pulping cleaves most β-O-4 ether bonds (~40 % of native) and
  enriches the condensed C–C linkages (5-5, β-5, β-β).
* **Organosolv** pulping retains more β-O-4 bonds (~55 %) and produces
  lower MW lignin with narrower dispersity.
* **Soda / sulfite** pulping yield intermediate modifications.
* **DES (deep eutectic solvent)** extraction is mild and preserves the
  native linkage pattern most faithfully.

These trends are encoded in the literature reference tables in
:mod:`lignoforge.priors.literature_data` and are automatically applied
by :class:`~lignoforge.priors.estimator.LigninPriorEstimator`.

Characterisation observables
-----------------------------

LignoForge computes and optimises the following experimentally
accessible metrics:

.. list-table::
   :header-rows: 1
   :widths: 30 60

   * - Observable
     - Definition
   * - S / G / H fractions
     - Molar percent of each monomer type
   * - Linkage fractions
     - Molar percent of each of the 7 linkage types
   * - M\ :sub:`n`
     - Number-average molecular weight
   * - M\ :sub:`w`
     - Weight-average molecular weight
   * - PDI
     - Polydispersity index M\ :sub:`w` / M\ :sub:`n`
   * - Branching coefficient
     - Fraction of monomers connected to ≥ 3 other monomers
   * - Free OH count
     - Unoccupied phenolic and γ-hydroxy groups per chain
   * - OCH₃ count
     - Methoxy groups per chain (correlated with G + 2×S content)
