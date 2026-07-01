=========
Changelog
=========

All notable changes to LignoForge are documented here.

---------
v0.2.1
---------

*2026-07 patch.*

Patch updates
~~~~~~~~~~~~~

* **Improved OPLS-AA sp³ carbon types in** ``lignin_ff/lignin.rtp`` **and**
  ``tools/assign_chain_types.py``: Cα (sp³ secondary alcohol) now uses
  ``opls_219`` (benzyl-alcohol type, q = +0.260 e) instead of the generic
  ``opls_157``; Cβ / Cα-ether (sp³ C bearing one ether bond) now uses
  ``opls_183`` (isopropyl-ether type, q = +0.170 e).  C1 (Cipso) is
  promoted from ``opls_145`` to the more specific ``opls_221`` (substituted
  aryl C).  All LJ parameters are unchanged (σ = 3.50 Å for CT, σ = 3.55 Å
  for CA); only partial charges differ.  Per-residue charge neutrality is
  maintained in all nine residue types via ``balance_charges()``, which
  adjusts q(C1) to absorb the sp³ imbalance (GYU/HPU/SYU: C1 = −0.085 e;
  HNM/GNM/SNM: C1 = −0.055 e; GYM/HPM/SYM: C1 = 0.000 e).  Validated by
  ``pdb2gmx`` (total charge = 0.000 e) and energy minimisation convergence.

---------
v0.2.0
---------

*First stable release by the DAIMON Team (2026).*

New features
~~~~~~~~~~~~

* **Modular package layout** — ``core``, ``simulation``, ``pipeline``,
  ``priors``, ``io``, and ``structure`` sub-packages provide clean,
  independently usable interfaces.

* **LigninPipeline** — high-level entry point that handles the full
  workflow (JSON input → prior estimation → simulation → export) in a
  single call.

* **LigninPriorEstimator** — three-level precedence estimator (direct
  input → constrained prior → unconditional prior) backed by structured
  literature tables.

* **LigninExporter** — unified exporter class writing JSON topologies,
  SMILES, SDF, PDB, and an interactive HTML viewer.

* **Three-letter PDB residue codes** — ``GYU`` (guaiacyl), ``SYU``
  (syringyl), ``HPU`` (p-hydroxyphenyl) are used in all structure
  outputs for compatibility with standard molecular-dynamics toolchains.

* **Hybrid 3-D embedding strategy** — heavy atoms embedded first, then
  hydrogens added with coordinate inference; MMFF94 → UFF fallback.

* **Ring-closure control** — ``i_max_ring`` and
  ``branching_propensity`` exposed as first-class parameters.

* **``LigninResults`` dataclass** — structured return value carrying
  priors, simulation kwargs, polymer list, and artefact paths.

* **JSON Schema input validation** — all inputs validated against
  the bundled ``lignin_info_schema.json`` at load time.

Patch updates
~~~~~~~~~~~~~

* **``lignoforge-chain`` CLI tool** — standalone command registered as a
  ``console_scripts`` entry point (``lignoforge.cli.build_chain:main``).
  Accepts a JSON input file plus per-run overrides for chain size
  (``--n-monomers`` / ``--mw-target``), monomer composition
  (``--monomer-type``, ``--G-fraction``, …), all seven linkage fractions,
  branching propensity, random seed, number of chains, and output formats.
  Each run produces per-chain files, a ``_stats.json`` summary, and a
  ``_manifest.json`` for full reproducibility.
  See :doc:`/user_guide/cli_tools` for the full option reference.

* **Bug fix — aromatic ring perception**: ``_graph_to_ordered_rdkit_mol``
  in ``lignoforge.core.utils`` now correctly passes ``BondType.AROMATIC``
  for intra-monomer ring bonds and ``SetIsAromatic(True)`` for aromatic
  atoms, eliminating the cyclohexane representation that caused wrong
  hydrogen counts in SMILES and 3-D structures.

* **Bug fix — α-OH on β-O-4 Cα**: ``connect_C1_C2`` in
  ``lignoforge.core.polymer`` now adds the α-hydroxy group on Cα when
  forming a β-O-4 linkage.  This corrects the molecular formula of all
  β-O-4-containing structures (validated against PubChem 517057 for the
  G-G dimer: C\ :sub:`20`\ H\ :sub:`24`\ O\ :sub:`7`).

