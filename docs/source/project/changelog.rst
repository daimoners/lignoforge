=========
Changelog
=========

All notable changes to LignoForge are documented here.

---------
v0.2.2
---------

*Unreleased.*

New features
~~~~~~~~~~~~

* **Graph-based GROMACS topology generation** (``lignoforge.forcefield``,
  ``lignoforge-chain --format gromacs``).  Writes ``.top``, ``.gro`` and
  ``em.mdp`` directly from the bond graph, with OPLS-AA types assigned from
  atom names and linkage context.  No ``pdb2gmx`` and no assumption about
  residue order, so linear, branched and ring-closed chains are all valid.
  All angles, proper dihedrals, 1-4 pairs and aromatic impropers are derived
  from the graph; bonded parameters are resolved by ``grompp`` from
  ``oplsaa.ff``.
* **Charge scheme**: OPLS-AA fragment charges with a **local correction per junction**
  instead of per-residue neutralisation.  Residues are not forced to be neutral (a
  β-O-4 ether leaves them at ∓0.2 e, which is physical); each junction type (linkage
  bond or attached derivative) leaves a fixed excess that is removed on the junction
  atoms and their bonded neighbours, so the chain total is exactly 0 for any topology.
  Against HF/6-31G* ESPs of G–G dimers this lowers the RMS error from 5–7 to 2.7–3.4
  kcal/mol/e for β-O-4, α-O-4, β-5, 4-O-5 and β-β.  The charge report lists the
  junction corrections (``topo.charge_corrections``).  *Behaviour change:* charges of
  existing chains differ slightly from earlier versions.

* **Atomistic MD workflow** (``lignoforge.md``, ``--format md``): box,
  solvation (TIP3P/SPC/E or vacuum), EM, NVT, NPT and production ``.mdp``
  files plus an executable ``run.sh``.
* **Coarse-grained model** (``lignoforge.cg``, ``--format cg``): one bead per
  monomer, linkage-keyed harmonic bonds and angles, LJ beads, stochastic
  dynamics, multi-chain boxes via ``gmx insert-molecules``.  Parameters are
  serialisable (``CGParameters``) and can be derived from atomistic MD with
  Boltzmann inversion (``lignoforge-cg-fit``).  Shipped defaults are
  provisional.

Removed
~~~~~~~

* ``lignin_ff/`` (``lignin.rtp``, ``residuetypes_lignin.dat``,
  ``tools/assign_chain_types.py``, ``tools/test_assign.py`` and notes): the
  ``pdb2gmx`` / RTP route linked only residues *i* and *i+1* and is replaced
  by ``lignoforge.forcefield``.  The files remain in the git history (v0.2.1).

* **Virtual chemical derivatisation** (``lignoforge.core.derivatization``,
  ``lignoforge-chain --derivatize SITE:REAGENT[:DEGREE]``, ``derivatization`` in
  the chain specification of the Python/REST API): methylation, acetylation and
  longer acyl esterification of phenolic, benzylic (α) and γ hydroxyls, with a
  controllable degree of substitution.  The modified molecules use OPLS-AA
  ester/ether types (opls_465/466/467, aryl ester 472/473, 490/491, 180) with
  planarity impropers on the carbonyl carbons; four ester torsions absent from
  ``ffbonded.itp`` are added as analogues.  Carboxylic acids and carbonyls are
  not part of the monomer model, so their reactions are not available.
  The web interface exposes it in the *Build* page (card "Chemical modification") and
  reports the achieved degree of substitution on the chain page.
* **CG parametrisation tool** (``lignoforge.cgfit``, ``lignoforge-cg-param``,
  *CG parametrisation* page of the web interface): builds and runs an atomistic
  OPLS-AA bulk reference of lignin with GROMACS, maps it onto beads (two beads
  per monomer or one), fits bonded terms by Boltzmann inversion and Lennard-Jones
  non-bonded terms by relative-entropy minimisation with a pressure correction,
  validates the result at constant pressure and exports a CG bulk system.  Includes
  an optimiser self-test on a synthetic reference.  See :ref:`cg_parametrisation`.
* **Project API** (``lignoforge.api``): UI-independent service layer
  (``Project``) with on-disk, versioned records for chains, local GROMACS runs,
  analyses and CG parameter sets; shared by CLI, notebooks and the web
  interface.  See :ref:`project_api`.
* **Trajectory analysis** (``lignoforge.analysis``, MDAnalysis): radius of
  gyration, end-to-end distance, RMSD, density, RDF, inter-monomer contact
  map and energy terms, for atomistic and CG runs.
* **Web interface** (React + TypeScript, ``frontend/``; bundle shipped in
  ``lignoforge/web/static``): project management, chain builder with live
  specification preview, 3-D viewer (3Dmol.js; colour by monomer / element /
  OPLS type / charge / residue, atomistic and CG views), force-field and
  charge-report browsers, linkage graph, local simulation launcher with live
  log and progress, analysis charts (time series, RDF, contact map) with
  export, light/dark themes, environment and model-validation status page.
  Covered by Playwright end-to-end tests against the real backend.
* **Web backend** (``lignoforge.web``, FastAPI) and the ``lignoforge``
  command (``chain``, ``cg-fit``, ``gui``, ``check``).  Optional extras:
  ``pip install 'lignoforge[gui]'`` and ``'lignoforge[analysis]'``.

Bug fixes
~~~~~~~~~

* MD runs are now reproducible: the velocity seed is explicit and recorded.
* ``nstenergy`` / ``nstlog`` are always multiples of ``nstcalcenergy`` (GROMACS
  refused runs with output intervals such as 0.5 ps).
* ``--beta-1`` was accepted but the bond is never built; its fraction is now
  set to 0 with a warning, so the reported linkage distribution matches the
  chain.
* sp\ :sup:`2` vinyl carbons now get ``improper_Z_CM_X_Y`` planarity terms.
* Charge renormalisation is reported per atom and per OPLS type
  (``*_charge_report.txt/.json``).
* ``--branching 0`` (and the default) now yields strictly linear chains.
  Previously 0 was treated as "unrestricted" in the exact-size growth path,
  so any monomer with a free C4/C5/Cβ site could branch.
* Cα that already carries an α-OH (β-O-4 acceptor) can no longer accept a
  second oxygen (α-O-4 / β-5 / β-β), which produced chemically invalid acetals.
* Chain growth that cannot reach the requested size now emits a
  ``RuntimeWarning`` instead of silently returning a shorter chain.
* The legacy RTP lacked ``all_dihedrals`` and generated only 57 of 332
  dihedrals for a 5-mer.  Added missing OPLS torsions ``CA-CA-CM-HC`` and
  ``CA-CT-CT-CA``.

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

