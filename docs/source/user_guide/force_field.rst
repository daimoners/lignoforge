.. _force_field:

=======================================
Force-Field Setup (OPLS-AA / GROMACS)
=======================================

LignoForge builds an OPLS-AA topology for any generated chain directly from
its bond graph (package ``lignoforge.forcefield``).  There is no ``pdb2gmx``
step and no residue-template file: atom types, charges, angles, dihedrals,
1-4 pairs and impropers are all derived from the graph, so linear, branched
and ring-closed (β-5, β-β) chains are handled identically.

.. contents:: On this page
   :depth: 2
   :local:

----

Quick start
-----------

.. code-block:: bash

   lignoforge-chain --n-monomers 10 --format gromacs,pdb
   cd chain_output/gromacs/chain_0
   gmx grompp -f em.mdp -c chain_0.gro -p chain_0.top -o em.tpr
   gmx mdrun -deffnm em

Files written per chain:

.. list-table::
   :widths: 35 65

   * - ``<name>.top``
     - Self-contained topology; includes ``oplsaa.ff/forcefield.itp``
       (GROMACS must find ``oplsaa.ff``: system install or ``GMXLIB``).
       Bonded parameters are left blank and resolved by ``grompp`` from
       ``ffbonded.itp``, so a missing parameter is an explicit error.
   * - ``<name>.gro``
     - Coordinates (nm) centred in a cubic box.
   * - ``em.mdp``
     - Energy-minimisation parameters.
   * - ``<name>_topology_report.json``
     - Linkages and per-residue charge summary.
   * - ``<name>_charge_report.txt / .json``
     - Junction charge corrections, the atoms and OPLS types changed, and by how
       much (see :ref:`charge-report`).

For the complete simulation workflow (solvation, equilibration, production,
coarse-grained models) see :ref:`md_workflows`.

----

How types are assigned
----------------------

Typing uses only the atom names of each residue and the *linkage context*
(which atoms take part in an inter-monomer bond, and of which linkage).
Every inter-monomer bond of the graph must match one of the valid atom pairs
below; anything else raises ``TypingError`` instead of being guessed.

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Linkage
     - Bonded atom pair(s)
     - Notes
   * - β-O-4
     - ``CB`` – ``O4H``
     - Cα receives an α-OH
   * - α-O-4
     - ``CA`` – ``O4H``
     - Cβ left as CH\ :sub:`2`
   * - 4-O-5
     - ``C5`` – ``O4H``
     -
   * - 5-5
     - ``C5`` – ``C5``
     -
   * - β-5
     - ``CB`` – ``C5`` and ``CA`` – ``O4H``
     - second bond closes the coumaran ring
   * - β-β
     - ``CB`` – ``CB`` and ``CA`` – ``OG`` (×2)
     - resinol: two furan-ring closures

Atom names follow the generator: ring ``C1 C2 H2 C3 C4 C5 C6``; phenol
``O4H HO4``; methoxy ``OM3 CM3 HM31-33`` and ``OM5 CM5 HM51-53``; side chain
``CA HA OA HOA CB HB CG HG1 HG2 OG HOG``.  The full type table with
justifications is in :ref:`opls_aa_parametrization`.

----

.. _charge-report:

Charge report
-------------

Fragment charges are OPLS-AA database values.  Residues are **not** forced to
be neutral: where a bond joins two residues the neutral OPLS group is split
between them, so a β-O-4 ether leaves the two residues at −0.2 e and +0.2 e.
These charges are physical.  What the OPLS fragments do get wrong is small and
specific to each junction type: every inter-monomer bond (and every attached
derivative group) leaves a fixed *excess* of charge, additive over the chain
(β-O-4 and 5-5: 0; α-O-4: −0.025 e; 4-O-5: −0.115 e; β-5: −0.0125 e per bond;
β-β: +0.103 e per bond; derivatives: γ-esters +0.06 e, benzylic −0.025 e).
LignoForge removes each excess **locally**, spreading it uniformly over the
atoms of the junction and their directly bonded neighbours, so the chain is
exactly neutral for any topology and no atom far from a junction is touched
(largest single-atom change in random test chains: 0.046 e; most chains change
by less than 0.02 e).

Why not neutralise each residue?  Because that moves the ±0.2 e of a split
ether group onto the whole residue, adding a spurious 4–6 D dipole across a
β-O-4 dimer.  Against HF/6-31G\* electrostatic potentials of G–G dimers the
per-residue scheme gave RMS errors of 5–7 kcal/mol/e for β-O-4, α-O-4, β-5,
4-O-5 and β-β; the junction scheme gives 2.7–3.4 kcal/mol/e, the same as the
raw OPLS charges where those are already neutral (β-O-4, 5-5).  The approach
follows the modular-linkage-charge practice of GLYCAM and LignAmb25, where
linking atoms carry fixed charges instead of residues being neutralised.

Every run documents the changes:

.. code-block:: text

   Junction corrections:
   junction               atoms (core)             excess (e) atoms hit  shift/atom (e)
   link: 4-O-5            C5(1) – O4H(3)              -0.1150         5        +0.02300
   link: 4-O-5            O4H(1) – C5(2)              -0.1150         5        +0.02300

    res name   raw net (e)  final net (e)  max |Δ| (e)   linkages
      1 GYU        -0.1150        +0.0000      0.04600   4-O-5
      2 GYU        +0.0850        +0.1540      0.02300   4-O-5
      3 GYU        -0.2000        -0.1540      0.02300   4-O-5
   ...
   Largest single-atom change: +0.04600 e on C4 (index 4, residue 1 GYU, opls_199; +0.0850 → +0.1310)

The JSON file lists the junction corrections, the raw charge, final charge and
shift of **every atom** and the per-type aggregation.  The same information is available in Python:

.. code-block:: python

   from lignoforge.forcefield import build_chain_topology
   topo = build_chain_topology(chain)          # chain = atomistic topology dict
   topo.charge_corrections                     # per junction
   topo.charge_adjustments()                   # per atom
   topo.charge_adjustments_by_type()           # per OPLS type
   print(topo.charge_report_text())

----

Torsions and impropers
----------------------

* All proper dihedrals (Ryckaert-Bellemans, func 3) and 1-4 pairs are
  generated for every central bond.
* Ring carbons use ``improper_Z_CA_X_Y``; sp\ :sup:`2` vinyl carbons
  (``opls_142``) use ``improper_Z_CM_X_Y``.
* Two torsions are absent from OPLS-AA and are added to the ``.top`` as
  analogues: ``CA CA CM HC`` (zero, as for ethylbenzene ``CA CA CT HC``) and
  ``CA CT CT CA`` (hydrocarbon values, as ``CA CT CT CT``).

----

Validating a setup
------------------

``tools/validate_pipeline.py`` runs growth → chemistry → clash check →
topology → ``grompp`` → EM → short MD → CG for a set of reference and random
chains::

   python tools/validate_pipeline.py --md --random 10
   python -m pytest tests -v
