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
     - Which charges were renormalised, on which atoms and types, and by how
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

Charge renormalisation report
-----------------------------

Fragment charges are OPLS-AA database values; C–O–C bridges (``opls_199`` /
``opls_179`` and ``opls_183`` / ``opls_185``) are neutral by construction.
A residue that carries other linkages can still end with a small net charge,
which is removed by shifting *every atom of that residue* by the same amount
(residual / number of atoms).  The shifts are small (largest in the test
chains: 0.013 e) and the chain is exactly neutral for any topology.

Every run documents the changes:

.. code-block:: text

   res name  atoms  raw net (e)  shift/atom (e)  max |Δ| (e)   linkages
     4 SYU      31      +0.2000        -0.00645      0.00650   beta-O-4
   ...
   By OPLS type (all residues):
   type        atoms   mean Δ (e)      min Δ      max Δ
   opls_145       10     -0.00192   -0.00650   +0.00170
   ...
   Largest single-atom change: -0.00650 e on C2 (index 86, residue 4 SYU, opls_145; -0.1150 → -0.1215)

The JSON file lists the raw charge, final charge and shift of **every atom**
and the per-type aggregation.  The same information is available in Python:

.. code-block:: python

   from lignoforge.forcefield import build_chain_topology
   topo = build_chain_topology(chain)          # chain = atomistic topology dict
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
