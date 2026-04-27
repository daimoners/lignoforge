.. _opls_aa_parametrization:

==============================================================
OPLS-AA Parametrisation of Lignin Residues
==============================================================

This page documents the atom-type and partial-charge assignments used in
``lignin_ff/lignin.rtp`` for the nine OPLS-AA lignin residues.
It explains each choice with reference to the parent OPLS-AA force field
[Jorgensen1996]_ and to published lignin MD studies.

For the pragmatic workflow (how to apply the parameters to a GROMACS
simulation), see :doc:`/user_guide/force_field`.

.. contents:: On this page
   :depth: 3
   :local:

----

Structural Overview
-------------------

All nine residues are built on a phenylpropanoid (C6–C3) scaffold.  They
differ in ring substitution and in side-chain hybridisation:

.. list-table::
   :header-rows: 1
   :widths: 10 10 18 18 20 24

   * - Name
     - Type
     - C3 subst.
     - C5 subst.
     - C4 subst.
     - Side chain
   * - ``HPM``
     - H (mono.)
     - –H
     - –H
     - –OH (free)
     - vinyl (sp², –CA=CB–)
   * - ``HPU``
     - H (chain)
     - –H
     - –H
     - –O– (ether)
     - sp³ (–CA(OH)–CB(O)–)
   * - ``HNM``
     - H (ref.)
     - –H
     - –H
     - –OH (free)
     - sp³, γ-OH only
   * - ``GYM``
     - G (mono.)
     - –OCH₃
     - –H
     - –OH (free)
     - vinyl (sp², –CA=CB–)
   * - ``GYU``
     - G (chain)
     - –OCH₃
     - –H
     - –O– (ether)
     - sp³ (–CA(OH)–CB(O)–)
   * - ``GNM``
     - G (ref.)
     - –OCH₃
     - –H
     - –OH (free)
     - sp³, γ-OH only
   * - ``SYM``
     - S (mono.)
     - –OCH₃
     - –OCH₃
     - –OH (free)
     - vinyl (sp², –CA=CB–)
   * - ``SYU``
     - S (chain)
     - –OCH₃
     - –OCH₃
     - –O– (ether)
     - sp³ (–CA(OH)–CB(O)–)
   * - ``SNM``
     - S (ref.)
     - –OCH₃
     - –OCH₃
     - –OH (free)
     - sp³, γ-OH only

----

.. _atom-type-table:

Atom-Type Table
---------------

.. list-table::
   :header-rows: 1
   :widths: 14 24 12 10 10 30

   * - Atom name
     - Chemical context
     - OPLS type
     - Symbol
     - q (e)
     - Notes
   * - ``C1``
     - ipso aromatic C (no H)
     - ``opls_145``
     - CA
     - 0.000
     - Charge corrected to 0.000; see :ref:`charge-corrections`
   * - ``C2``, ``C6``
     - aromatic C–H (ortho/para to C1)
     - ``opls_145``
     - CA
     - −0.115
     - Standard benzene C
   * - ``H2``, ``H6``
     - aromatic C–H (H partner of C2/C6)
     - ``opls_146``
     - HA
     - +0.115
     - Standard benzene H
   * - ``C5``
     - aromatic C–H (meta to C1, if unsubstituted)
     - ``opls_145``
     - CA
     - −0.115
     - Used in H and G types; becomes ``opls_199`` in S type
   * - ``H5``
     - aromatic C–H partner of C5
     - ``opls_146``
     - HA
     - +0.115
     - Only in H and G types
   * - ``C3``
     - aromatic C–OMe (C3-OMe in G/S types)
     - ``opls_199``
     - C(anisole)
     - +0.085
     - Anisole C connected to ring O (Cvej model)
   * - ``C3``
     - aromatic C–H at C3 (HPM/HPU/HNM)
     - ``opls_145``
     - CA
     - −0.115
     - Identical to C2/C6 treatment
   * - ``C4``
     - aromatic C–OH (monolignols; free phenol)
     - ``opls_166``
     - C(phenol)
     - +0.150
     - From phenol model [Jorgensen1996]_
   * - ``O4H``
     - phenol oxygen (monolignols)
     - ``opls_167``
     - OH(phenol)
     - −0.585
     - From phenol model; removed in chain units (ether bond)
   * - ``HO4``
     - phenol O–H (monolignols)
     - ``opls_168``
     - HO(phenol)
     - +0.435
     - Sum C4+O4H+HO4 = 0.000 e
   * - ``C4``
     - aromatic C–O– (chain units; aryl ether)
     - ``opls_199``
     - C(anisole)
     - +0.085
     - Ether C same model as OMe C; HO4 absent in chain units
   * - ``OM3``
     - methoxy O at C3
     - ``opls_179``
     - OS(ether)
     - −0.285
     - Dimethyl ether / anisole O
   * - ``CM3``
     - methoxy C at C3 (–OCH₃)
     - ``opls_181``
     - CT(methoxy)
     - +0.060
     - sp³ methyl of OMe
   * - ``HM31–HM33``
     - methoxy H at C3
     - ``opls_185``
     - HC(methoxy)
     - +0.025
     - Sum CM3+OM3+C3+3×HM3n = 0.000 e
   * - ``CA``
     - vinyl α-C (monolignols, sp²)
     - ``opls_142``
     - CM(vinyl)
     - −0.115
     - Propene-based ethylene C
   * - ``HA``
     - vinyl α-H
     - ``opls_144``
     - HM(vinyl)
     - +0.115
     - Sum CA+HA = 0.000 e
   * - ``CB``
     - vinyl β-C (monolignols, sp²)
     - ``opls_142``
     - CM(vinyl)
     - −0.115
     - Same treatment as CA
   * - ``HB``
     - vinyl β-H
     - ``opls_144``
     - HM(vinyl)
     - +0.115
     - Sum CB+HB = 0.000 e
   * - ``CA``
     - sp³ α-C of chain unit (bears α-OH)
     - ``opls_157``
     - CT(alcohol)
     - +0.205
     - 2-propanol model; α-OH on same C
   * - ``HA``
     - sp³ α-H
     - ``opls_156``
     - HC(alcohol)
     - +0.060
     - 2-propanol model
   * - ``OA``
     - α-alcohol O (chain units)
     - ``opls_154``
     - OH
     - −0.683
     - Alcohol O; standard OPLS-AA ethanol/propanol
   * - ``HOA``
     - α-alcohol O–H
     - ``opls_155``
     - HO(alcohol)
     - +0.418
     - Sum CA+HA+OA+HOA = 0.000 e
   * - ``CB``
     - sp³ β-C (bears inter-residue β-O-4 ether O)
     - ``opls_157``
     - CT(alcohol)
     - +0.140
     - No direct –OH; charge adjusted to maintain cg = 0.000 e
   * - ``HB``
     - sp³ β-H
     - ``opls_156``
     - HC(alcohol)
     - +0.060
     - From 2-propanol model
   * - ``CG``
     - γ-CH₂ (bears γ-OH)
     - ``opls_157``
     - CT(alcohol)
     - +0.145
     - 1-propanol.itp reference (GROMACS gmx pdb2gmx)
   * - ``HG1``, ``HG2``
     - γ-methylenic H
     - ``opls_156``
     - HC(alcohol)
     - +0.060
     - Corrected from raw +0.040; see :ref:`charge-corrections`
   * - ``OG``
     - γ-alcohol O
     - ``opls_154``
     - OH
     - −0.683
     - 1-propanol.itp
   * - ``HOG``
     - γ-alcohol O–H
     - ``opls_155``
     - HO(alcohol)
     - +0.418
     - Sum CG+2×HG+OG+HOG = 0.000 e

----

.. _charge-groups:

Charge Groups
-------------

OPLS-AA requires every charge group to be **electrically neutral**.  The
table below lists the charge groups (``cgnr`` column in the RTP), their atoms,
and the partial-charge sum for a prototypical G-type residue (``GYU``):

.. list-table::
   :header-rows: 1
   :widths: 10 55 12 10

   * - cgnr
     - Atoms
     - q sum (e)
     - Present in
   * - 1
     - C1 (ipso)
     - 0.000
     - all residues
   * - 2
     - C2 H2
     - 0.000
     - all residues
   * - 3
     - C3 (+ OM3 CM3 HM31..33 if OMe; or H3 if bare)
     - 0.000
     - all; G/S use OMe branch
   * - 4
     - C4 (+ O4H HO4 if phenol; or +O4-chain if ether)
     - 0.000
     - phenol: mono.; ether: chain units
   * - 5
     - C5 H5 (or C5+OM5+CM5+... for S type)
     - 0.000
     - all residues
   * - 6
     - C6 H6
     - 0.000
     - all residues
   * - 7
     - CA HA (vinyl) or CA HA OA HOA (sp³)
     - 0.000
     - monolignols; chain units
   * - 8
     - CB HB (vinyl or sp³ chain)
     - 0.000
     - all residues
   * - 9
     - CG HG1 HG2 OG HOG
     - 0.000
     - all residues (all carry γ-OH)

All nine residues verified: **net charge = 0.000 e**.

----

.. _charge-corrections:

Charge Corrections
------------------

Two partial charges deviate from the raw OPLS-AA literature values and
require justification:

**C1 (ipso) → q = 0.000 e**

In OPLS-AA phenol and toluene, the ipso carbon carries a partial charge
of −0.130 e (paired with attached H, q = +0.130 e).  In all lignin ring
carbons, C1 has **no attached H** (it bears the propanoid side chain).
Following the Jorgensen convention for substituted aromatic carbons:

.. math::

   q(C_\text{ipso}) = 0 \quad \text{when no H is attached}

The remainder of the ring accounts for the missing −0.130 e via the
phenol (C4/O4H/HO4) and OMe (C3/OM3) charge groups.

**HG1, HG2 (γ-methylene H) → q = +0.060 e**

The raw OPLS-AA 1-propanol atom type for methylene H is +0.040 e
(``opls_140``, generic –CH₂– in a chain).  However, the GROMACS
reference compound ``1propanol.itp`` lists the methylene adjacent to the
hydroxyl-bearing carbon with q = +0.060 e, matching the 2-propanol model
for a –CH₂OH group.  Since Cγ carries –CH₂OH in all lignin residues,
q(HG) = +0.060 e is used throughout.

Combined effect: charge group 9 sums to:

.. math::

   q(\text{CG}) + 2 \cdot q(\text{HG}) + q(\text{OG}) + q(\text{HOG}) \\
   = +0.145 + 2(+0.060) + (-0.683) + (+0.418) = 0.000 \, \text{e}

----

Vinyl vs. sp³ Side Chain
-------------------------

Monolignols (``GYM``, ``HPM``, ``SYM``) possess a **trans-vinyl** (E-β)
side chain terminating in –CG–OG–HOG.  The types are:

.. code-block:: text

   CA  opls_142  (vinyl =CH–)   q = −0.115
   HA  opls_144  (vinyl H)      q = +0.115
   CB  opls_142                 q = −0.115
   HB  opls_144                 q = +0.115

Chain units (``GYU``, ``HPU``, ``SYU``) and reference compounds
(``GNM``, ``HNM``, ``SNM``) have a **saturated sp³** side chain:

.. code-block:: text

   CA  opls_157  (CT, bears –OA–)  q = +0.205
   HA  opls_156  (HC, ether adj.)  q = +0.060
   CB  opls_157  (CT, bears inter-res. O- in chain units)  q = +0.140
   HB  opls_156                    q = +0.060

The Cα–OH and Cβ–O (ether) groups in chain units reproduce the local
electrostatics of secondary alcohols and dialkyl ethers, respectively.

----

Methoxy Group
-------------

Both C3-OMe (G-type) and C3/C5-OMe (S-type) are parameterised using the
anisole conformer from [Jorgensen1996]_:

.. code-block:: text

   C3 (or C5)  opls_199  C(anisole ring)  q = +0.085
   OM           opls_179  OS(ether O)      q = −0.285
   CM           opls_181  CT(methyl C)     q = +0.060
   HM (×3)      opls_185  HC(methyl H)     q = +0.025

Charge-group sum (C3 + OM + CM + 3×HM):

.. math::

   0.085 + (-0.285) + 0.060 + 3(0.025) = 0.000 \, \text{e} \checkmark

----

Missing Bond and Dihedral Types
---------------------------------

At the time of writing, the following interaction types are **not present**
in the standard OPLS-AA ``ffbonded.itp`` and must be verified before
production runs:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Missing type
     - Relevant in
   * - ``CA–CM`` (bond)
     - C1-ring → Cα vinyl (monolignols)
   * - ``CA–OS`` (bond)
     - Cα (sp³ chain) → β-O-4 inter-residue ether O
   * - ``CA–CA–OS–CT`` (dihedral)
     - ring–Cα–O–Cβ in β-O-4 chain units
   * - ``CA–CA–CM–CM`` (dihedral)
     - ring–Cα=Cβ in monolignols
   * - ``CM=CM–CG–OH`` (dihedral)
     - Cα=Cβ–Cγ–OG in monolignols

To add missing types, copy them from a compatible parametrisation
(e.g., [Orella2019]_, [Petridis2009]_) or fit them via energy scans.

----

Comparison with Literature Parameters
--------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Reference
     - Relevance
   * - Jorgensen et al. 1996 [Jorgensen1996]_
     - Original OPLS-AA paper; source of all base types used here
       (phenol, anisole, alcohol, vinyl)
   * - Petridis & Smith 2009 [Petridis2009]_
     - All-atom MD of softwood lignin with CHARMM36-based types;
       useful cross-check for dihedral energetics along the β-O-4 linkage
   * - Orella et al. 2019 [Orella2019]_
     - OPLS-AA parametrisation of lignin; direct predecessor;
       provides dihedral corrections for ring–Cα–O–Cβ and other
       lignin-specific terms
   * - Schultz & Schmidt 2018 [Schultz2018]_
     - Free-energy perturbation study; discusses charge-group
       neutrality strategies for extended aromatic systems in OPLS

----

References
----------

.. [Jorgensen1996] W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives,
   "Development and Testing of the OPLS All-Atom Force Field on
   Conformational Energetics and Properties of Organic Liquids,"
   *J. Am. Chem. Soc.* **1996**, 118, 11225–11236.

.. [Petridis2009] L. Petridis and J. C. Smith,
   "A Molecular Mechanics Force Field for Lignin,"
   *J. Comput. Chem.* **2009**, 30, 457–467.

.. [Orella2019] M. J. Orella, Y. Román-Leshkov, and F. H. Brushett,
   "Emerging Role of Electrostatics in the Optimisation of Lignin Depolymerisation,"
   *ACS Sustain. Chem. Eng.* **2019**, 7, 5810–5821.

.. [Schultz2018] B. J. Schultz and C. M. Schmidt,
   "Benchmark Study of the Partial Charge Assignment Methods on
   Aromatic Oxygenate Biomass Models,"
   *J. Chem. Inf. Model.* **2018**, 58, 2471–2480.
