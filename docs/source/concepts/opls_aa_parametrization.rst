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
     - ipso aromatic C (no H, bears propanoid chain)
     - ``opls_221``
     - CA
     - see note
     - Substituted aryl C; q adjusted per residue context by ``balance_charges()``;
       see :ref:`charge-corrections`
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
     - sp³ α-C of chain unit (bears α-OH; Ar–CHOH–)
     - ``opls_219``
     - CT(benzyl-alc.)
     - +0.260
     - Benzyl-alcohol type [Jorgensen1996]_; most specific for secondary alcohol α to aryl
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
     - Sum CA+HA+OA+HOA = +0.055 e; absorbed into C1
   * - ``CB``
     - sp³ β-C (bears inter-residue β-O-4 ether O; –CH–O–Ar)
     - ``opls_183``
     - CT(i-Pr ether)
     - +0.170
     - Isopropyl-ether type [Jorgensen1996]_; most specific for sp³ C bearing one ether O
   * - ``HB``
     - sp³ β-H
     - ``opls_156``
     - HC(alcohol)
     - +0.060
     - From 2-propanol model; CB+HB = +0.230 e; absorbed into C1
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
     - C1 (ipso, ``opls_221``)
     - −0.085 / −0.055 / 0.000
     - GYU/HPU/SYU / HNM/GNM/SNM / GYM/HPM/SYM; C1 absorbs sp³ side-chain imbalance
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
     - CA HA (vinyl) or CA HA OA HOA (sp³, ``opls_219``)
     - 0.000 (vinyl) / +0.055 (sp³, absorbed into C1)
     - monolignols / chain units + ref. monomers
   * - 8
     - CB HB (vinyl or sp³ ether, ``opls_183``)
     - 0.000 (vinyl) / +0.230 (sp³, absorbed into C1)
     - all residues
   * - 9
     - CG HG1 HG2 OG HOG
     - 0.000
     - all residues (all carry γ-OH)

All nine residues verified: **net charge = 0.000 e** (individual cgnr 1, 7, 8 may
be non-zero in chain units; the imbalance is absorbed by C1 via ``balance_charges()``).

----

.. _charge-corrections:

Charge Corrections
------------------

**C1 (ipso) → type ``opls_221``, q context-dependent**

``opls_221`` is the OPLS-AA type for a **substituted aryl carbon** (no attached H),
with a standard database charge of −0.055 e and LJ parameters identical to
``opls_145`` (aromatic CH).  In lignin, C1 is always trisubstituted (ring bond ×2
plus propanoid chain), making ``opls_221`` more specific than ``opls_145``.

Because the sp³ carbons in the propanoid side chain (Cα ``opls_219``, Cβ ``opls_183``)
carry charges that differ from the vinyl types used in monolignols, a small
per-residue charge imbalance arises that must be absorbed by C1 to maintain
per-residue neutrality.  The function ``balance_charges()`` in
``assign_chain_types.py`` computes this at run time:

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - Residue context
     - C1 charge (e)
     - Explanation
   * - Chain units: GYU / HPU / SYU
     - −0.085
     - Cα (+0.260) + Cβ (+0.170) vs. vinyl (−0.115 each): net imbalance +0.085
   * - Ref. monomers: HNM / GNM / SNM
     - −0.055
     - Equals raw ``opls_221`` charge; self-consistent (no Cβ ether, Cα alone)
   * - Monolignols: GYM / HPM / SYM
     - 0.000
     - Vinyl side chain (``opls_142``, q = −0.115 + H +0.115 = 0); no imbalance

**HG1, HG2 (γ-methylene H) → q = +0.060 e**

The raw OPLS-AA generic CH₂ H type carries +0.040 e.  The GROMACS reference
``1propanol.itp`` uses +0.060 e for the methylene adjacent to the hydroxyl-bearing
carbon (the –CH₂OH group).  Since Cγ in all lignin residues is –CH₂OH, q(HG) =
+0.060 e is used throughout.

Combined effect: charge group 9 (Cγ–OH) sums to:

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

   CA  opls_219  (CT, benzyl-alcohol secondary C; Ar–CHOH–)  q = +0.260
   HA  opls_156  (HC)                                        q = +0.060
   CB  opls_183  (CT, i-Pr ether secondary C; –CH–O–Ar)     q = +0.170  [chain units only]
   HB  opls_156  (HC)                                        q = +0.060

``opls_219`` (benzyl-alcohol type) is more specific than the generic
secondary-alcohol type ``opls_157`` for Cα because Cα is adjacent to the
aromatic ring.  ``opls_183`` (isopropyl-ether type) is the most specific
available OPLS-AA type for a secondary sp³ carbon bearing one ether oxygen.
Both types share the same LJ parameters as ``opls_157`` (σ = 3.50 Å,
ε = 0.276 kJ/mol); only the partial charges differ.

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
