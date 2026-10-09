.. _opls_aa_parametrization:

==============================================================
OPLS-AA Parametrisation of Lignin Chains
==============================================================

This page documents the atom-type and partial-charge choices implemented in
``lignoforge/forcefield/opls.py``.  All types are taken from the GROMACS
``oplsaa.ff`` database [Jorgensen1996]_ and cross-checked against the TYR /
PHE entries of ``aminoacids.rtp``.  For the practical workflow see
:doc:`/user_guide/force_field`.

.. contents:: On this page
   :depth: 2
   :local:

----

Atom-type table
---------------

.. list-table::
   :header-rows: 1
   :widths: 22 30 14 12 22

   * - Atom(s)
     - Context
     - OPLS type
     - q (e)
     - Notes
   * - ``C1``
     - ipso C, sp³ side chain
     - ``opls_221``
     - −0.055
     - substituted aryl C
   * - ``C1``
     - ipso C, vinyl side chain
     - ``opls_221``
     - 0.000
     - self-neutral with vinyl C–H
   * - ``C2 C6``, bare ``C3 C5``
     - aromatic C–H
     - ``opls_145``
     - −0.115
     -
   * - ``H2 H3 H5 H6``
     - aromatic H
     - ``opls_146``
     - +0.115
     -
   * - ``C3 C5`` (OMe), ``C4`` (ether), ``C5`` (4-O-5 acceptor)
     - ring C bearing O
     - ``opls_199``
     - +0.085
     - anisole C
   * - ``C5`` (5-5, β-5)
     - ring C bearing C
     - ``opls_145``
     - 0.000
     - H5 removed
   * - ``C4``
     - free phenol
     - ``opls_166``
     - +0.150
     -
   * - ``O4H`` / ``HO4``
     - free phenol
     - ``opls_167`` / ``opls_168``
     - −0.585 / +0.435
     -
   * - ``O4H``, ``OM3``, ``OM5``
     - aryl ether O
     - ``opls_179``
     - −0.285
     -
   * - ``OG``
     - β-β furan ether O
     - ``opls_179``
     - −0.285
     - no ``HOG``
   * - ``CM3 CM5`` / ``HM3x HM5x``
     - methoxy
     - ``opls_181`` / ``opls_185``
     - +0.110 / +0.030
     -
   * - ``CA`` / ``HA``
     - sp³ benzylic alcohol
     - ``opls_219`` / ``opls_156``
     - +0.260 / +0.060
     - β-O-4 acceptor (``OA`` present)
   * - ``OA`` / ``HOA``, ``OG`` / ``HOG``
     - alcohol
     - ``opls_154`` / ``opls_155``
     - −0.683 / +0.418
     -
   * - ``CA`` / ``HA``
     - sp³ ether (α-O-4, β-5, β-β)
     - ``opls_183`` / ``opls_156``
     - +0.170 / +0.060
     -
   * - ``CB`` / ``HB``
     - sp³ ether (β-O-4 acceptor)
     - ``opls_183`` / ``opls_185``
     - +0.170 / +0.030
     - sums to +0.200 (see below)
   * - ``CB`` / ``HB``
     - sp³ CH, C–C linkage (β-5, β-β)
     - ``opls_137`` / ``opls_140``
     - −0.060 / +0.060
     -
   * - ``CB`` / ``HB1 HB2``
     - sp³ CH\ :sub:`2` (after α-O-4)
     - ``opls_136`` / ``opls_140``
     - −0.120 / +0.060
     -
   * - ``CA CB`` / ``HA HB``
     - vinyl (end units)
     - ``opls_142`` / ``opls_144``
     - −0.115 / +0.115
     -
   * - ``CG`` / ``HG1 HG2``
     - CH\ :sub:`2`–OH
     - ``opls_157`` / ``opls_156``
     - +0.145 / +0.060
     - 1-propanol H charge

Charge neutrality of ether bridges
----------------------------------

Turning a phenol into an aryl ether removes ``HO4`` and changes ``C4`` /
``O4H`` from +0.150 / −0.585 / +0.435 to +0.085 / −0.285, a net change of
−0.200 e.  The β-O-4 acceptor ``CB`` + ``HB`` (``opls_183`` +0.170,
``opls_185`` +0.030) change by +0.200 e, so every β-O-4 and α-O-4 bond is
exactly neutral.  The 4-O-5, 5-5, β-5 and β-β linkages are not neutral by
construction and rely on the per-residue renormalisation.

Per-residue renormalisation
---------------------------

Each residue is made exactly neutral by shifting all of its atoms by the same
amount (residual / n atoms), rounded to 10\ :sup:`-4` e.  This is preferred to
recomputing the charges: when the residual is modest it perturbs the polar
groups least, and it keeps the OPLS fragment values recognisable.  The
resulting changes are reported for every chain (:ref:`charge-report`).

Torsions absent from OPLS-AA
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Torsion
     - Treatment
   * - ``CA CA CM HC``
     - zero, as OPLS does for ethylbenzene ``CA CA CT HC``
   * - ``CA CT CT CA``
     - hydrocarbon values, as OPLS does for ``CA CT CT CT`` (β-5 coumaran)

All other bonds, angles and torsions of every generated chain (all linkages,
linear and branched) are found in ``ffbonded.itp``; this is checked by
``grompp`` in ``tools/validate_pipeline.py``.

References
----------

.. [Jorgensen1996] Jorgensen, W. L.; Maxwell, D. S.; Tirado-Rives, J.
   *J. Am. Chem. Soc.* **1996**, 118, 11225.
.. [Petridis2009] Petridis, L.; Smith, J. C. *J. Comput. Chem.* **2009**, 30, 457.
.. [Orella2019] Orella, M. J. et al. *ACS Sustain. Chem. Eng.* **2019**, 7, 9979.
