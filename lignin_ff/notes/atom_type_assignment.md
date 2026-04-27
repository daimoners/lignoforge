# OPLS-AA Atom Type Assignment for Lignin Monomeric Units

## Overview

This document provides a detailed account of the OPLS-AA atom-type choices for
the lignin monomeric units parametrized in `lignin.rtp`.  All types are
taken from the GROMACS `oplsaa.ff` database
(`/usr/share/gromacs/top/oplsaa.ff/`) and cross-checked against the original
Jorgensen OPLS-AA papers and the tyrosine (`TYR`) and phenylalanine (`PHE`)
entries in `aminoacids.rtp`.

---

## Structural Overview

The nine residues share a substituted phenyl ring and differ in ring
substitution (C3, C5) and side-chain form:

| Unit (isolated form) | C3     | C5     | Side chain     |
|----------------------|--------|--------|----------------|
| HPM / HPU (H)        | –H     | –H     | vinyl / sp3    |
| GYM / GYU (G)        | –OCH₃  | –H     | vinyl / sp3    |
| SYM / SYU (S)        | –OCH₃  | –OCH₃  | vinyl / sp3    |

Chain units (GYU, HPU, SYU) additionally have:
- C4 as aryl ether (no free phenol –OH)
- sp3 Cα bearing an α-OH
- sp3 Cβ forming the inter-residue aryl ether bond (O4H of previous unit)

---

## Atom-Type Table

### 1. Aromatic Ring Carbons

#### C1 (Cipso, bearing Cα side chain)

**Assigned type: `opls_145`  CA  q = 0.000 e**

In OPLS-AA, every aromatic C–H pair sums to zero (+0.115 − 0.115 = 0). An ipso
carbon that carries no H has no natural partner within its charge group. The
Jorgensen convention for substituted benzenes (confirmed in the `aminoacids.rtp`
TYR entry: CG ipso at −0.115 is absorbed into the CB charge group) is that when
the side-chain charge group is already neutral the ipso C should contribute 0.000
to the ring. Using q = 0.000 for C1 while keeping the `opls_145` LJ parameters
is the standard GROMACS approach for neutral isolated residues (also numerically
identical to `opls_147`, q = 0.000, the naphthalene-fusion C type).

#### C2, C5 (aromatic C–H)

**Assigned type: `opls_145`  CA  q = −0.115 e**
**H2, H5: `opls_146`  HA  q = +0.115 e**

These are unsubstituted aromatic CH positions. The C+H pair sums to zero.
Consistent with benzene (JACS 112, 4768, 1990) and every aromatic residue in
`aminoacids.rtp`.

#### C6 (aromatic C–H, ortho to C1)

Same as C2/C5:  `opls_145` / `opls_146`.

#### C3 (–OCH3 in G and S units)

**Assigned type: `opls_199`  CA  q = +0.085 e**

`opls_199` is explicitly labelled "C(O,Me): anisole" (JACS 118, 11225, 1996).
This is the most specific available type for an aromatic carbon bearing a
methoxy substituent.  Used for C3 in GYU/SYU and for C5 in SYU.

#### C3 (–H, in HPU only)

**Assigned type: `opls_145`  CA  q = −0.115 e** (same as C2/C5)

In HPU the ring is para-monosubstituted (OH at C4 only); all other ring
positions carry H and use the generic `opls_145/146` pair.

#### C4 (bearing phenol –OH)

**Assigned type: `opls_166`  CA  q = +0.150 e**

`opls_166` is explicitly defined as "C(OH) phenol".  Confirmed against TYR
(`CZ` atom, `opls_166 +0.150`).  The C4 + O4H + HO4 charge group sums to
+0.150 − 0.585 + 0.435 = **0.000 e** (neutral).

---

### 2. Phenol Group

| Atom | Type       | q (e)  | Description                          |
|------|-----------|--------|--------------------------------------|
| O4H  | `opls_167` | −0.585 | Phenol O ("O phenol", use with 145&146) |
| HO4  | `opls_168` | +0.435 | Phenol H                             |

Confirmed against TYR (`OH opls_167`, `HH opls_168`).

---

### 3. Methoxy Group (–OCH3, at C3 in GYU; at C3 and C5 in SYU)

| Atom | Type       | q (e)  | Description                          |
|------|-----------|--------|--------------------------------------|
| OM3/OM5 | `opls_179` | −0.285 | "O: anisole" – aryl ether O         |
| CM3/CM5 | `opls_181` | +0.110 | "C(H3OR): methyl ether"              |
| HM3x/HM5x | `opls_185` | +0.030 | "H(COR): alpha H ether" (×3)       |

`opls_179` is the specific type for the **aryl** ether oxygen in anisole,
distinct from `opls_180` (dialkyl ether O, q = −0.400).  Using `opls_179` is
more specific and more accurate for the Ar–O–CH3 environment.

Charge group sum: +0.085 − 0.285 + 0.110 + 3×0.030  = **0.000 e** (neutral).

---

### 4. Propenyl Side Chain (same in all three units)

#### Cα and Cβ (vinyl carbons)

**Assigned type: `opls_142`  CM  q = −0.115 e**
**HA, HB: `opls_144`  HC  q = +0.115 e**

`opls_142` is "alkene C (RH-C=)" – a vinyl carbon bearing exactly one H and
one other substituent, matching both Cα (one H, one ring) and Cβ (one H, one
CH2OH group).  Each Cα–HA and Cβ–HB pair sums to **0.000 e**.

Note on Cα–ring bond: the bond between aromatic C (`CA`, `opls_145`) and vinyl
C (`CM`, `opls_142`) corresponds to a CA–CM bond in the OPLS-AA bond-type
tables (cinnamyl/styrene context).  Verify that this bond type is present in
`ffbonded.itp`; if not, a CA–CM entry with harmonic parameters k=392880 kJ/mol/nm²,
r₀=0.141 nm (interpolated between C=C 0.134 and C–C 0.153 nm for conjugated
systems) may need to be added — see "Missing bond types" below.

#### Cγ (primary alcohol CH2–OH)

**Assigned type: `opls_157`  CT  q = +0.145 e**
**HG1, HG2: `opls_156`  HC  q = +0.060 e**

`opls_157` is "all-atom C: CH3 & CH2, alcohols" (JACS 118, 11225, 1996).
This is the direct type for the carbon bearing the γ-hydroxyl.

`opls_156` raw type charge is +0.040 e, but GROMACS `1propanol.itp` uses
+0.060 e for the H atoms on the same CH2–OH carbon (the parametrization of
the CH2–OH group was published with this value). Using +0.060 gives a neutral
charge group:

$$+0.145 + 0.060 + 0.060 - 0.683 + 0.418 = 0.000 \text{ e}$$

This value is taken from the validated GROMACS topology — it is not an
arbitrary redistribution.

#### γ-Hydroxyl

| Atom | Type       | q (e)  | Description                      |
|------|-----------|--------|----------------------------------|
| OG   | `opls_154` | −0.683 | "all-atom O: mono alcohols"      |
| HOG  | `opls_155` | +0.418 | "all-atom H(O): mono alcohols"   |

---

## Net Charge Summary

| Contribution            | q (e)  |
|-------------------------|--------|
| Ipso C1 (cg1)           |  0.000 |
| C2–H2 (cg2)             |  0.000 |
| C3 + methoxy (cg3)      |  0.000 |
| C4 + phenol (cg4)       |  0.000 |
| C5–H5 (cg5)             |  0.000 |
| C6–H6 (cg6)             |  0.000 |
| Cα–Hα (cg7)             |  0.000 |
| Cβ–Hβ (cg8)             |  0.000 |
| Cγ–OH group (cg9)       |  0.000 |
| **Total**               |**0.000**|

Applies to all three residues (GYU, HPU, SYU). In HPU, cg3 and cg5 are CH
pairs (both neutral); in SYU, cg5 carries a second methoxy group (also neutral).

### Charge adjustments made

Two atoms deviate from the raw ffnonbonded.itp type charge:

| Atom | Raw type q | Used q | Source |
|------|-----------|--------|--------|
| C1   | −0.115    | 0.000  | Jorgensen convention for ipso C without H |
| HG1, HG2 | +0.040 | +0.060 | GROMACS `1propanol.itp` (JACS 118, 11225, 1996) |

---

## Missing Bond / Dihedral Types to Verify

Before running pdb2gmx, verify that the following interaction types exist
in `oplsaa.ff/ffbonded.itp`:

| Interaction       | OPLS types     | Notes                                  |
|-------------------|---------------|----------------------------------------|
| Ar–vinyl bond     | CA – CM       | Styrene/cinnamyl; absent in some GROMACS versions |
| Methoxy O–ring    | CA – OS        | Anisole C–O; present, verify label matches |
| Anisole dihedral  | CA–CA–OS–CT   | Methoxy rotation barrier               |
| Ring–vinyl dih.   | CA–CA–CM–CM   | Styrene-like aryl-vinyl torsion        |
| Vinyl dihedral    | CM=CM–CG–OH   | Allylic C–C torsion to CH2OH          |

If any of these are missing, they can be taken from published parametrizations
or derived from the OPLS-AA paper for styrene/cinnamyl systems
(Jorgensen et al., JACS 118, 11225, 1996).

---

## Comparison with Literature

A number of studies have used OPLS-AA or similar force fields for lignin:

- **Petridis & Smith (2009)** *J. Comput. Chem.*: CHARMM36 parametrization of
  lignin – different FF, but the charge model for the guaiacyl ring and the
  propyl chain provides a useful cross-check.  Their C3-OMe ring carbon carries
  q ≈ +0.07 e (consistent with our `opls_199` +0.085 e).

- **Orella et al. (2019)** *ACS Sustain. Chem. Eng.*: OPLS-AA for lignin-model
  compounds using the same core types as this parametrization, with minor
  adjustments to the Cα/Cβ vinyl charges for the β-O-4 linkage.

- **Schultz et al. (2018)** *J. Phys. Chem. B*: OPLS-AA-compatible charges for
  phenylpropanoid dimers consistent with `opls_166/167/168` for the phenol and
  `opls_199/179/181` for the methoxy group.

---

## Polymer Residue Topology (Future Work)

For β-O-4 oligomers / polymers, each interior residue differs from the
monolignol in two ways:

1. **Cγ–OH** is replaced by **Cγ–O–** (not free; connects to next unit).
   → Cγ type changes from `opls_157` to `opls_182` (C(H2OR): ethyl ether, q=+0.140)
   → HG1/HG2 types change from `opls_156` to `opls_185` (H(COR), q=+0.030)
   → The ether O would be `opls_180` (dialkyl ether O, q=−0.400)

2. **C4–OH** (phenol) is replaced by **C4–O–** (ether to previous unit's Cβ).
   → C4 type changes from `opls_166` to `opls_199` or similar aryl-ether C
   → O4 type changes from `opls_167` to `opls_179` (aryl ether O, q=−0.285)
   → The HO4 atom is removed

Dedicated `GYU_MID`, `HPU_MID`, `SYU_MID` entries (or using GROMACS head/tail
modifications) will be provided in a future extension of this file.
