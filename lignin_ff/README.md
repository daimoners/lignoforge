# OPLS-AA Force-Field Parameters for Lignin — GROMACS

This directory contains OPLS-AA force-field parameter files for all-atom
molecular dynamics simulations of lignin with GROMACS, together with the
Python tooling that assigns per-atom types to any PDB produced by LignoForge.

---

## Directory Layout

```
lignin_ff/
├── README.md                    ← this file
├── lignin.rtp                   ← GROMACS residue topology (.rtp)
│                                   9 residue definitions (see below)
├── residuetypes_lignin.dat      ← pdb2gmx residuetype entries
├── notes/
│   └── atom_type_assignment.md  ← per-atom type justification & charge accounting
└── tools/
    ├── assign_chain_types.py    ← assigns OPLS-AA types to a LignoForge PDB
    └── test_assign.py           ← validation tests for assign_chain_types
```

---

## Residues Defined

### β-O-4 internal chain units (sp3 side chain, two inter-residue bonds)

| GROMACS name | Monomer type        | Description                                       |
|-------------|---------------------|---------------------------------------------------|
| `GYU`       | G (guaiacyl)        | –OCH₃ at C3; sp3 α-OH (CA) + α-ether (CB)        |
| `HPU`       | H (p-hydroxyphenyl) | no OMe; sp3 α-OH (CA) + α-ether (CB)             |
| `SYU`       | S (syringyl)        | –OCH₃ at C3 and C5; sp3 α-OH (CA) + α-ether (CB) |

Each chain unit carries two inter-residue bonds in `[ bonds ]`:

```
O4H  -CB     ; aryl ether O of this unit → CB of the previous residue
CB   +O4H    ; CB of this unit → aryl ether O of the next residue
```

### Isolated free monolignols (vinyl side chain, free phenol — no inter-residue bonds)

| GROMACS name | Compound            | Formula    | MW (g/mol) |
|-------------|---------------------|------------|------------|
| `GYM`       | coniferyl alcohol   | C₁₀H₁₂O₃  | 180.20     |
| `HPM`       | p-coumaryl alcohol  | C₉H₁₀O₂   | 150.17     |
| `SYM`       | sinapyl alcohol     | C₁₁H₁₄O₄  | 210.23     |

### Neutral saturated references (sp3 side chain, free phenol — no inter-residue bonds)

| GROMACS name | Compound                        |
|-------------|----------------------------------|
| `GNM`       | dihydroconiferyl alcohol (G)     |
| `HNM`       | dihydro-p-coumaryl alcohol (H)   |
| `SNM`       | dihydrosinapyl alcohol (S)       |

**Net charge:** 0.000 e for interior units (GYU/HPU/SYU) and all isolated monomers. Terminal β-O-4 head/tail units carry ±0.240 e that cancel at chain level (total chain = 0.000 e).

---

## Force Field

**Base**: OPLS-AA (Jorgensen et al.)  
**Implementation**: GROMACS `oplsaa.ff` (tested with GROMACS 2022/2023)  
**Location**: `/usr/share/gromacs/top/oplsaa.ff/`

All atom types and LJ parameters are taken directly from `ffnonbonded.itp`.
Two minor charge corrections are applied to achieve neutral charge groups;
see `notes/atom_type_assignment.md` for the full per-atom justification.

---

## Key OPLS-AA Types

| Group                         | Position                 | Type       | q (e)        | Description                                    |
|-------------------------------|--------------------------|------------|--------------|------------------------------------------------|
| Aromatic ring                 | Cipso (no H, C1)         | `opls_221` | −0.055 (fixed) | Raw OPLS-AA charge; same for all residue types |
|                               | C–H                      | `opls_145` | −0.115       | Benzene C                                      |
|                               | H on C–H                 | `opls_146` | +0.115       | Benzene H                                      |
|                               | C–OH (free phenol)       | `opls_166` | +0.150       | Phenol ring C                                  |
|                               | C–O– (aryl ether, terminal) | `opls_199` | +0.085    | Standard aryl ether C / anisole C              |
|                               | C–O– (aryl ether, interior) | `opls_199` | +0.095    | Adjusted +0.010 for interior units             |
| Free phenol –OH               | O                        | `opls_167` | −0.585       | Phenol O                                       |
|                               | H                        | `opls_168` | +0.435       | Phenol H                                       |
| Aryl ether –O– (β-O-4)       | O (terminal)             | `opls_179` | −0.285       | Standard aryl ether O                          |
|                               | O (interior)             | `opls_179` | −0.275       | Adjusted +0.010 for interior units             |
| Aryl ether –O– (OMe)         | O                        | `opls_179` | −0.285       | Methoxy / anisole ether O                      |
| Methoxy –OCH₃                 | CH₃                      | `opls_181` | +0.110       | Methyl ether C                                 |
|                               | H on CH₃                 | `opls_185` | +0.030       | Ether α-H                                      |
| Vinyl chain (monolignols)     | Cα sp² (isolated)        | `opls_142` | −0.115       | Standard alkene C                              |
|                               | Cα sp² (chain tail donor) | `opls_142` | −0.105      | Adjusted so tail net = −0.240 e                |
|                               | Hα (isolated)            | `opls_144` | +0.115       | Standard vinyl H                               |
|                               | Hα (chain tail donor)    | `opls_144` | +0.120       | Adjusted for tail neutrality                   |
|                               | Cβ sp²                   | `opls_142` | −0.115       | Alkene C                                       |
|                               | Hβ                       | `opls_144` | +0.115       | Vinyl H                                        |
| sp³ chain (chain units)       | CA – secondary alc. (Ar–CHOH) | `opls_219` | +0.260  | Benzyl-alcohol C; most specific for Cα         |
|                               | CB – sp³ ether (–CH–O–Ar) | `opls_183` | +0.180      | Isopropyl-ether C; +0.180 (adjusted from OPLS raw) |
|                               | HB on CB ether (interior) | `opls_185` | +0.040      | H on sp³ ether C (interior units)              |
|                               | HB on CB ether (terminal) | `opls_185` | +0.060      | H on sp³ ether C (terminal units)              |
|                               | HA on Cα                 | `opls_156` | +0.060       | H on sp³ alcohol C                             |
| γ-Alcohol                     | CG                       | `opls_157` | +0.145       | sp³ CH₂–OH carbon                             |
|                               | HG1, HG2 (interior)      | `opls_156` | +0.040       | H on Cγ — adjusted for interior units          |
|                               | HG1, HG2 (terminal)      | `opls_156` | +0.060       | H on Cγ — standard/terminal value              |
|                               | OG, OA                   | `opls_154` | −0.683       | Aliphatic alcohol O                            |
|                               | HOG, HOA                 | `opls_155` | +0.418       | Aliphatic alcohol H                            |

---

## Type-Assignment Tool

For PDB structures generated by LignoForge, OPLS-AA types must account for
which linkages are present — e.g. C4 changes from phenol (`opls_166`) to
aryl ether (`opls_199`) when involved in a β-O-4 bond.
`tools/assign_chain_types.py` handles this automatically for all 7 lignin
linkage types (β-O-4, 5-5, 4-O-5, β-5, β-β, α-O-4, β-1).

### Usage

```bash
# Assign types; write custom per-chain RTP and a residue-renamed PDB
python lignin_ff/tools/assign_chain_types.py chain.pdb

# Explicit output paths
python lignin_ff/tools/assign_chain_types.py chain.pdb \
    -o chain_custom.rtp \
    --renamed-pdb chain_renamed.pdb
```

The script:
1. Reads the PDB and detects all inter-residue bonds.
2. Identifies each linkage type from the bonded atom names.
3. Assigns OPLS-AA types from the template residues in `lignin.rtp`, then
   applies linkage-specific modifications from the built-in `LINKAGE_MODS` table.
4. Verifies that each residue has net charge exactly 0.000 e.
5. Writes a per-chain custom RTP for use with GROMACS `pdb2gmx`.
6. Writes an optional renamed PDB with ≤3-char GROMACS residue names.

### Running the tests

```bash
python lignin_ff/tools/test_assign.py
# or
python -m pytest lignin_ff/tools/test_assign.py -v
```

13 tests cover isolated H/G/S monomers, GGG/HH/SSS β-O-4 chains, 5-5,
4-O-5, and β-1 dimers. All tests print `ALL TESTS PASSED` on success.

---

## Installation into GROMACS

```bash
# System-wide (requires write access)
cp lignin_ff/lignin.rtp /usr/share/gromacs/top/oplsaa.ff/
cat lignin_ff/residuetypes_lignin.dat >> /usr/share/gromacs/top/residuetypes.dat
```

Project-local force field (recommended for reproducibility):

```bash
cp -r /usr/share/gromacs/top/oplsaa.ff ./oplsaa.ff
cp lignin_ff/lignin.rtp oplsaa.ff/
cat lignin_ff/residuetypes_lignin.dat >> oplsaa.ff/residuetypes.dat
```

Running pdb2gmx with the custom chain topology:

```bash
gmx pdb2gmx -f chain_renamed.pdb -ff ./oplsaa -water tip3p \
    -o processed.gro -p topol.top \
    -rtpres chain_custom.rtp
```

---

## PDB Atom-Name Conventions

| Position               | Atom names                         |
|------------------------|------------------------------------|
| Ring                   | `C1 C2 H2 C3 C4 C5 C6`            |
| Free phenol            | `O4H HO4`                         |
| OMe at C3              | `OM3 CM3 HM31 HM32 HM33`          |
| OMe at C5 (S-type)     | `OM5 CM5 HM51 HM52 HM53`          |
| Vinyl (monolignols)    | `CA HA CB HB`                      |
| sp3 chain (chain units)| `CA HA OA HOA CB HB`               |
| γ-alcohol              | `CG HG1 HG2 OG HOG`               |

---

## References

1. Jorgensen, W. L.; Maxwell, D. S.; Tirado-Rives, J.  
   *Development and Testing of the OPLS All-Atom Force Field.*  
   J. Am. Chem. Soc. **1996**, 118, 11225–11236.

2. Jorgensen, W. L.; Schyman, P.  
   *Treatment of Halogen Bonding in OPLS-AA.*  
   J. Chem. Theory Comput. **2012**, 8, 3895–3901.  
   (aromatic parameters: JACS **1990**, 112, 4768)

3. Petridis, L.; Smith, J. C.  
   *A molecular mechanics force field for lignin.*  
   J. Comput. Chem. **2009**, 30, 457–467.

4. Orella, M. J. et al.  
   *Lignin-KMC: A Toolkit for Simulating Lignin Biosynthesis.*  
   ACS Sustain. Chem. Eng. **2019**, 7, 9979.
