"""
Literature-sourced reference tables for structural prior estimation.

All values are derived from published experimental and NMR/GPC measurements.
Ranges are expressed as (mean, std) unless otherwise noted.

References (selected):
    - Rinaldi et al., Angew. Chem. 2016 (linkage distribution overview)
    - Ragauskas et al., Science 2014 (biomass composition)
    - Chio et al., Bioresour. Technol. 2019 (kraft / organosolv contrast)
    - Constant et al., Green Chem. 2016 (technical lignin MW)
    - Lupoi et al., Front. Bioeng. Biotechnol. 2015 (S/G ratios)
"""

from __future__ import annotations

# ── S:G:H molar fractions ─────────────────────────────────────────────────────
# Format: {biomass_type: {process_type | "native": (S_mean, S_std, G_mean, G_std, H_mean, H_std)}}
# Normalisation: S+G+H ~ 1 after sampling from independent normals.

SGH_BY_BIOMASS_PROCESS: dict = {
    "hardwood": {
        "native":        (0.55, 0.08, 0.43, 0.08, 0.02, 0.01),
        "kraft":         (0.48, 0.08, 0.50, 0.08, 0.02, 0.01),
        "organosolv":    (0.55, 0.07, 0.43, 0.07, 0.02, 0.01),
        "soda":          (0.50, 0.09, 0.48, 0.09, 0.02, 0.01),
        "sulfite":       (0.45, 0.09, 0.53, 0.09, 0.02, 0.01),
        "des":           (0.52, 0.07, 0.46, 0.07, 0.02, 0.01),
        "steam_explosion": (0.50, 0.08, 0.48, 0.08, 0.02, 0.01),
    },
    "softwood": {
        "native":        (0.02, 0.02, 0.97, 0.03, 0.01, 0.01),
        "kraft":         (0.02, 0.02, 0.95, 0.04, 0.03, 0.01),
        "organosolv":    (0.02, 0.02, 0.96, 0.03, 0.02, 0.01),
        "soda":          (0.02, 0.02, 0.94, 0.04, 0.04, 0.02),
        "sulfite":       (0.02, 0.02, 0.95, 0.03, 0.03, 0.01),
        "des":           (0.02, 0.02, 0.96, 0.03, 0.02, 0.01),
        "steam_explosion": (0.02, 0.02, 0.92, 0.05, 0.06, 0.02),
    },
    "agricultural_residue": {
        "native":        (0.25, 0.08, 0.55, 0.08, 0.20, 0.06),
        "kraft":         (0.22, 0.08, 0.55, 0.08, 0.23, 0.07),
        "organosolv":    (0.27, 0.07, 0.55, 0.07, 0.18, 0.06),
        "soda":          (0.23, 0.08, 0.53, 0.08, 0.24, 0.07),
        "sulfite":       (0.20, 0.07, 0.57, 0.08, 0.23, 0.07),
        "des":           (0.25, 0.07, 0.55, 0.07, 0.20, 0.06),
        "steam_explosion": (0.20, 0.08, 0.55, 0.08, 0.25, 0.07),
    },
    "mixed": {
        "native":        (0.35, 0.10, 0.55, 0.10, 0.10, 0.04),
        "kraft":         (0.30, 0.09, 0.57, 0.09, 0.13, 0.05),
        "organosolv":    (0.35, 0.09, 0.55, 0.09, 0.10, 0.04),
        "soda":          (0.32, 0.09, 0.55, 0.09, 0.13, 0.05),
        "sulfite":       (0.30, 0.09, 0.57, 0.09, 0.13, 0.05),
        "des":           (0.33, 0.09, 0.55, 0.09, 0.12, 0.04),
        "steam_explosion": (0.30, 0.10, 0.55, 0.10, 0.15, 0.05),
    },
}

# ── Linkage distributions ─────────────────────────────────────────────────────
# Format: {biomass_type: {process_type: {linkage_name: (mean, std)}}}
# Linkage order: beta-O-4, 4-O-5, alpha-O-4, 5-5, beta-5, beta-beta, beta-1

LINKAGE_BY_BIOMASS_PROCESS: dict = {
    "hardwood": {
        "native": {
            "beta-O-4":   (0.60, 0.05),
            "4-O-5":      (0.04, 0.02),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.04, 0.02),
            "beta-5":     (0.06, 0.02),
            "beta-beta":  (0.12, 0.03),
            "beta-1":     (0.08, 0.02),
        },
        "kraft": {
            "beta-O-4":   (0.40, 0.07),   # kraft cleaves beta-O-4
            "4-O-5":      (0.08, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.12, 0.04),   # condensed C-C bonds increase
            "beta-5":     (0.14, 0.04),
            "beta-beta":  (0.12, 0.03),
            "beta-1":     (0.09, 0.03),
        },
        "organosolv": {
            "beta-O-4":   (0.55, 0.06),
            "4-O-5":      (0.05, 0.02),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.07, 0.02),
            "beta-5":     (0.08, 0.02),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.08, 0.02),
        },
        "soda": {
            "beta-O-4":   (0.45, 0.07),
            "4-O-5":      (0.06, 0.02),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.10, 0.03),
            "beta-5":     (0.12, 0.03),
            "beta-beta":  (0.12, 0.03),
            "beta-1":     (0.10, 0.03),
        },
        "sulfite": {
            "beta-O-4":   (0.35, 0.07),
            "4-O-5":      (0.09, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.14, 0.04),
            "beta-5":     (0.15, 0.04),
            "beta-beta":  (0.13, 0.03),
            "beta-1":     (0.09, 0.03),
        },
        "des": {
            "beta-O-4":   (0.58, 0.06),
            "4-O-5":      (0.05, 0.02),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.06, 0.02),
            "beta-5":     (0.07, 0.02),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.07, 0.02),
        },
        "steam_explosion": {
            "beta-O-4":   (0.42, 0.07),
            "4-O-5":      (0.07, 0.02),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.12, 0.03),
            "beta-5":     (0.14, 0.04),
            "beta-beta":  (0.12, 0.03),
            "beta-1":     (0.08, 0.03),
        },
    },
    "softwood": {
        "native": {
            "beta-O-4":   (0.50, 0.06),
            "4-O-5":      (0.08, 0.03),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.10, 0.03),
            "beta-5":     (0.12, 0.04),
            "beta-beta":  (0.08, 0.03),
            "beta-1":     (0.06, 0.02),
        },
        "kraft": {
            "beta-O-4":   (0.30, 0.07),
            "4-O-5":      (0.12, 0.04),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.18, 0.05),
            "beta-5":     (0.18, 0.05),
            "beta-beta":  (0.10, 0.03),
            "beta-1":     (0.07, 0.03),
        },
        "organosolv": {
            "beta-O-4":   (0.48, 0.06),
            "4-O-5":      (0.08, 0.03),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.12, 0.03),
            "beta-5":     (0.12, 0.03),
            "beta-beta":  (0.08, 0.03),
            "beta-1":     (0.06, 0.02),
        },
        "soda": {
            "beta-O-4":   (0.35, 0.07),
            "4-O-5":      (0.10, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.15, 0.04),
            "beta-5":     (0.17, 0.04),
            "beta-beta":  (0.10, 0.03),
            "beta-1":     (0.08, 0.03),
        },
        "sulfite": {
            "beta-O-4":   (0.25, 0.07),
            "4-O-5":      (0.14, 0.04),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.20, 0.05),
            "beta-5":     (0.20, 0.05),
            "beta-beta":  (0.10, 0.03),
            "beta-1":     (0.06, 0.03),
        },
        "des": {
            "beta-O-4":   (0.47, 0.06),
            "4-O-5":      (0.09, 0.03),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.11, 0.03),
            "beta-5":     (0.13, 0.03),
            "beta-beta":  (0.08, 0.03),
            "beta-1":     (0.06, 0.02),
        },
        "steam_explosion": {
            "beta-O-4":   (0.32, 0.07),
            "4-O-5":      (0.11, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.17, 0.04),
            "beta-5":     (0.18, 0.04),
            "beta-beta":  (0.10, 0.03),
            "beta-1":     (0.07, 0.03),
        },
    },
    "agricultural_residue": {
        "native": {
            "beta-O-4":   (0.55, 0.06),
            "4-O-5":      (0.05, 0.02),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.06, 0.02),
            "beta-5":     (0.10, 0.03),
            "beta-beta":  (0.10, 0.03),
            "beta-1":     (0.08, 0.02),
        },
        "kraft": {
            "beta-O-4":   (0.38, 0.07),
            "4-O-5":      (0.09, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.13, 0.04),
            "beta-5":     (0.15, 0.04),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.09, 0.03),
        },
        "organosolv": {
            "beta-O-4":   (0.53, 0.06),
            "4-O-5":      (0.06, 0.02),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.08, 0.02),
            "beta-5":     (0.10, 0.03),
            "beta-beta":  (0.10, 0.03),
            "beta-1":     (0.07, 0.02),
        },
        "soda": {
            "beta-O-4":   (0.42, 0.07),
            "4-O-5":      (0.07, 0.02),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.12, 0.03),
            "beta-5":     (0.14, 0.04),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.09, 0.03),
        },
        "sulfite": {
            "beta-O-4":   (0.32, 0.07),
            "4-O-5":      (0.10, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.16, 0.04),
            "beta-5":     (0.17, 0.04),
            "beta-beta":  (0.12, 0.03),
            "beta-1":     (0.08, 0.03),
        },
        "des": {
            "beta-O-4":   (0.55, 0.06),
            "4-O-5":      (0.05, 0.02),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.07, 0.02),
            "beta-5":     (0.09, 0.03),
            "beta-beta":  (0.10, 0.03),
            "beta-1":     (0.08, 0.02),
        },
        "steam_explosion": {
            "beta-O-4":   (0.40, 0.07),
            "4-O-5":      (0.08, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.13, 0.03),
            "beta-5":     (0.15, 0.04),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.08, 0.03),
        },
    },
    "mixed": {
        "native": {
            "beta-O-4":   (0.55, 0.06),
            "4-O-5":      (0.06, 0.02),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.07, 0.02),
            "beta-5":     (0.09, 0.03),
            "beta-beta":  (0.10, 0.03),
            "beta-1":     (0.07, 0.02),
        },
        "kraft": {
            "beta-O-4":   (0.38, 0.07),
            "4-O-5":      (0.10, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.14, 0.04),
            "beta-5":     (0.15, 0.04),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.07, 0.03),
        },
        "organosolv": {
            "beta-O-4":   (0.52, 0.06),
            "4-O-5":      (0.07, 0.02),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.08, 0.02),
            "beta-5":     (0.09, 0.03),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.07, 0.02),
        },
        "soda": {
            "beta-O-4":   (0.43, 0.07),
            "4-O-5":      (0.08, 0.02),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.12, 0.03),
            "beta-5":     (0.13, 0.03),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.08, 0.03),
        },
        "sulfite": {
            "beta-O-4":   (0.33, 0.07),
            "4-O-5":      (0.11, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.17, 0.04),
            "beta-5":     (0.17, 0.04),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.06, 0.03),
        },
        "des": {
            "beta-O-4":   (0.54, 0.06),
            "4-O-5":      (0.06, 0.02),
            "alpha-O-4":  (0.06, 0.02),
            "5-5":        (0.07, 0.02),
            "beta-5":     (0.09, 0.03),
            "beta-beta":  (0.10, 0.03),
            "beta-1":     (0.08, 0.02),
        },
        "steam_explosion": {
            "beta-O-4":   (0.40, 0.07),
            "4-O-5":      (0.09, 0.03),
            "alpha-O-4":  (0.05, 0.02),
            "5-5":        (0.13, 0.03),
            "beta-5":     (0.14, 0.04),
            "beta-beta":  (0.11, 0.03),
            "beta-1":     (0.08, 0.03),
        },
    },
}

# ── Molecular weight (Mn, Mw in g/mol) ───────────────────────────────────────
# Format: {process_type: (Mn_mean, Mn_std, Mw_mean, Mw_std)}

MW_BY_PROCESS: dict = {
    "kraft":          (2000, 500,  5500,  1300),
    "soda":           (2500, 600,  6000,  1400),
    "organosolv":     (3000, 700,  7000,  1500),
    "sulfite":        (1500, 400,  4000,  1000),
    "des":            (3500, 800,  8000,  1800),
    "steam_explosion":(1200, 350,  3500,   900),
}

# ── Branching index ───────────────────────────────────────────────────────────
# Fraction of monomers with ≥3 connections.
# Format: {process_type: (mean, std)}

BRANCHING_BY_PROCESS: dict = {
    "kraft":          (0.35, 0.07),
    "soda":           (0.30, 0.06),
    "organosolv":     (0.18, 0.05),
    "sulfite":        (0.40, 0.08),
    "des":            (0.15, 0.04),
    "steam_explosion":(0.38, 0.08),
}

# ── Condensation degree ───────────────────────────────────────────────────────
# Ratio of C-C bonds to total inter-monomer bonds.
# Format: {process_type: (mean, std)}

CONDENSATION_BY_PROCESS: dict = {
    "kraft":          (0.62, 0.08),
    "soda":           (0.55, 0.07),
    "organosolv":     (0.30, 0.06),
    "sulfite":        (0.70, 0.08),
    "des":            (0.25, 0.05),
    "steam_explosion":(0.60, 0.08),
}

# ── Average monomer molecular weights (g/mol) ──────────────────────────────────
MONOMER_MW = {"H": 150.17, "G": 180.20, "S": 210.23}
