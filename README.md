# LignoForge

> **Top-down stochastic generation of lignin structural models**  
> From experimental constraints to atomistic polymer ensembles

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![Documentation](https://img.shields.io/badge/docs-readthedocs-green.svg)](https://lignoforge.readthedocs.io/en/latest/)

---

## Overview

LignoForge is a Python framework that generates statistically valid ensembles of lignin polymer structures directly from experimental characterisation data (NMR, GPC, elemental analysis). The workflow proceeds entirely without requiring atomic-resolution input: a JSON file describing the biomass type, extraction process, and any known compositional data is sufficient to produce a library of atomistic chains.

**Workflow summary:**

```
JSON input (biomass type, process, optional constraints)
        │
        ▼
Prior estimation  ←  curated literature database (S/G/H, linkages, MW)
        │
        ▼
Stochastic simulation  ←  hierarchical kMC / MCMC algorithm
        │
        ▼
Population of lignin polymers  →  SMILES · SDF · PDB · JSON topologies · HTML viewer
```

**Key features:**
- Literature-backed prior estimation for 4 biomass types × 7 extraction processes
- Hierarchical Markov-Chain Monte Carlo (hMMC) that matches target S/G/H and linkage distributions
- Full atomistic 3-D coordinate generation via RDKit (ETKDG + MMFF94/UFF)
- PDB export with PDB-compatible residue codes (`GYU` / `SYU` / `HPU`)
- Coarse-grained topology JSON + self-contained interactive HTML viewer
- Reproducible runs via integer seeds
- JSON Schema-validated input

---

## Installation

### Requirements

- Python ≥ 3.9
- RDKit (see note below)

### With pip

```bash
pip install lignoforge
```

### From source (recommended for development)

```bash
git clone https://github.com/daimoners/lignoforge.git
cd lignoforge
pip install -e .
```

### RDKit

RDKit is required for 3-D coordinate generation and SDF/PDB export.  
The easiest installation route is via conda:

```bash
conda install -c conda-forge rdkit
```

Or via pip (PyPI wheel):

```bash
pip install rdkit
```

### Full dependency list

| Package | Minimum version |
|---------|-----------------|
| networkx | ≥ 2.5 |
| numpy | ≥ 1.19 |
| scipy | ≥ 1.3 |
| pandas | ≥ 0.25 |
| matplotlib | ≥ 3.1 |
| pysmiles | ≥ 1.0.1 |
| rdkit | ≥ 2021.09 |
| jsonschema | ≥ 4.0 |

---

## Quick Start

### Three-line pipeline

```python
from lignoforge.pipeline import LigninPipeline

pipeline = LigninPipeline.from_json("examples/lignin_input_example.json")
results  = pipeline.run()
```

### With a Python dictionary

```python
from lignoforge.pipeline import LigninPipeline

results = LigninPipeline.from_dict({
    "material_origin": {
        "lignin_type": "kraft_lignin",
        "biomass_type": "softwood"
    },
    "extraction_process": {
        "process_type": "kraft"
    }
}).run()
```

### With explicit constraints

```python
results = LigninPipeline.from_dict({
    "material_origin": {
        "biomass_type": "hardwood",
        "S_fraction": 0.52,
        "G_fraction": 0.45,
        "H_fraction": 0.03
    },
    "extraction_process": {"process_type": "organosolv"},
    "molecular_weight": {
        "Mn": 1200,
        "Mw": 3000,
        "PDI": 2.5
    },
    "simulation_config": {
        "n_population": 20
    }
}).run()
```

---

## Command-Line Tools

LignoForge ships a standalone command-line tool for building individual chains without writing any Python code.

### `lignoforge-chain`

Build one or more lignin chains by passing a JSON input file and morphology flags directly on the command line.  The JSON input is optional — when omitted, built-in hardwood-kraft defaults are used.

```bash
# Minimal: 10-monomer chain, SMILES + PDB output
lignoforge-chain examples/lignin_input_example.json --n-monomers 10

# No JSON needed — uses built-in defaults
lignoforge-chain --n-monomers 10

# Target MW instead of monomer count
lignoforge-chain input.json --mw-target 3000

# Pure softwood G-type, β-O-4 rich, 3 independent chains
lignoforge-chain input.json --monomer-type G --beta-O-4 0.80 --n-chains 3 --seed 7

# All output formats, 4 parallel workers
lignoforge-chain input.json --n-monomers 15 --format all --n-workers 4 --output my_chain/

# Demo-equivalent: full pipeline output, statistics, priors, 2 chains
lignoforge-chain examples/lignin_input_example.json \
  --n-chains 2 --seed 42 \
  --format all \
  --population-stats --export-priors --export-sim-params \
  --n-workers 4 --output demo_cli/
```

**Key options:**

| Option | Description |
|--------|-------------|
| `--n-monomers N` | Target degree of polymerisation |
| `--mw-target MW` | Approximate target MW in g/mol (auto-converts to DP) |
| `--monomer-type H\|G\|S\|mix` | Force a single monomer type |
| `--S-fraction F`, `--G-fraction F`, `--H-fraction F` | Override S/G/H fractions (0–1) |
| `--beta-O-4 F`, `--5-5 F`, … | Override individual linkage fractions (re-normalised) |
| `--branching F` | Branching propensity per MC step (0 = linear) |
| `--n-chains N` | Number of independent chains (default: 1) |
| `--seed N` | Random seed (default: 42) |
| `--n-workers N` | Parallel CPU workers for 3-D embedding (default: all cores; `1` = sequential) |
| `--format FMT` | `smiles`, `population-smiles`, `sdf`, `pdb`, `pdb-cg`, `json-atomistic`, `json-cg`, `html`, `all` (default: `smiles,pdb`) |
| `--no-3d` | Skip 3-D coordinate generation entirely |
| `--no-optimize` | Embed 3-D coordinates but skip MMFF/UFF optimisation (faster) |
| `--no-H` | Strip explicit hydrogens from 3-D output |
| `--export-priors` | Write `*_estimated_priors.json` |
| `--export-sim-params` | Write `*_simulation_parameters.json` |
| `--population-stats` | Write `*_population_stats.json` (ensemble mean/std/min/max) |
| `--output DIR` | Output directory (default: `chain_output/`) |

Each run writes per-chain output files, a `<name>_chain_stats.json` with full characterisation metrics, and a `<name>_manifest.json` recording all run parameters for reproducibility.

Run `lignoforge-chain --help` for the full option list.

---

## Running the Demo

A ready-to-run demo script is provided in `examples/`:

```bash
# Quick demo: 2 chains (fast)
python examples/demo_run.py

# Custom number of chains and seed
python examples/demo_run.py --n-chains 10 --seed 42

# Custom output directory
python examples/demo_run.py --n-chains 5 --output my_results/

# Limit parallel workers (default: all CPU cores)
python examples/demo_run.py --n-chains 20 --workers 4
```

The script loads `examples/lignin_input_example.json` (hardwood kraft lignin), runs the full pipeline,
and prints a summary of all generated artefacts.

> **Note on runtime:** 3-D coordinate generation with RDKit runs in parallel across all available
> CPU cores by default. On an 8-core machine, generating 20 chains takes roughly the same wall-clock
> time as generating 2-3 chains sequentially. Use `--workers N` to cap the number of processes.
> For very large populations (`n_population > 100`), consider disabling 3-D generation during
> initial exploration and enabling it only for the final export.

---

## Output Files

Running the pipeline creates the following files in the output directory:

| File | Description |
|------|-------------|
| `input_high_level.json` | Validated copy of the input data |
| `estimated_priors.json` | Structural priors estimated from the literature database |
| `simulation_parameters.json` | Translated simulation kwargs passed to the kMC engine |
| `population.smi` | One SMILES string per chain |
| `population.sdf` | SDF file with all chains (requires 3-D enabled) |
| `chain_statistics.json` | Per-chain characterisation metrics (MW, DP, linkage fracs, …) |
| `population_statistics.json` | Ensemble mean / std / min / max for every metric |
| `atomistic_topology.json` | Full atom-level topology with 3-D coordinates (one record per atom) |
| `coarse_grained_topology.json` | Bead-level topology (one bead = one monomer) with 3-D coordinates derived as the centre-of-mass of each monomer's atoms |
| `coarse_grained_topology_viewer.html` | Self-contained interactive 3-D Plotly viewer; beads coloured by type (H/G/S), bonds coloured by linkage |
| `pdb/polymer_A.pdb`, … | Atomistic PDB file per chain (HETATM + CONECT records) |
| `pdb_cg/polymer_cg_A.pdb`, … | Coarse-grained PDB (one pseudo-atom `BB` per monomer, VMD-compatible) |

---

## Visualisation

### PDB structures (atomistic)

PDB files are exported with HETATM records and standard three-letter residue codes:

| Monomer | Residue code |
|---------|-------------|
| *p*-Hydroxyphenyl (H) | `HPU` |
| Guaiacyl (G)          | `GYU` |
| Syringyl (S)          | `SYU` |

Open PDB files with any standard molecular viewer:

| Tool | Command / note |
|------|---------------|
| **PyMOL** | `pymol pdb/chain_A.pdb` |
| **UCSF ChimeraX** | `open pdb/chain_A.pdb` |
| **VMD** | `vmd pdb/chain_A.pdb` |
| **Avogadro** | File → Open |
| **RCSB Mol* Viewer** | Upload at [https://molstar.org/viewer/](https://molstar.org/viewer/) |
| **NGL Viewer** | Online at [https://nglviewer.org/ngl/](https://nglviewer.org/ngl/) |

### Interactive coarse-grained viewer

The pipeline exports a self-contained HTML file (`coarse_grained_topology_viewer.html`) that can be opened directly in any web browser — no installation or server required. Simply double-click the file, or drag it into an open browser window.

The viewer shows each monomer as a bead coloured by type (H = red, G = green, S = blue) and each inter-unit linkage as an edge. Backbone linkages (β-O-4, α-O-4) are drawn as solid lines; C–C cross-links (5-5, β-5, β–β, 4-O-5, β-1) as dashed lines. Bead positions are the centres of mass of the corresponding monomer atoms in the atomistic 3-D structure. Hovering over beads displays metadata (chain, bead ID, type).

### Graph visualisation (Python)

The `lignoforge.core.utils` module provides matplotlib-based graph drawing functions that can be called directly on any `Polymer` object:

```python
from lignoforge.core.utils import draw_big_graph, draw_graph, draw_atomic_graph

polymer = results.polymers[0]   # first chain in the population

# Coarse-grained monomer graph (fast, good for overview)
draw_big_graph(polymer.bigG)

# Full atom-level connectivity graph
draw_graph(polymer.G)

# Atom-level graph with element labels
draw_atomic_graph(polymer.G)
```

### SMILES and 2-D structure rendering

```python
from lignoforge.core.utils import graph_to_smile
from rdkit import Chem
from rdkit.Chem import Draw

polymer = results.polymers[0]
smiles  = graph_to_smile(polymer.G)
mol     = Chem.MolFromSmiles(smiles)
Draw.MolToFile(mol, "chain_A_2d.png", size=(1200, 400))
```

Online SMILES viewers: [ChemDraw Online](https://chemdrawdirect.perkinelmer.cloud/), [Mol2000](https://www.molinspiration.com/cgi-bin/properties).

### Converting to other molecular formats

```python
from lignoforge.core.utils import graph_to_mol

mol = graph_to_mol(
    polymer.G,
    generate_3d=True,
    add_hs_for_3d=True,
    optimize_3d=True,
)

# Write as SDF for use with GROMACS, LAMMPS, or OpenMM
from rdkit.Chem import AllChem
writer = AllChem.SDWriter("chain_A.sdf")
writer.write(mol)
writer.close()
```

---

## Advanced Usage

### Controlling the simulation

```python
results = pipeline.run(
    random_seed=42,
    library_name="pine_kraft",
    simulation_overrides={
        "n_population":    50,
        "Tmetro":          8.0,
        "Tmetro_out":      15.0,
        "i_max":           3000,
        "i_max_ring":      1000,
        "branching_propensity": 0.03,
    }
)
```

### Export options

Export options are passed inside `simulation_overrides` (the pipeline
extracts them automatically before forwarding the rest to the simulator):

```python
results = pipeline.run(
    simulation_overrides={
        "generate_3d":       True,
        "include_hydrogens": True,
        "optimize_3d":       True,
        "max_uff_iterations": 300,   # MMFF geometry optimisation steps
        "n_workers":         None,   # None = all CPU cores; int to cap
    }
)
```

Disable 3-D generation for faster runs (topology and SMILES are still produced):

```python
results = pipeline.run(
    simulation_overrides={"generate_3d": False}
)
```

### Using the low-level API

```python
from lignoforge.core.monomer import monomer_G, monomer_S
from lignoforge.core.polymer import Polymer

# Build a chain manually
p = Polymer(monomer_G())
p.add_specific_monomer(monomer_S(), "beta-O-4")
p.add_specific_monomer(monomer_G(), "5-5")
p.add_random_monomer()
p.add_specific_ring(0, 3)   # ring closure between monomers 0 and 3
```

---

## Input JSON Schema

The input file is validated against the JSON Schema bundled at `lignoforge/io/lignin_info_schema.json`. The minimal valid input is:

```json
{
  "material_origin": {
    "biomass_type": "hardwood"
  },
  "extraction_process": {
    "process_type": "kraft"
  }
}
```

Supported `biomass_type` values: `"hardwood"`, `"softwood"`, `"grass"`, `"miscanthus"`  
Supported `process_type` values: `"kraft"`, `"organosolv"`, `"soda"`, `"lignosulfonate"`, `"steam_explosion"`, `"enzymatic"`, `"milled_wood"`

See `examples/lignin_input_example.json` for a complete annotated example and `lignoforge/io/lignin_info_schema.json` for the full schema.

---

## Documentation

Full documentation, including API reference, concept explanations, and user guides, is available at:

**[https://lignoforge.readthedocs.io](https://lignoforge.readthedocs.io)**

To build the documentation locally:

```bash
pip install -r docs/requirements.txt
cd docs
make html
# Open docs/build/html/index.html in a browser
```

---

## Project Structure

```
lignoforge/
├── core/               # Monomer, PolymerGraph, Polymer, Characterize, rules, utils
├── simulation/         # Trajectory (inner loop) + Simulation (outer population loop)
├── pipeline/           # LigninPipeline, LigninResults, ParameterTranslator
├── priors/             # LigninPriorEstimator + literature data tables
├── io/                 # InputSchemaValidator, LigninExporter, lignin_info_schema.json
├── structure/          # MolecularStructureGenerator, PDBStructureWriter
└── cli/                # Command-line entry points (lignoforge-chain)

examples/
├── demo_run.py                  # End-to-end demo script (population pipeline)
└── lignin_input_example.json    # Minimal hardwood kraft input example

lignin_ff/                       # OPLS-AA force-field files for GROMACS
├── lignin.rtp                   # GROMACS residue topology (9 residue definitions)
├── residuetypes_lignin.dat      # pdb2gmx residuetype entries
├── notes/
│   └── atom_type_assignment.md  # per-atom type justification and charge accounting
└── tools/
    ├── assign_chain_types.py    # assigns OPLS-AA types to a LignoForge PDB
    └── test_assign.py           # validation tests

docs/                   # Sphinx documentation source (Diátaxis structure)
```

---

## Force-Field Parameters (GROMACS / OPLS-AA)

LignoForge ships OPLS-AA force-field parameters for all-atom MD simulations
with GROMACS in the `lignin_ff/` directory.

### Residue definitions

Nine residue topologies are defined in `lignin_ff/lignin.rtp`:

| Group | Names | Description |
|-------|-------|-------------|
| β-O-4 chain units | `GYU`, `HPU`, `SYU` | Internal polymer residues; sp3 side chain; two inter-residue bonds via O4H↔CB |
| Isolated monolignols | `GYM`, `HPM`, `SYM` | Free coniferyl / p-coumaryl / sinapyl alcohol; vinyl side chain; free phenol |
| Neutral saturated refs | `GNM`, `HNM`, `SNM` | Dihydro- analogues for parametrization checks |

All residues carry **net charge 0.000 e** (verified per charge group).

### Automatic type assignment

The script `lignin_ff/tools/assign_chain_types.py` reads any LignoForge PDB,
detects inter-residue bonds, and writes a per-chain custom RTP with correct
OPLS-AA types at all linkage sites:

```python
from lignoforge.core.monomer import Monomer
from lignoforge.core.polymer import Polymer
from lignoforge.structure.pdb import PDBStructureWriter
import subprocess, sys

# 1. Generate PDB with LignoForge
m = Monomer("G", monomer_index=0); m.create()
p = Polymer(m, verbose=False)
p.add_specific_monomer("G", "beta-O-4")
p.add_specific_monomer("G", "beta-O-4")
PDBStructureWriter().write_polymer_pdb(p, "chain.pdb", optimize_3d=True)

# 2. Assign OPLS-AA types → chain_custom.rtp + chain_renamed.pdb
subprocess.run([sys.executable,
    "lignin_ff/tools/assign_chain_types.py", "chain.pdb"], check=True)

# 3. Run GROMACS pdb2gmx
# gmx pdb2gmx -f chain_renamed.pdb -ff ./oplsaa -water tip3p \
#             -o processed.gro -p topol.top -rtpres chain_custom.rtp
```

Or directly from the command line:

```bash
python lignin_ff/tools/assign_chain_types.py chain.pdb
```

See `lignin_ff/README.md` for the full GROMACS workflow.

---

## Citation

If you use LignoForge in academic work, please cite:

```bibtex
@software{lignoforge2026,
  title   = {{LignoForge}: top-down stochastic generation of lignin structural models},
  author  = {DAIMON Team},
  year    = {2026},
  url     = {https://github.com/daimoners/lignoforge},
  version = {0.2.0}
}
```

---

## License

LignoForge is released under the **MIT License**. See [LICENSE](LICENSE) for details.
