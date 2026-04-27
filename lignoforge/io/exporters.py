"""
Export utilities for generated lignin populations.
"""

from __future__ import annotations

import json
import os
from typing import Iterable

import numpy as np

from lignoforge.core.utils import graph_to_smile, graph_to_mol
from lignoforge.core.characterization import Characterize
from lignoforge.structure.generator import MolecularStructureGenerator
from lignoforge.structure.pdb import PDBStructureWriter
from lignoforge.structure.visualization import write_cg_topology_viewer_html


class LigninExporter:
    """
    Export polymer populations to common formats.
    """

    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        self.pdb_writer = PDBStructureWriter()
        self.topology_generator = MolecularStructureGenerator()

    def export_smiles(self, polymers: Iterable, filename: str = "population.smi") -> str:
        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            for p in polymers:
                smi = graph_to_smile(p.G)
                if smi:
                    f.write(smi + "\n")
        return path

    def export_sdf(self, polymers: Iterable, folder: str = "sdf") -> str:
        out = os.path.join(self.output_dir, folder)
        os.makedirs(out, exist_ok=True)
        for i, p in enumerate(polymers):
            mol = graph_to_mol(p.G)
            if mol is not None:
                mol_path = os.path.join(out, f"polymer_{i}.mol")
                from rdkit import Chem
                Chem.MolToMolFile(mol, mol_path)
        return out

    def export_summary(self, priors: dict, sim_kwargs: dict, filename: str = "summary.json") -> str:
        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            json.dump({"priors": priors, "simulation": sim_kwargs}, f, indent=2)
        return path

    def export_input(self, input_data: dict, filename: str = "input_high_level.json") -> str:
        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            json.dump(input_data, f, indent=2)
        return path

    def export_priors(self, priors: dict, filename: str = "estimated_priors.json") -> str:
        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            json.dump(priors, f, indent=2)
        return path

    def export_simulation_parameters(
        self,
        sim_kwargs: dict,
        filename: str = "simulation_parameters.json",
    ) -> str:
        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            json.dump(sim_kwargs, f, indent=2)
        return path

    def export_chain_statistics(
        self,
        polymers: Iterable,
        filename: str = "chain_statistics.json",
    ) -> str:
        rows = []
        for i, polymer in enumerate(polymers):
            ch = Characterize(polymer)
            row = ch.summary()
            row["polymer_id"] = i
            rows.append(row)

        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            json.dump(rows, f, indent=2)
        return path

    def export_population_statistics(
        self,
        polymers: Iterable,
        filename: str = "population_statistics.json",
    ) -> str:
        rows = []
        for polymer in polymers:
            ch = Characterize(polymer)
            rows.append(ch.summary())

        if not rows:
            stats = {
                "n_polymers": 0,
                "numeric": {},
            }
        else:
            numeric_keys = [
                key for key, val in rows[0].items()
                if isinstance(val, (int, float, np.integer, np.floating))
            ]
            numeric_stats = {}
            for key in numeric_keys:
                arr = np.array([float(r[key]) for r in rows], dtype=float)
                numeric_stats[key] = {
                    "mean": float(np.mean(arr)),
                    "std": float(np.std(arr)),
                    "min": float(np.min(arr)),
                    "max": float(np.max(arr)),
                }
            stats = {
                "n_polymers": len(rows),
                "numeric": numeric_stats,
            }

        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            json.dump(stats, f, indent=2)
        return path

    def export_coarse_grained_topology_json(
        self,
        atomistic_population: dict,
        filename: str = "coarse_grained_topology.json",
    ) -> str:
        """
        Derive a CG topology JSON from an already-computed atomistic population
        topology dict (bead positions = centres of mass).

        .. deprecated::
            Prefer :meth:`export_coarse_grained_from_atomistic` which also
            produces the HTML viewer and CG PDB files in a single call.
        """
        cg_topology = self.topology_generator.population_cg_from_atomistic(
            atomistic_population
        )
        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            json.dump(cg_topology, f, indent=2)
        return path

    def export_atomistic_topology_json(
        self,
        polymers: Iterable,
        filename: str = "atomistic_topology.json",
        random_seed: int = 42,
        include_hydrogens: bool = True,
        optimize_3d: bool = True,
        max_uff_iterations: int = 500,
        n_workers: int | None = None,
        return_topology: bool = False,
    ):
        topology = self.topology_generator.population_atomistic_topology(
            polymers,
            random_seed=random_seed,
            include_hydrogens=include_hydrogens,
            optimize_3d=optimize_3d,
            max_uff_iterations=max_uff_iterations,
            n_workers=n_workers,
        )
        path = os.path.join(self.output_dir, filename)
        with open(path, "w") as f:
            json.dump(topology, f, indent=2)
        if return_topology:
            return path, topology
        return path

    def export_coarse_grained_from_atomistic(
        self,
        atomistic_population: dict,
        json_filename: str = "coarse_grained_topology.json",
        html_filename: str = "coarse_grained_topology_viewer.html",
        pdb_folder: str = "pdb_cg",
        pdb_basename: str = "polymer_cg",
    ) -> dict:
        """
        Derive the complete CG layer from an atomistic population topology.

        Each monomer bead is placed at the centre of mass of its atoms;
        inter-monomer bonds (chain, branch, ring-closing) are taken directly
        from the atomistic bond list.

        Outputs
        -------
        Writes three artefacts to ``self.output_dir``:

        - **JSON** – structural / topological representation.
        - **HTML** – interactive 3D Plotly viewer (self-contained, open in browser).
        - **PDB** – pseudo-atom PDB (one bead per monomer, VMD-compatible);
          each bead is placed at the monomer centre-of-mass.

        Returns
        -------
        dict with keys ``cg_json``, ``cg_html``, ``cg_pdb``.
        """
        print("          [CG] Deriving beads from atomistic centres-of-mass...")
        cg_topology = self.topology_generator.population_cg_from_atomistic(
            atomistic_population
        )

        # JSON
        json_path = os.path.join(self.output_dir, json_filename)
        with open(json_path, "w") as f:
            json.dump(cg_topology, f, indent=2)
        print(f"          [CG] JSON  → {json_path}")

        # HTML interactive viewer
        html_path = os.path.join(self.output_dir, html_filename)
        try:
            html_path = write_cg_topology_viewer_html(json_path, html_path)
            print(f"          [CG] HTML  → {html_path}")
        except ImportError:
            html_path = None

        # PDB (VMD-compatible pseudo-atoms)
        pdb_dir      = os.path.join(self.output_dir, pdb_folder)
        cg_pdb_paths = self.pdb_writer.write_cg_population_pdbs_from_cg_population(
            cg_topology, pdb_dir, basename=pdb_basename
        )
        print(f"          [CG] PDB   → {pdb_dir}/ ({len(cg_pdb_paths)} file(s))")

        return {
            "cg_json": json_path,
            "cg_html": html_path,
            "cg_pdb":  cg_pdb_paths,
        }

    def export_full_pipeline_bundle(
        self,
        input_data: dict,
        priors: dict,
        sim_kwargs: dict,
        polymers: Iterable,
        export_options: dict | None = None,
    ) -> dict:
        # Convert to list to allow multiple iterations
        polymers_list = list(polymers)
        n_chains = len(polymers_list)
        print(f"          Exporting {n_chains} chain(s)...")

        print(f"          [Exporting] Input data...")
        input_path = self.export_input(input_data)

        print(f"          [Exporting] Estimated priors...")
        priors_path = self.export_priors(priors)

        print(f"          [Exporting] Simulation parameters...")
        sim_path = self.export_simulation_parameters(sim_kwargs)

        print(f"          [Exporting] Summary...")
        summary_path = self.export_summary(priors, sim_kwargs)

        print(f"          [Exporting] SMILES strings...")
        smiles_path = self.export_smiles(polymers_list)

        options             = export_options or {}
        generate_3d         = bool(options.get("generate_3d", True))
        include_hydrogens   = bool(options.get("include_hydrogens", True))
        optimize_3d         = bool(options.get("optimize_3d", True))
        max_uff_iterations  = int(options.get("max_uff_iterations", 300))
        n_workers           = options.get("n_workers", None)   # None → all CPU cores
        if n_workers is not None:
            n_workers = int(n_workers)

        atomistic_json_path  = None
        atomistic_topology   = None
        cg_artifacts: dict   = {"cg_json": None, "cg_html": None, "cg_pdb": []}
        pdb_paths:    list   = []

        if generate_3d:
            print(f"          [Exporting] Atomistic topology (3D coordinates)...")
            atomistic_json_path, atomistic_topology = self.export_atomistic_topology_json(
                polymers_list,
                include_hydrogens=include_hydrogens,
                optimize_3d=optimize_3d,
                max_uff_iterations=max_uff_iterations,
                n_workers=n_workers,
                return_topology=True,
            )

            print(f"          [Exporting] Coarse-grained topology (JSON + viewer + PDB)...")
            cg_artifacts = self.export_coarse_grained_from_atomistic(atomistic_topology)

            print(f"          [Exporting] Atomistic PDB files...")
            pdb_paths = self.pdb_writer.write_population_pdbs_from_atomistic_population(
                atomistic_topology,
                os.path.join(self.output_dir, "pdb"),
                basename="polymer",
            )
        else:
            print(f"          [Exporting] 3D export skipped (generate_3d=False).")

        print(f"          [Exporting] Chain statistics...")
        chains_path = self.export_chain_statistics(polymers_list)

        print(f"          [Exporting] Population statistics...")
        population_path = self.export_population_statistics(polymers_list)

        return {
            "input":                     input_path,
            "priors":                    priors_path,
            "simulation":               sim_path,
            "summary":                  summary_path,
            "smiles":                   smiles_path,
            "atomistic_topology":       atomistic_json_path,
            "coarse_grained_topology":  cg_artifacts["cg_json"],
            "coarse_grained_viewer":    cg_artifacts["cg_html"],
            "coarse_grained_pdb":       cg_artifacts["cg_pdb"],
            "chains":                   chains_path,
            "population":               population_path,
            "pdb":                      pdb_paths,
        }

    def export_single_pdb(
        self,
        polymer,
        filename: str = "polymer.pdb",
        random_seed: int = 42,
        optimize_3d: bool = True,
        max_uff_iterations: int = 200,
        explicit_hydrogens: bool = True,
    ) -> str:
        path = os.path.join(self.output_dir, filename)
        return self.pdb_writer.write_polymer_pdb(
            polymer,
            path,
            random_seed=random_seed,
            optimize_3d=optimize_3d,
            max_uff_iterations=max_uff_iterations,
            explicit_hydrogens=explicit_hydrogens,
        )

    def export_population_pdbs(
        self,
        polymers: Iterable,
        folder: str = "pdb",
        basename: str = "polymer",
        random_seed: int = 42,
        optimize_3d: bool = True,
        max_uff_iterations: int = 200,
        explicit_hydrogens: bool = True,
    ) -> list[str]:
        out = os.path.join(self.output_dir, folder)
        return self.pdb_writer.write_population_pdbs(
            polymers,
            out,
            basename=basename,
            random_seed=random_seed,
            optimize_3d=optimize_3d,
            max_uff_iterations=max_uff_iterations,
            explicit_hydrogens=explicit_hydrogens,
        )
