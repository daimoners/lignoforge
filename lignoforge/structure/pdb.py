"""
PDB export helpers for lignin polymer structures.
"""

from __future__ import annotations

import os
from typing import Iterable, List

from lignoforge.structure.generator import MolecularStructureGenerator


class PDBStructureWriter:
    """
    Generate and export complete PDB structures (with explicit hydrogens).
    """

    def __init__(self):
        self.generator = MolecularStructureGenerator()

    @staticmethod
    def _format_pdb_atom_line(
        serial: int,
        atom_name: str,
        residue_name: str,
        chain_id: str,
        residue_seq: int,
        x: float,
        y: float,
        z: float,
        element: str,
    ) -> str:
        atom_name = (atom_name or "X")[:4]
        residue_name = (residue_name or "LIG")[:3]
        chain_id = (chain_id or "A")[:1]
        element = (element or "X")[:2]
        return (
            f"HETATM{serial:5d} {atom_name:<4s} {residue_name:>3s} {chain_id:1s}"
            f"{residue_seq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
        )

    def _write_atomistic_topology_as_pdb(self, chain_topology: dict, output_path: str) -> str:
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

        chain_id = chain_topology["chain_id"]

        atoms = []
        for monomer in chain_topology["monomers"]:
            monomer_id = int(monomer["monomer_id"])
            residue_name = monomer.get("residue_name", monomer.get("monomer_type", "LIG"))
            for atom in monomer["atoms"]:
                atoms.append(
                    {
                        "atom_id": int(atom["atom_id"]),
                        "atom_name": atom.get("atom_name", atom.get("element", "X")),
                        "element": atom.get("element", "X"),
                        "residue_name": residue_name,
                        "residue_seq": monomer_id + 1,
                        "chain_id": chain_id,
                        "x": float(atom["x"]),
                        "y": float(atom["y"]),
                        "z": float(atom["z"]),
                    }
                )

        atoms = sorted(atoms, key=lambda a: (a["residue_seq"], a["atom_id"]))
        serial_map = {a["atom_id"]: i + 1 for i, a in enumerate(atoms)}

        lines = [
            "HEADER    LIGNOFORGE GENERATED STRUCTURE",
            f"TITLE     CHAIN {chain_id} | MONOMERS {chain_topology.get('n_monomers', 0)}",
        ]

        for atom in atoms:
            lines.append(
                self._format_pdb_atom_line(
                    serial=serial_map[atom["atom_id"]],
                    atom_name=atom["atom_name"],
                    residue_name=atom["residue_name"],
                    chain_id=atom["chain_id"],
                    residue_seq=atom["residue_seq"],
                    x=atom["x"],
                    y=atom["y"],
                    z=atom["z"],
                    element=atom["element"],
                )
            )

        seen = set()
        for bond in chain_topology.get("bonds", []):
            a1 = int(bond["atom1_id"])
            a2 = int(bond["atom2_id"])
            if a1 not in serial_map or a2 not in serial_map:
                continue
            key = tuple(sorted((a1, a2)))
            if key in seen:
                continue
            seen.add(key)
            s1, s2 = serial_map[a1], serial_map[a2]
            lines.append(f"CONECT{s1:5d}{s2:5d}")

        lines.append("END")

        with open(output_path, "w") as f:
            f.write("\n".join(lines) + "\n")

        return output_path

    def write_polymer_pdb(
        self,
        polymer,
        output_path: str,
        chain_index: int = 0,
        random_seed: int = 42,
        optimize_3d: bool = False,
        max_uff_iterations: int = 120,
        explicit_hydrogens: bool = True,
    ) -> str:
        topology = self.generator.polymer_atomistic_topology(
            polymer,
            chain_index=chain_index,
            random_seed=random_seed,
            include_hydrogens=explicit_hydrogens,
            optimize_3d=optimize_3d,
            max_uff_iterations=max_uff_iterations,
        )
        return self._write_atomistic_topology_as_pdb(topology, output_path)

    def write_population_pdbs(
        self,
        polymers: Iterable,
        output_dir: str,
        basename: str = "polymer",
        random_seed: int = 42,
        optimize_3d: bool = False,
        max_uff_iterations: int = 120,
        explicit_hydrogens: bool = True,
    ) -> List[str]:
        os.makedirs(output_dir, exist_ok=True)
        polymers_list = list(polymers)
        written = []
        for i, polymer in enumerate(polymers_list):
            print(f"            Writing PDB: chain {i+1}/{len(polymers_list)}")
            path = os.path.join(output_dir, f"{basename}_{i}.pdb")
            self.write_polymer_pdb(
                polymer,
                path,
                chain_index=i,
                random_seed=random_seed + i,
                optimize_3d=optimize_3d,
                max_uff_iterations=max_uff_iterations,
                explicit_hydrogens=explicit_hydrogens,
            )
            written.append(path)
        return written

    def write_population_pdbs_from_atomistic_population(
        self,
        atomistic_population: dict,
        output_dir: str,
        basename: str = "polymer",
    ) -> List[str]:
        os.makedirs(output_dir, exist_ok=True)
        chains = atomistic_population.get("chains", [])
        written = []

        for i, chain_topology in enumerate(chains):
            print(f"            Writing PDB: chain {i+1}/{len(chains)}")
            path = os.path.join(output_dir, f"{basename}_{i}.pdb")
            self._write_atomistic_topology_as_pdb(chain_topology, path)
            written.append(path)

        return written

    # ── Coarse-grained PDB (VMD-compatible pseudo-atoms) ─────────────────────

    @staticmethod
    def _write_cg_topology_as_pdb(cg_chain_topology: dict, output_path: str) -> str:
        """
        Write a CG chain topology as a PDB file suitable for VMD (and other
        molecular viewers).

        Each monomer bead is written as a single HETATM pseudo-atom named
        ``BB``.  CONECT records encode every inter-monomer bond so that chain
        connectivity, branch points, and ring closures all appear correctly.
        The coordinate space is identical to the atomistic PDB (bead positions
        are centres-of-mass of monomer atoms), so both files can be loaded
        together in VMD and displayed simultaneously.
        """
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)

        chain_id = cg_chain_topology.get("chain_id", "A")
        beads    = cg_chain_topology.get("beads", [])
        links    = cg_chain_topology.get("links", [])

        # bead_id → 1-based PDB serial
        serial_map = {int(b["bead_id"]): i + 1 for i, b in enumerate(beads)}

        lines = [
            "HEADER    LIGNOFORGE COARSE-GRAINED STRUCTURE",
            f"TITLE     CHAIN {chain_id} | {len(beads)} BEADS | COM-DERIVED FROM ATOMISTIC",
            "REMARK    Each bead = centre-of-mass of one monomer (H / G / S unit).",
            "REMARK    Atom name: BB  |  Element: C  (pseudo-atom, not a real carbon).",
            "REMARK    CONECT records encode ALL inter-monomer linkages.",
            "REMARK    Load alongside the atomistic PDB in VMD to overlay both models.",
        ]

        for i, bead in enumerate(beads):
            serial   = i + 1
            res_name = (bead.get("bead_type") or "UNK")[:3]
            res_seq  = int(bead["bead_id"]) + 1
            x, y, z  = float(bead["x"]), float(bead["y"]), float(bead["z"])
            lines.append(
                f"HETATM{serial:5d}  BB  {res_name:>3s} {chain_id:1s}{res_seq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
            )

        seen: set = set()
        for link in links:
            bi1 = int(link["source_bead_id"])
            bi2 = int(link["target_bead_id"])
            if bi1 not in serial_map or bi2 not in serial_map:
                continue
            key = tuple(sorted((bi1, bi2)))
            if key in seen:
                continue
            seen.add(key)
            s1, s2 = serial_map[bi1], serial_map[bi2]
            lines.append(f"CONECT{s1:5d}{s2:5d}")

        lines.append("END")

        with open(output_path, "w") as f:
            f.write("\n".join(lines) + "\n")

        return output_path

    def write_cg_population_pdbs_from_cg_population(
        self,
        cg_population: dict,
        output_dir: str,
        basename: str = "polymer_cg",
    ) -> List[str]:
        """Write one CG PDB file per chain in a CG population topology dict."""
        os.makedirs(output_dir, exist_ok=True)
        chains  = cg_population.get("chains", [])
        written = []
        for i, chain in enumerate(chains):
            print(f"            Writing CG PDB: chain {i + 1}/{len(chains)}")
            path = os.path.join(output_dir, f"{basename}_{i}.pdb")
            self._write_cg_topology_as_pdb(chain, path)
            written.append(path)
        return written
