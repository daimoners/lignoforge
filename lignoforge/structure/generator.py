"""
Molecular structure generation utilities for lignin polymers.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from typing import Iterable, List, Dict

from lignoforge.core.utils import graph_to_mol
from lignoforge.core.rules import MONOMER_RESIDUE_CODE


def _atomistic_topology_worker(args: tuple) -> dict:
    """
    Module-level worker for ProcessPoolExecutor.

    Parameters arrive as a single tuple so that ``executor.map`` can be used
    with a single iterable.  Each worker creates its own
    ``MolecularStructureGenerator`` instance (the class is stateless).
    """
    polymer, chain_index, random_seed, include_hydrogens, optimize_3d, max_uff_iterations = args
    gen = MolecularStructureGenerator()
    return gen.polymer_atomistic_topology(
        polymer,
        chain_index=chain_index,
        random_seed=random_seed,
        include_hydrogens=include_hydrogens,
        optimize_3d=optimize_3d,
        max_uff_iterations=max_uff_iterations,
    )


# ── RTP-compatible atom name helpers ─────────────────────────────────────────

def _ome_ring_side(o_node: int, G) -> int:
    """Return the 1-based index of the aromatic ring C bonded to an OCH3 oxygen node."""
    for nb in G.neighbors(o_node):
        if G.nodes[nb].get("aromatic"):
            return int(G.nodes[nb]["index"])
    return 0


def _rtp_heavy_name(node: int, G) -> str:
    """Map a heavy-atom graph node to its GROMACS/RTP atom name."""
    attrs = G.nodes[node]
    elem = attrs["element"]
    idx  = int(attrs.get("index", 0))
    grp  = attrs.get("group")
    aro  = bool(attrs.get("aromatic", False))

    if aro:
        return f"C{idx}"
    if grp == "4OH":
        return "O4H"
    if grp == "9OH":
        return "OG"
    if grp == "alpha_OH":
        return "OA"
    if grp == "OCH3":
        if elem == "O":
            side = _ome_ring_side(node, G)
        else:
            o_nb = next(
                nb for nb in G.neighbors(node)
                if G.nodes[nb].get("group") == "OCH3"
                and G.nodes[nb]["element"] == "O"
            )
            side = _ome_ring_side(o_nb, G)
        suffix = "3" if side == 3 else "5"
        return f"OM{suffix}" if elem == "O" else f"CM{suffix}"
    # Non-aromatic side chain (index 7 = Cα, 8 = Cβ, 9 = Cγ)
    sc = {7: "CA", 8: "CB", 9: "CG"}
    return sc.get(idx, f"{elem}{idx}")


def _rtp_h_name(parent_rtp: str, h_index: int, h_total: int) -> str:
    """Return the GROMACS/RTP name for the h_index-th (1-based) H on *parent_rtp*."""
    if parent_rtp == "C2":  return "H2"
    if parent_rtp == "C3":  return "H3"
    if parent_rtp == "C5":  return "H5"
    if parent_rtp == "C6":  return "H6"
    if parent_rtp == "CA":  return "HA"
    if parent_rtp == "CB":  return "HB" if h_total == 1 else f"HB{h_index}"
    if parent_rtp == "CG":  return f"HG{h_index}"
    if parent_rtp == "O4H": return "HO4"
    if parent_rtp == "OG":  return "HOG"
    if parent_rtp == "OA":  return "HOA"
    if parent_rtp == "CM3": return f"HM3{h_index}"
    if parent_rtp == "CM5": return f"HM5{h_index}"
    return f"H{h_index}"


class MolecularStructureGenerator:
    """
    Convert generated polymer graphs to RDKit molecules.
    """

    @staticmethod
    def _chain_id_from_index(chain_index: int) -> str:
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        if chain_index < len(alphabet):
            return alphabet[chain_index]
        major = (chain_index // len(alphabet)) - 1
        minor = chain_index % len(alphabet)
        return alphabet[major % len(alphabet)] + alphabet[minor]

    @staticmethod
    def _atomic_mass(element: str) -> float:
        """Approximate atomic mass (a.m.u.) for the elements found in lignin."""
        return {"C": 12.011, "H": 1.008, "O": 15.999}.get(
            element.strip().capitalize(), 12.011
        )

    def polymer_cg_from_atomistic(self, atomistic_chain: dict) -> dict:
        """
        Derive a coarse-grained topology from an atomistic chain topology dict.

        Each CG bead is placed at the centre of mass of all atoms that belong
        to that monomer (weighted by atomic mass).  CG bonds are taken directly
        from the inter-monomer bonds in the atomistic topology, so chain
        connectivity, branch points, and ring closures (e.g. beta-beta / 5-5)
        are all faithfully preserved.

        Parameters
        ----------
        atomistic_chain : dict
            Single-chain atomistic topology as returned by
            ``polymer_atomistic_topology()``.

        Returns
        -------
        dict
            Single-chain CG topology (JSON schema
            ``lignoforge-coarse-grained-v1``).
        """
        chain_id    = atomistic_chain.get("chain_id", "A")
        chain_index = atomistic_chain.get("chain_index", 0)

        # ── centre-of-mass per monomer ────────────────────────────────────────
        monomer_com:  Dict[int, List[float]] = {}
        monomer_type: Dict[int, str]         = {}

        for monomer in atomistic_chain.get("monomers", []):
            mi    = int(monomer["monomer_id"])
            mtype = monomer.get("residue_name") or monomer.get("monomer_type") or "UNK"
            monomer_type[mi] = mtype

            total_mass = 0.0
            cx = cy = cz = 0.0
            for atom in monomer.get("atoms", []):
                m   = self._atomic_mass(atom.get("element", "C"))
                cx += float(atom["x"]) * m
                cy += float(atom["y"]) * m
                cz += float(atom["z"]) * m
                total_mass += m
            if total_mass > 0:
                cx /= total_mass
                cy /= total_mass
                cz /= total_mass
            monomer_com[mi] = [cx, cy, cz]

        beads = [
            {
                "bead_id":   mi,
                "chain_id":  chain_id,
                "bead_type": monomer_type.get(mi, "UNK"),
                "x":         monomer_com[mi][0],
                "y":         monomer_com[mi][1],
                "z":         monomer_com[mi][2],
            }
            for mi in sorted(monomer_com.keys())
        ]

        # ── inter-monomer bonds → CG links ────────────────────────────────────
        seen:  set        = set()
        links: List[dict] = []
        for bond in atomistic_chain.get("bonds", []):
            if bond.get("scope") != "inter_monomer":
                continue
            mi1 = int(bond["monomer1_id"])
            mi2 = int(bond["monomer2_id"])
            key = (min(mi1, mi2), max(mi1, mi2))
            if key in seen:
                continue
            seen.add(key)
            links.append({
                "chain_id":       chain_id,
                "source_bead_id": mi1,
                "target_bead_id": mi2,
                "linkage_type":   bond.get("linkage_type") or "unknown",
            })

        return {
            "schema":      "lignoforge-coarse-grained-v1",
            "chain_id":    chain_id,
            "chain_index": chain_index,
            "n_beads":     len(beads),
            "beads":       beads,
            "links":       links,
        }

    def population_cg_from_atomistic(self, atomistic_population: dict) -> dict:
        """
        Derive CG population topology from an atomistic population topology dict.

        Parameters
        ----------
        atomistic_population : dict
            Population-level atomistic topology as returned by
            ``population_atomistic_topology()``.

        Returns
        -------
        dict
            Population-level CG topology (schema
            ``lignoforge-coarse-grained-population-v1``).
        """
        chains_at = atomistic_population.get("chains", [atomistic_population])
        chains_cg = []
        for i, chain in enumerate(chains_at):
            print(
                f"            Deriving CG from atomistic: "
                f"chain {i + 1}/{len(chains_at)}"
            )
            chains_cg.append(self.polymer_cg_from_atomistic(chain))
        return {
            "schema":   "lignoforge-coarse-grained-population-v1",
            "n_chains": len(chains_cg),
            "chains":   chains_cg,
        }

    def polymer_atomistic_topology(
        self,
        polymer,
        chain_index: int = 0,
        random_seed: int = 42,
        include_hydrogens: bool = True,
        optimize_3d: bool = True,
        max_uff_iterations: int = 300,
    ) -> dict:
        chain_id = self._chain_id_from_index(chain_index)
        G = polymer.G
        bigG = polymer.bigG
        mol, idx_to_node = graph_to_mol(
            G,
            generate_3d=True,
            random_seed=random_seed,
            optimize_3d=optimize_3d,
            max_uff_iterations=max_uff_iterations,
            add_hs_for_3d=include_hydrogens,
            return_index_map=True,
            chain_label=f"[chain {chain_id}] ",
        )
        if mol is None:
            raise ValueError("Unable to generate RDKit molecule for atomistic topology")

        conf = mol.GetConformer()
        idx_to_node = idx_to_node or {}

        monomer_records: Dict[int, dict] = {}
        atom_to_monomer: Dict[int, int] = {}

        # Pre-pass: for each heavy-atom graph node, collect the RDKit indices of
        # its explicit H neighbours, so we can distinguish HB (single) from
        # HB1/HB2 (two hydrogens) etc.
        h_per_heavy: Dict[int, List[int]] = {}
        for _a in mol.GetAtoms():
            if _a.GetSymbol() != "H":
                continue
            _parent_rd = next(
                (_nb.GetIdx() for _nb in _a.GetNeighbors() if _nb.GetIdx() in idx_to_node),
                None,
            )
            if _parent_rd is not None:
                _hn = idx_to_node[_parent_rd]
                h_per_heavy.setdefault(_hn, []).append(_a.GetIdx())

        # Per-node hydrogen sequence counter (reset for every node): tracks which
        # H on a given heavy atom we are currently naming (1-based).
        h_seq_counter: Dict[int, int] = {}

        for atom in mol.GetAtoms():
            rd_idx = int(atom.GetIdx())
            atom_id = rd_idx + 1

            heavy_node = idx_to_node.get(rd_idx)
            if heavy_node is not None:
                node_attrs = G.nodes[heavy_node]
                mi = int(node_attrs["mi"])
                mtype = bigG.nodes[mi].get("mtype", "UNK") if mi in bigG.nodes else "UNK"
                chem_index = int(node_attrs.get("index", heavy_node))
                group = node_attrs.get("group")
                bonding = bool(node_attrs.get("bonding", False))
            else:
                neigh = [n.GetIdx() for n in atom.GetNeighbors()]
                heavy_neighbor = next((n for n in neigh if n in idx_to_node), None)
                if heavy_neighbor is None:
                    mi = 0
                    mtype = "UNK"
                else:
                    heavy_node = idx_to_node[heavy_neighbor]
                    node_attrs = G.nodes[heavy_node]
                    mi = int(node_attrs["mi"])
                    mtype = bigG.nodes[mi].get("mtype", "UNK") if mi in bigG.nodes else "UNK"
                chem_index = 0
                group = "explicit_H"
                bonding = False

            atom_to_monomer[atom_id] = mi

            # ── RTP-compatible atom name ──────────────────────────────────────
            if group == "explicit_H":
                if heavy_node is not None:
                    _pname  = _rtp_heavy_name(heavy_node, G)
                    _htotal = len(h_per_heavy.get(heavy_node, []))
                    h_seq_counter[heavy_node] = h_seq_counter.get(heavy_node, 0) + 1
                    atom_name = _rtp_h_name(_pname, h_seq_counter[heavy_node], _htotal)
                else:
                    atom_name = "H"
            else:
                atom_name = _rtp_heavy_name(heavy_node, G)

            if mi not in monomer_records:
                residue_name = MONOMER_RESIDUE_CODE.get(str(mtype), str(mtype)[:3])
                monomer_records[mi] = {
                    "monomer_id": int(mi),
                    "chain_id": chain_id,
                    "monomer_type": residue_name,
                    "residue_name": residue_name,
                    "n_atoms": 0,
                    "atoms": [],
                }

            pos = conf.GetAtomPosition(rd_idx)
            monomer_records[mi]["atoms"].append(
                {
                    "atom_id": atom_id,
                    "chain_id": chain_id,
                    "monomer_id": int(mi),
                    "residue_name": monomer_records[mi]["residue_name"],
                    "atom_name": atom_name,
                    "element": atom.GetSymbol(),
                    "chemical_index": chem_index,
                    "group": group,
                    "bonding": bonding,
                    "x": float(pos.x),
                    "y": float(pos.y),
                    "z": float(pos.z),
                }
            )

        monomers = [monomer_records[mi] for mi in sorted(monomer_records.keys())]
        for monomer in monomers:
            monomer["atoms"] = sorted(monomer["atoms"], key=lambda a: a["atom_id"])
            monomer["n_atoms"] = len(monomer["atoms"])

        edge_lookup = {}
        for n1, n2, attrs in G.edges(data=True):
            edge_lookup[tuple(sorted((int(n1), int(n2))))] = attrs

        bonds = []
        for bond in mol.GetBonds():
            a1 = int(bond.GetBeginAtomIdx()) + 1
            a2 = int(bond.GetEndAtomIdx()) + 1
            mi1 = atom_to_monomer.get(a1, 0)
            mi2 = atom_to_monomer.get(a2, 0)

            linkage_type = None
            if bond.GetBeginAtomIdx() in idx_to_node and bond.GetEndAtomIdx() in idx_to_node:
                n1 = idx_to_node[bond.GetBeginAtomIdx()]
                n2 = idx_to_node[bond.GetEndAtomIdx()]
                eattrs = edge_lookup.get(tuple(sorted((int(n1), int(n2)))))
                if eattrs is not None:
                    linkage_type = eattrs.get("btype")

            bonds.append(
                {
                    "chain_id": chain_id,
                    "atom1_id": a1,
                    "atom2_id": a2,
                    "bond_order": int(round(float(bond.GetBondTypeAsDouble()))),
                    "linkage_type": linkage_type,
                    "scope": "intra_monomer" if mi1 == mi2 else "inter_monomer",
                    "monomer1_id": int(mi1),
                    "monomer2_id": int(mi2),
                }
            )

        return {
            "schema": "lignoforge-atomistic-topology-v1",
            "chain_id": chain_id,
            "chain_index": int(chain_index),
            "n_monomers": len(monomers),
            "n_atoms": int(mol.GetNumAtoms()),
            "n_bonds": int(mol.GetNumBonds()),
            "monomers": monomers,
            "bonds": bonds,
        }

    def population_atomistic_topology(
        self,
        polymers: Iterable,
        random_seed: int = 42,
        include_hydrogens: bool = True,
        optimize_3d: bool = True,
        max_uff_iterations: int = 300,
        n_workers: int | None = None,
    ) -> dict:
        """
        Generate atomistic topologies for a population of polymers.

        Parameters
        ----------
        n_workers : int or None
            Number of parallel worker processes to use for 3-D embedding.
            ``None`` (default) uses all logical CPU cores.
            Set to ``1`` to disable parallelism (sequential, useful for
            debugging or when the population is small).
        """
        polymers_list = list(polymers)
        n = len(polymers_list)

        args_list = [
            (p, i, random_seed + i, include_hydrogens, optimize_3d, max_uff_iterations)
            for i, p in enumerate(polymers_list)
        ]

        if n_workers is None:
            n_workers = min(n, os.cpu_count() or 1)

        if n_workers <= 1 or n <= 1:
            chains = []
            for i, args in enumerate(args_list):
                print(f"            Generating atomistic topology: chain {i + 1}/{n}")
                chains.append(_atomistic_topology_worker(args))
        else:
            print(
                f"            Generating atomistic topology: "
                f"{n} chains on {n_workers} parallel workers..."
            )
            with ProcessPoolExecutor(max_workers=n_workers) as executor:
                chains = list(executor.map(_atomistic_topology_worker, args_list))
            print(f"            All {n} chains completed.")

        return {
            "schema": "lignoforge-atomistic-population-v1",
            "n_chains": len(chains),
            "chains": chains,
        }

    def polymer_to_mol(
        self,
        polymer,
        generate_3d: bool = True,
        random_seed: int = 42,
        optimize_3d: bool = True,
        max_uff_iterations: int = 200,
        explicit_hydrogens: bool = True,
    ):
        return graph_to_mol(
            polymer.G,
            generate_3d=generate_3d,
            random_seed=random_seed,
            optimize_3d=optimize_3d,
            max_uff_iterations=max_uff_iterations,
            add_hs_for_3d=explicit_hydrogens,
        )

    def population_to_mols(
        self,
        polymers: Iterable,
        generate_3d: bool = True,
        random_seed: int = 42,
        optimize_3d: bool = True,
        max_uff_iterations: int = 200,
        explicit_hydrogens: bool = True,
    ) -> List:
        mols = []
        for i, polymer in enumerate(polymers):
            mol = self.polymer_to_mol(
                polymer,
                generate_3d=generate_3d,
                random_seed=random_seed + i,
                optimize_3d=optimize_3d,
                max_uff_iterations=max_uff_iterations,
                explicit_hydrogens=explicit_hydrogens,
            )
            if mol is not None:
                mols.append(mol)
        return mols
