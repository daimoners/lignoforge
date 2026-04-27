#!/usr/bin/env python
"""
lignoforge-chain — single-chain lignin builder
===============================================

Build one or more individual lignin chains from a JSON input file, with
full control over the chain morphology from the command line.  Each chain
is grown by the same kMC/MCMC trajectory engine used by the population
pipeline, but here chains are generated one at a time and exported
immediately.

Usage
-----
    lignoforge-chain input.json [options]

    # or without installation:
    python -m lignoforge.cli.build_chain input.json [options]

Chain size
----------
    --n-monomers N      Target degree of polymerisation (monomer units)
    --mw-target MW      Approximate target MW in g/mol  (overrides --n-monomers)

Monomer composition
-------------------
    --S-fraction F      Syringyl fraction   (0-1; overrides JSON)
    --G-fraction F      Guaiacyl fraction   (0-1; overrides JSON)
    --H-fraction F      p-Hydroxyphenyl     (0-1; overrides JSON)
    --monomer-type T    Force a single type: H | G | S | mix (default: mix)

Linkage distribution
--------------------
    --beta-O-4 F        beta-O-4 fraction (0-1)
    --alpha-O-4 F       alpha-O-4 fraction
    --4-O-5 F           4-O-5 fraction
    --5-5 F             5-5 fraction
    --beta-5 F          beta-5 fraction
    --beta-beta F       beta-beta fraction
    --beta-1 F          beta-1 fraction
    (fractions are re-normalised automatically)

Branching
---------
    --branching F       Branching propensity per MC step (0=linear, default 0)

Simulation
----------
    --seed N            Random seed (default: 42)
    --n-chains N        Number of independent chains (default: 1)
    --Tmetro T          Metropolis temperature (default: from priors)
    --max-steps N       Max MC steps per chain (default: 2000)

Output
------
    --output DIR        Output directory (default: chain_output/)
    --name NAME         Base name for output files (default: chain)
    --format FMT        Comma-separated output formats:
                        smiles, sdf, pdb, pdb-cg, json-atomistic, json-cg, all
                        (default: smiles,pdb)
    --no-3d             Skip 3-D coordinates (fast; SMILES and graph-only JSON)
    --no-H              Strip explicit hydrogens from 3-D output
    --max-iter N        Max MMFF/UFF optimisation iterations (default: 300)
    --verbose           Print per-step chain-growth progress

Examples
--------
    # Hardwood kraft G/S chain, ~15 monomers
    lignoforge-chain input.json --n-monomers 15

    # Softwood pure-G, beta-O-4 rich, ~3000 g/mol
    lignoforge-chain input.json --mw-target 3000 --G-fraction 1.0 \\
        --beta-O-4 0.80 --output run1/ --format pdb

    # Three branched chains, all formats
    lignoforge-chain input.json --n-monomers 20 --branching 0.05 \\
        --n-chains 3 --format all --seed 7
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from pathlib import Path
from typing import Optional

import numpy as np


# ── Built-in defaults (used when no JSON input file is provided) ───────────────

_DEFAULT_INPUT = {
    "material_origin": {"biomass_type": "hardwood"},
    "extraction_process": {"process_type": "kraft"},
}


# ── Input validators ───────────────────────────────────────────────────────────

def _fraction(v: str) -> float:
    x = float(v)
    if not (0.0 <= x <= 1.0):
        raise argparse.ArgumentTypeError(f"must be in [0, 1]; got {x}")
    return x


def _positive_int(v: str) -> int:
    x = int(v)
    if x < 1:
        raise argparse.ArgumentTypeError(f"must be >= 1; got {x}")
    return x


def _positive_float(v: str) -> float:
    x = float(v)
    if x <= 0:
        raise argparse.ArgumentTypeError(f"must be > 0; got {x}")
    return x


# ── Argument parser ────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="lignoforge-chain",
        description=(
            "Build lignin chains with explicit morphology control.\n"
            "Reads a LignoForge JSON input file and allows per-run "
            "overrides of all structural parameters."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument(
        "input",
        metavar="INPUT_JSON",
        nargs="?",
        default=None,
        help="Path to a LignoForge JSON input file (optional; uses built-in defaults if omitted).",
    )

    # ── Size ────────────────────────────────────────────────────────────────────
    sg = p.add_argument_group("Chain size")
    sex = sg.add_mutually_exclusive_group()
    sex.add_argument("--n-monomers",  metavar="N",  type=_positive_int,   default=None,
                     help="Target number of monomeric units.")
    sex.add_argument("--mw-target",   metavar="MW", type=_positive_float, default=None,
                     help="Approximate target MW (g/mol).")

    # ── Composition ─────────────────────────────────────────────────────────────
    cg = p.add_argument_group("Monomer composition")
    cg.add_argument("--S-fraction",   metavar="F", type=_fraction, default=None, dest="S_fraction")
    cg.add_argument("--G-fraction",   metavar="F", type=_fraction, default=None, dest="G_fraction")
    cg.add_argument("--H-fraction",   metavar="F", type=_fraction, default=None, dest="H_fraction")
    cg.add_argument("--monomer-type", metavar="T", default=None,
                    choices=["H", "G", "S", "mix"],
                    help="Force a single monomer type or 'mix' (default: mix).")

    # ── Linkages ─────────────────────────────────────────────────────────────────
    lg = p.add_argument_group("Linkage distribution (fractions are re-normalised)")
    lg.add_argument("--beta-O-4",   metavar="F", type=_fraction, default=None, dest="lk_beta_O_4")
    lg.add_argument("--alpha-O-4",  metavar="F", type=_fraction, default=None, dest="lk_alpha_O_4")
    lg.add_argument("--4-O-5",      metavar="F", type=_fraction, default=None, dest="lk_4_O_5")
    lg.add_argument("--5-5",        metavar="F", type=_fraction, default=None, dest="lk_5_5")
    lg.add_argument("--beta-5",     metavar="F", type=_fraction, default=None, dest="lk_beta_5")
    lg.add_argument("--beta-beta",  metavar="F", type=_fraction, default=None, dest="lk_beta_beta")
    lg.add_argument("--beta-1",     metavar="F", type=_fraction, default=None, dest="lk_beta_1")

    # ── Branching ────────────────────────────────────────────────────────────────
    bg = p.add_argument_group("Branching / topology")
    bg.add_argument("--branching", metavar="F", type=_fraction, default=None,
                    help="Branching propensity per MC step (0=linear).")

    # ── Simulation ──────────────────────────────────────────────────────────────
    simg = p.add_argument_group("Simulation")
    simg.add_argument("--seed",       metavar="N", type=int,            default=42)
    simg.add_argument("--n-chains",   metavar="N", type=_positive_int,  default=1,
                      help="Number of independent chains (default: 1).")
    simg.add_argument("--Tmetro",     metavar="T", type=_positive_float, default=None,
                      help="Metropolis temperature (default: from priors).")
    simg.add_argument("--max-steps",  metavar="N", type=_positive_int,  default=2000,
                      help="Max MC steps per chain (default: 2000).")
    simg.add_argument("--n-workers",  metavar="N", type=int, default=None, dest="n_workers",
                      help="Parallel worker processes for 3-D embedding "
                           "(default: all CPU cores; 1=sequential).")

    # ── Output ──────────────────────────────────────────────────────────────────
    og = p.add_argument_group("Output")
    og.add_argument("--output",   metavar="DIR",  default="chain_output",
                    help="Output directory (default: chain_output/).")
    og.add_argument("--name",     metavar="NAME", default="chain",
                    help="Base name for output files (default: chain).")
    og.add_argument(
        "--format", metavar="FMT", default="smiles,pdb", dest="formats",
        help=(
            "Comma-separated list of output formats (default: smiles,pdb). "
            "Available: smiles, population-smiles, sdf, pdb, pdb-cg, "
            "json-atomistic, json-cg, html, all."
        ),
    )
    og.add_argument("--no-3d",         action="store_true",
                    help="Skip 3-D coordinate generation entirely.")
    og.add_argument("--no-optimize",   action="store_true",
                    help="Embed 3-D coordinates but skip force-field optimisation "
                         "(faster, less accurate geometry).")
    og.add_argument("--no-H",          action="store_true",
                    help="Strip explicit hydrogens from 3-D output.")
    og.add_argument("--max-iter",      metavar="N", type=_positive_int, default=300, dest="max_iter",
                    help="Max MMFF/UFF optimisation iterations (default: 300).")
    og.add_argument("--export-priors", action="store_true",
                    help="Write estimated_priors.json to the output directory.")
    og.add_argument("--export-sim-params", action="store_true",
                    help="Write simulation_parameters.json to the output directory.")
    og.add_argument("--population-stats", action="store_true",
                    help="Write population_statistics.json with ensemble-level stats.")
    og.add_argument("--verbose",       action="store_true",
                    help="Print per-step chain-growth progress.")

    return p


# ── MW → DP estimation ─────────────────────────────────────────────────────────

_MONOMER_MW = {"H": 150.17, "G": 180.20, "S": 210.23}
_BOND_MW_LOSS = 18.015  # approximate loss per inter-monomer condensation bond


def _dp_from_mw(mw_target: float, monomer_dist: list) -> int:
    """Convert MW target to approximate target DP."""
    avg_mw = sum(f * _MONOMER_MW[t] for f, t in zip(monomer_dist, ("H", "G", "S")))
    denom = avg_mw - _BOND_MW_LOSS
    if denom <= 0:
        denom = avg_mw
    return int(max(2, round((mw_target + _BOND_MW_LOSS) / denom)))


# ── Composition / linkage overrides ───────────────────────────────────────────

def _apply_monomer_overrides(base: list, args: argparse.Namespace) -> list:
    """Return normalised [H, G, S] vector after applying CLI flags."""
    H, G, S = base
    if args.monomer_type in ("H", "G", "S"):
        v = [0.0, 0.0, 0.0]
        v[{"H": 0, "G": 1, "S": 2}[args.monomer_type]] = 1.0
        return v
    if args.H_fraction is not None:
        H = args.H_fraction
    if args.G_fraction is not None:
        G = args.G_fraction
    if args.S_fraction is not None:
        S = args.S_fraction
    v = np.array([H, G, S], dtype=float)
    total = v.sum()
    return (v / total if total > 0 else np.array([0.2, 0.5, 0.3])).tolist()


def _apply_linkage_overrides(base: list, args: argparse.Namespace) -> list:
    """Return normalised 7-element linkage vector after applying CLI flags.

    Order matches ``lignoforge.core.rules.linkage_names``:
        0: 4-O-5, 1: alpha-O-4, 2: beta-O-4,
        3: 5-5,   4: beta-5,    5: beta-beta, 6: beta-1
    """
    cli_vals = {
        0: args.lk_4_O_5, 1: args.lk_alpha_O_4, 2: args.lk_beta_O_4,
        3: args.lk_5_5,   4: args.lk_beta_5,     5: args.lk_beta_beta,
        6: args.lk_beta_1,
    }
    v = list(base)
    for idx, val in cli_vals.items():
        if val is not None:
            v[idx] = val
    arr = np.array(v, dtype=float)
    total = arr.sum()
    return (arr / total if total > 0 else np.ones(len(arr)) / len(arr)).tolist()


# ── Chain growth ───────────────────────────────────────────────────────────────

def _grow_chain(trajectory_kwargs: dict, seed: int, i_max: int) -> object:
    """Grow one polymer and return it (stochastic Trajectory engine)."""
    from lignoforge.simulation.trajectory import Trajectory
    traj = Trajectory(**trajectory_kwargs)
    polymer, _dist, _n, _steps = traj.run_MCMC(rseed=seed, i_max=i_max)
    return polymer


def _grow_chain_exact(
    n: int,
    monomer_dist: list,
    linkage_dist: list,
    seed: int,
    branching: Optional[float],
) -> object:
    """Grow a chain of exactly *n* monomers, bypassing the Trajectory engine.

    Uses ``Polymer.add_random_monomer`` in a direct loop so the final chain
    always contains exactly *n* monomers (or fewer if no compatible linkage
    can be found after 20 attempts for a given step).
    """
    from lignoforge.core.monomer import Monomer
    from lignoforge.core.polymer import Polymer
    from lignoforge.core.utils import (
        set_random_state,
        generate_random_monomer,
        generate_random_branching_state,
    )

    rstate   = set_random_state(seed)
    m_dist   = np.asarray(monomer_dist)
    l_dist   = np.asarray(linkage_dist)

    # Initialise with the first monomer
    mtype   = generate_random_monomer(m_dist, rstate)
    m_init  = Monomer(mtype)
    polymer = Polymer(m_init)

    # Grow n-1 additional monomers
    for _ in range(n - 1):
        # Determine branching state for this step
        if branching is not None and branching > 0.0:
            b_state = generate_random_branching_state(branching, rstate)
        else:
            b_state = None

        for _attempt in range(20):
            if polymer.add_random_monomer(
                monomer_distribution=m_dist,
                linkage_distribution=l_dist,
                branching_state=b_state,
                random_state=rstate,
            ):
                break

    return polymer


# ── Format parsing ─────────────────────────────────────────────────────────────

_ALL_FORMATS = {
    "smiles", "population-smiles", "sdf",
    "pdb", "pdb-cg", "json-atomistic", "json-cg", "html",
}


def _parse_formats(fmt_str: str) -> set:
    parts = {s.strip().lower() for s in fmt_str.split(",")}
    if "all" in parts:
        return set(_ALL_FORMATS)
    unknown = parts - _ALL_FORMATS
    if unknown:
        print(f"  [warn] Unknown format(s) ignored: {', '.join(sorted(unknown))}")
    return parts & _ALL_FORMATS


# ── Population export ──────────────────────────────────────────────────────────

def _export_population(
    polymers:          list,
    output_dir:        str,
    base_name:         str,
    formats:           set,
    seeds:             list,
    generate_3d:       bool,
    optimize_3d:       bool,
    include_hydrogens: bool,
    max_iter:          int,
    n_workers:         Optional[int],
    verbose:           bool,
) -> dict:
    """Export an entire polymer population.  Returns {fmt: path | [paths]}."""
    from lignoforge.core.utils import graph_to_smile, graph_to_mol
    from lignoforge.structure.generator import MolecularStructureGenerator
    from lignoforge.structure.pdb import PDBStructureWriter

    gen        = MolecularStructureGenerator()
    pdb_writer = PDBStructureWriter()
    written: dict = {}
    n = len(polymers)

    # ── Per-chain SMILES ───────────────────────────────────────────────────────
    if "smiles" in formats:
        smi_paths = []
        for i, polymer in enumerate(polymers):
            smi    = graph_to_smile(polymer.G)
            suffix = f"_{i:03d}" if i > 0 else ""
            path   = os.path.join(output_dir, f"{base_name}{suffix}.smi")
            Path(path).write_text(smi + "\n")
            smi_paths.append(path)
        written["smiles"] = smi_paths[0] if len(smi_paths) == 1 else smi_paths

    # ── Combined population SMILES (all chains in one file) ───────────────────
    if "population-smiles" in formats:
        path = os.path.join(output_dir, f"{base_name}_population.smi")
        with open(path, "w") as fh:
            for i, polymer in enumerate(polymers):
                smi = graph_to_smile(polymer.G)
                fh.write(f"{smi}\t{base_name}_{i}\n")
        written["population-smiles"] = path
        print(f"    [population-smiles] {path}")

    # ── Formats that require 3-D topology ─────────────────────────────────────
    needs_3d = formats & {"pdb", "pdb-cg", "json-atomistic", "json-cg", "html"}
    if needs_3d:
        if not generate_3d:
            print("    [warn] 3-D formats requested but --no-3d was set; skipping.")
        else:
            workers_label = f"{n_workers} worker(s)" if n_workers else "all CPUs"
            print(f"    Generating atomistic topology  "
                  f"({n} chain(s), {workers_label}, "
                  f"{'no optimize' if not optimize_3d else f'max {max_iter} iter'})...")
            atomistic_pop = gen.population_atomistic_topology(
                polymers,
                random_seed=seeds[0] if seeds else 42,
                include_hydrogens=include_hydrogens,
                optimize_3d=optimize_3d,
                max_uff_iterations=max_iter,
                n_workers=n_workers,
            )

            # Atomistic PDB (per-chain in pdb/ subfolder)
            if "pdb" in formats:
                pdb_dir = os.path.join(output_dir, "pdb")
                paths   = pdb_writer.write_population_pdbs_from_atomistic_population(
                    atomistic_pop, pdb_dir, basename=base_name
                )
                written["pdb"] = paths
                print(f"    [pdb] {pdb_dir}/  ({len(paths)} file(s))")

            # Atomistic topology JSON (population-level)
            if "json-atomistic" in formats:
                path = os.path.join(output_dir, f"{base_name}_atomistic_topology.json")
                with open(path, "w") as fh:
                    json.dump(atomistic_pop, fh, indent=2)
                written["json-atomistic"] = path
                print(f"    [json-atomistic] {path}")

            # CG topology (JSON, HTML viewer, CG PDB)
            if "pdb-cg" in formats or "json-cg" in formats or "html" in formats:
                if verbose:
                    print("    Deriving coarse-grained topology from atomistic CoMs...")
                cg_pop = gen.population_cg_from_atomistic(atomistic_pop)

                cg_json_path = None
                if "json-cg" in formats or "html" in formats:
                    cg_json_path = os.path.join(
                        output_dir, f"{base_name}_cg_topology.json"
                    )
                    with open(cg_json_path, "w") as fh:
                        json.dump(cg_pop, fh, indent=2)
                    written["json-cg"] = cg_json_path
                    print(f"    [json-cg] {cg_json_path}")

                if "html" in formats and cg_json_path is not None:
                    try:
                        from lignoforge.structure.visualization import (
                            write_cg_topology_viewer_html,
                        )
                        html_path = os.path.join(
                            output_dir, f"{base_name}_cg_viewer.html"
                        )
                        write_cg_topology_viewer_html(cg_json_path, html_path)
                        written["html"] = html_path
                        print(f"    [html] {html_path}")
                    except ImportError:
                        print("    [warn] plotly not installed — html viewer skipped.")

                if "pdb-cg" in formats:
                    cg_pdb_dir = os.path.join(output_dir, "pdb_cg")
                    cg_paths   = pdb_writer.write_cg_population_pdbs_from_cg_population(
                        cg_pop, cg_pdb_dir, basename=f"{base_name}_cg"
                    )
                    written["pdb-cg"] = cg_paths
                    print(f"    [pdb-cg] {cg_pdb_dir}/  ({len(cg_paths)} file(s))")

    # ── SDF (per-chain via graph_to_mol, no population topology needed) ────────
    if "sdf" in formats:
        from rdkit import Chem
        sdf_paths = []
        for i, polymer in enumerate(polymers):
            seed_i = seeds[i] if i < len(seeds) else 42
            mol = graph_to_mol(
                polymer.G,
                generate_3d=generate_3d,
                random_seed=seed_i,
                optimize_3d=generate_3d and optimize_3d,
                max_uff_iterations=max_iter,
                add_hs_for_3d=include_hydrogens,
            )
            if mol is not None:
                suffix = f"_{i:03d}" if i > 0 else ""
                path = os.path.join(output_dir, f"{base_name}{suffix}.sdf")
                Chem.MolToMolFile(mol, path)
                sdf_paths.append(path)
        written["sdf"] = sdf_paths
        print(f"    [sdf] {len(sdf_paths)} file(s) written")

    return written


# ── Main ───────────────────────────────────────────────────────────────────────

def main(argv: Optional[list] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    # ── Load and validate JSON input ─────────────────────────────────────────
    if args.input is None:
        input_data  = _DEFAULT_INPUT
        input_label = "<built-in defaults>"
    else:
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"[error] Input file not found: {input_path}", file=sys.stderr)
            return 1
        try:
            from lignoforge.io.schema import InputSchemaValidator
            input_data = InputSchemaValidator().validate_file(str(input_path))
        except Exception as e:
            print(f"[error] Invalid input JSON: {e}", file=sys.stderr)
            return 1
        input_label = str(input_path)

    # ── Estimate structural priors ────────────────────────────────────────────
    from lignoforge.priors.estimator import LigninPriorEstimator
    from lignoforge.pipeline.translator import ParameterTranslator

    priors     = LigninPriorEstimator(input_data, random_seed=args.seed).run()
    sim_kwargs = ParameterTranslator(priors, input_data).to_simulation_kwargs()

    # ── Apply user overrides ──────────────────────────────────────────────────
    monomer_dist = _apply_monomer_overrides(sim_kwargs["monomer_distribution_input"], args)
    linkage_dist = _apply_linkage_overrides(sim_kwargs["linkage_distribution_input"], args)

    if args.mw_target is not None:
        target_dp = _dp_from_mw(args.mw_target, monomer_dist)
        print(f"  MW target {args.mw_target:.0f} g/mol → estimated DP = {target_dp}")
    elif args.n_monomers is not None:
        target_dp = args.n_monomers
    else:
        target_dp = int(sim_kwargs.get("expected_size", 20))

    max_dp = max(target_dp + 5, int(round(target_dp * 1.5)))
    Tmetro = args.Tmetro if args.Tmetro is not None else float(sim_kwargs.get("Tmetro", 298.15))

    trajectory_kwargs = {
        "linkage_distribution_input": linkage_dist,
        "monomer_distribution_input": monomer_dist,
        "expected_size":              float(target_dp),
        "max_size":                   float(max_dp),
        "distribution_scaling":       float(sim_kwargs.get("distribution_scaling", 1.0)),
        "Tmetro":                     Tmetro,
        "branching_propensity":       args.branching,
        "verbose":                    args.verbose,
    }

    # ── Prepare output ────────────────────────────────────────────────────────
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    formats      = _parse_formats(args.formats)
    generate_3d  = not args.no_3d
    optimize_3d  = not args.no_optimize
    include_hs   = not args.no_H

    # ── Print run header ──────────────────────────────────────────────────────
    lnames = ["4-O-5", "a-O-4", "b-O-4", "5-5", "b-5", "b-b", "b-1"]
    print(f"\nlignoforge-chain  [{args.n_chains} chain(s)]")
    print(f"  Input      : {input_label}")
    print(f"  Target DP  : {target_dp}  (max={max_dp})")
    print(f"  Composition: H={monomer_dist[0]:.3f}  G={monomer_dist[1]:.3f}  S={monomer_dist[2]:.3f}")
    print(f"  Linkages   : " + "  ".join(f"{n}={v:.2f}" for n, v in zip(lnames, linkage_dist)))
    print(f"  Branching  : {args.branching or 0.0:.3f}")
    print(f"  Tmetro     : {Tmetro}")
    print(f"  3-D embed  : {'yes' if generate_3d else 'no'}  "
          f"{'(no optimize)' if generate_3d and not optimize_3d else ''}"
          f"  workers={args.n_workers or 'all CPUs'}")
    print(f"  Formats    : {', '.join(sorted(formats))}")
    print(f"  Output     : {output_dir}/")
    print()

    # ── Phase 1: grow all chains ───────────────────────────────────────────────
    from lignoforge.core.characterization import Characterize

    polymers:   list = []
    seeds:      list = []
    all_stats:  list = []

    for i in range(args.n_chains):
        chain_seed = args.seed + i
        print(f"  Chain {i + 1}/{args.n_chains}  (seed={chain_seed})")

        try:
            if args.n_monomers is not None or args.mw_target is not None:
                polymer = _grow_chain_exact(
                    n=target_dp,
                    monomer_dist=monomer_dist,
                    linkage_dist=linkage_dist,
                    seed=chain_seed,
                    branching=args.branching,
                )
            else:
                polymer = _grow_chain(
                    trajectory_kwargs=trajectory_kwargs,
                    seed=chain_seed,
                    i_max=args.max_steps,
                )
        except Exception as e:
            print(f"    [error] Chain growth failed: {e}", file=sys.stderr)
            traceback.print_exc()
            continue

        stats = Characterize(polymer).summary()
        stats["chain_index"] = i
        stats["seed"]        = chain_seed
        all_stats.append(stats)
        polymers.append(polymer)
        seeds.append(chain_seed)
        print(f"    n_monomers={stats['monomer_count']}  MW={stats['MW']:.0f} g/mol")

    if not polymers:
        print("[error] No chains could be grown.", file=sys.stderr)
        return 1

    # ── Phase 2: population-level export ──────────────────────────────────────
    print(f"\n  Exporting {len(polymers)} chain(s)...")
    try:
        _export_population(
            polymers=polymers,
            output_dir=str(output_dir),
            base_name=args.name,
            formats=formats,
            seeds=seeds,
            generate_3d=generate_3d,
            optimize_3d=optimize_3d,
            include_hydrogens=include_hs,
            max_iter=args.max_iter,
            n_workers=args.n_workers,
            verbose=args.verbose,
        )
    except Exception as e:
        print(f"    [error] Export failed: {e}", file=sys.stderr)
        traceback.print_exc()

    # ── Phase 3: statistics and metadata ─────────────────────────────────────
    stats_path = output_dir / f"{args.name}_chain_stats.json"
    with open(stats_path, "w") as f:
        json.dump(all_stats, f, indent=2)
    print(f"\n  [chain-stats]    {stats_path}")

    if args.population_stats:
        from lignoforge.io.exporters import LigninExporter
        exp = LigninExporter(str(output_dir))
        pop_stats_path = exp.export_population_statistics(
            polymers, filename=f"{args.name}_population_stats.json"
        )
        print(f"  [population-stats] {pop_stats_path}")

    if args.export_priors:
        priors_path = output_dir / f"{args.name}_estimated_priors.json"
        with open(priors_path, "w") as f:
            json.dump(priors, f, indent=2)
        print(f"  [priors]         {priors_path}")

    if args.export_sim_params:
        sim_path = output_dir / f"{args.name}_simulation_parameters.json"
        with open(sim_path, "w") as f:
            json.dump(sim_kwargs, f, indent=2)
        print(f"  [sim-params]     {sim_path}")

    # ── Run manifest ──────────────────────────────────────────────────────────
    manifest = {
        "input_file":   input_label,
        "target_dp":    target_dp,
        "monomer_dist": {"H": monomer_dist[0], "G": monomer_dist[1], "S": monomer_dist[2]},
        "linkage_dist": dict(zip(
            ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta", "beta-1"],
            linkage_dist,
        )),
        "branching":    args.branching or 0.0,
        "seed":         args.seed,
        "n_chains":     args.n_chains,
        "Tmetro":       Tmetro,
        "generate_3d":  generate_3d,
        "optimize_3d":  optimize_3d,
        "n_workers":    args.n_workers,
        "formats":      sorted(formats),
        "chains":       all_stats,
    }
    manifest_path = output_dir / f"{args.name}_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  [manifest]       {manifest_path}")

    n_ok = len(all_stats)
    print(f"\nDone. {n_ok}/{args.n_chains} chain(s) written to {output_dir}/\n")
    return 0 if n_ok > 0 else 1


def _module_main() -> None:
    sys.exit(main())


if __name__ == "__main__":
    _module_main()
