#!/usr/bin/env python
"""
LignoForge Demo Script
======================

Generate a small lignin chain ensemble to demonstrate the framework's capabilities.
Output files are saved to demo_output/ and do not interfere with repository files.

Usage:
    python examples/demo_run.py [--n-chains N] [--seed SEED] [--output DIR]

Options:
    --n-chains N    Number of chains to generate (default: 2)
    --seed SEED     Random seed for reproducibility (default: 42)
    --output DIR    Output directory (default: demo_output/)

Example:
    python examples/demo_run.py --n-chains 5 --seed 42
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Run LignoForge demo: generate lignin chain ensemble"
    )
    parser.add_argument("--n-chains", type=int, default=2,
                        help="Number of chains to generate (default: 2)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed (default: 42)")
    parser.add_argument("--output", type=str, default="demo_output",
                        help="Output directory (default: demo_output/)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Parallel worker processes for 3-D embedding "
                             "(default: all CPU cores)")

    args = parser.parse_args()

    # Setup output directory
    output_dir = Path(args.output)
    if output_dir.exists():
        print(f"Output directory '{output_dir}' exists. Cleaning...")
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"✓ Output directory created: {output_dir}/\n")

    # Import pipeline
    try:
        from lignoforge.pipeline import LigninPipeline
    except ImportError:
        print("ERROR: Could not import LignoForge. Make sure it's installed:")
        print("    pip install -e .")
        sys.exit(1)

    # Load input
    input_file = Path(__file__).parent / "lignin_input_example.json"
    if not input_file.exists():
        print(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)

    print(f"Loading input from: {input_file}")
    with open(input_file, "r") as f:
        demo_input = json.load(f)

    # Create pipeline with custom output directory
    print("\nInitializing pipeline...")
    try:
        pipeline = LigninPipeline.from_dict(demo_input, output_dir=str(output_dir))
    except Exception as e:
        print(f"ERROR: Could not create pipeline: {e}")
        sys.exit(1)

    # Run simulation with custom n_population
    print(f"Running simulation ({args.n_chains} chains, seed={args.seed})...")
    print("-" * 60)

    try:
        results = pipeline.run(
            random_seed=args.seed,
            library_name="demo_library",
            simulation_overrides={
                "n_population": args.n_chains,
                "n_workers":    args.workers,
            },
        )
    except Exception as e:
        print(f"ERROR during simulation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print("-" * 60)
    print("\n⏳ Note: 3-D coordinate generation uses all available CPU cores by default.")
    print("   To limit parallelism: python examples/demo_run.py --workers 4")
    print("   To disable 3-D generation: pass export_options={'generate_3d': False}\n")

    # Summary
    print("Generated artifacts:")
    print("-" * 60)
    if hasattr(results, 'artifacts') and results.artifacts:
        for artifact_type, file_path in results.artifacts.items():
            if file_path is None:
                continue
            if isinstance(file_path, list):
                for p in file_path:
                    if p and os.path.exists(p):
                        rel_path = os.path.relpath(p, output_dir)
                        file_size = os.path.getsize(p)
                        print(f"  {artifact_type:30s}: {rel_path:40s} ({file_size:,d} bytes)")
            elif os.path.exists(file_path):
                rel_path = os.path.relpath(file_path, output_dir)
                file_size = os.path.getsize(file_path)
                print(f"  {artifact_type:30s}: {rel_path:40s} ({file_size:,d} bytes)")
    else:
        print("  (no artifacts generated)")

    print("-" * 60)
    print(f"\nOutput directory: {output_dir}/")
    print("\nKey files to explore:")

    # Find PDB files
    pdb_files = list(output_dir.glob("pdb/*.pdb"))
    if pdb_files:
        print(f"  • PDB structures: {len(pdb_files)} chain(s)")
        print(f"    Visualization: PyMOL, UCSF Chimera, or online viewers")

    # Find JSON topology files
    atomistic_files = list(output_dir.glob("atomistic_topology_*.json"))
    if atomistic_files:
        print(f"  • Atomistic topologies: {len(atomistic_files)} chain(s)")
        print(f"    Structure: chains → monomers → atoms with 3D coordinates")

    cg_files = list(output_dir.glob("coarse_grained_topology_*.json"))
    if cg_files:
        print(f"  • Coarse-grained topologies: {len(cg_files)} chain(s)")
        print(f"    Structure: beads (monomers) and connectivity")

    # Find HTML viewer
    html_files = list(output_dir.glob("**/cg_topology_viewer_*.html"))
    if html_files:
        print(f"  • Interactive 3D viewers: {len(html_files)} chain(s)")
        print(f"    Open in browser: {html_files[0].name}")

    # Find statistics
    stats_file = output_dir / "population_statistics.json"
    if stats_file.exists():
        print(f"  • Population statistics: population_statistics.json")
        with open(stats_file) as f:
            stats = json.load(f)
            if "ensemble" in stats:
                num_chains = stats['ensemble'].get('n_chains', 'N/A')
                mean_mw = stats['ensemble'].get('mean_molecular_weight', 'N/A')
                print(f"    └─ Chains: {num_chains}, Mean MW: {mean_mw:.0f} g/mol" if isinstance(mean_mw, float) else f"    └─ Chains: {num_chains}")

    print("\n" + "=" * 60)
    print("Demo complete! All outputs are in '%s/' directory." % output_dir)
    print("=" * 60 + "\n")


if __name__ == "__main__":
    main()
