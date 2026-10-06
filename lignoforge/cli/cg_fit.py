#!/usr/bin/env python
"""
lignoforge-cg-fit — derive CG bond/angle parameters from atomistic MD
=====================================================================

Boltzmann inversion of the monomer-centre-of-mass distributions of one or
more atomistic runs produced by ``lignoforge-chain --format md``.

Usage
-----
    lignoforge-cg-fit CHAINS.json --run 0:frames0.pdb [--run 1:frames1.pdb ...]
                      [--temperature 300] [-o cg_parameters.json]

``CHAINS.json`` is the ``*_atomistic_topology.json`` written with
``--format json-atomistic``.  Each ``--run INDEX:FRAMES`` pairs a chain of that
file with a multi-model PDB containing only that chain, e.g.::

    echo "a 1-<n_atoms>
    q" | gmx make_ndx -f chain.gro -o chain.ndx
    echo "a_1-<n_atoms>" | gmx trjconv -s prod.tpr -f prod.xtc -n chain.ndx \\
                             -pbc mol -o frames.pdb

Samples are pooled by linkage (bonds) and linkage pair (angles) across all
runs; entries with too few samples keep the provisional defaults.  Pass the
output to ``lignoforge-chain --format cg --cg-params cg_parameters.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Optional


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(prog="lignoforge-cg-fit", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("chains_json")
    ap.add_argument("--run", action="append", required=True, metavar="INDEX:FRAMES")
    ap.add_argument("--temperature", type=float, default=300.0)
    ap.add_argument("--min-samples", type=int, default=200)
    ap.add_argument("-o", "--output", default="cg_parameters.json")
    args = ap.parse_args(argv)

    from lignoforge.cg import BoltzmannFitter, map_topology, read_pdb_frames
    from lignoforge.forcefield import build_chain_topology

    with open(args.chains_json) as fh:
        data = json.load(fh)
    chains = data.get("chains", [data])

    fitter = BoltzmannFitter()
    for spec in args.run:
        try:
            idx_s, path = spec.split(":", 1)
            chain = chains[int(idx_s)]
        except (ValueError, IndexError):
            print(f"[error] bad --run '{spec}' (expected INDEX:FRAMES, "
                  f"index < {len(chains)})", file=sys.stderr)
            return 1
        topo = build_chain_topology(chain)
        cg = map_topology(topo)
        frames = read_pdb_frames(path, len(topo.atoms))
        n = fitter.add(cg, frames)
        print(f"  chain {idx_s}: {n} frames, {len(cg.bonds)} bonds, {len(cg.angles)} angles")

    params, report = fitter.fit(args.temperature, min_samples=args.min_samples)
    params.save(args.output)
    print(f"\n  fitted : {len(report['fitted'])} entries")
    for k, v in sorted(report["fitted"].items()):
        print(f"    {k:34s} " + "  ".join(f"{a}={b}" for a, b in v.items()))
    if report["kept_default"]:
        print(f"  kept default ({len(report['kept_default'])}): "
              + ", ".join(report["kept_default"]))
    print(f"\n  wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
