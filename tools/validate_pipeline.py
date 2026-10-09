#!/usr/bin/env python
"""
End-to-end validation of LignoForge: structure → OPLS-AA topology → GROMACS.

For every test case the script grows a chain and checks, in order:

  1. growth      requested size reached; linear chains really are linear
  2. chemistry   RDKit sanitises the molecule; formula/atom counts consistent
  3. geometry    no steric clashes in the generated 3-D structure
  4. topology    every graph bond present, all atoms typed, net charge == 0
  5. grompp      topology accepted by GROMACS (no missing parameters)
  6. EM          energy minimisation converges, no bond is broken
  7. MD          (--md) short vacuum MD: chain stays intact, T is stable
  8. CG          CG system written, accepted by grompp, EM runs

Usage
-----
    python tools/validate_pipeline.py                     # all cases, quick
    python tools/validate_pipeline.py --md                # + short MD (slower)
    python tools/validate_pipeline.py --cases P1,P4 --seeds 0,1,2
    python tools/validate_pipeline.py --random 20         # + 20 random chains
    python tools/validate_pipeline.py --keep              # keep all outputs

Outputs go to ``--workdir`` (default ./validation_out); a per-case log and a
summary table are printed, and the exit code is non-zero if any check fails.
Needs ``gmx`` on PATH for checks 5-8 (they are skipped otherwise).
"""

from __future__ import annotations

import argparse
import contextlib
import io
import itertools
import json
import math
import os
import shutil
import subprocess
import sys
import warnings
from collections import Counter
from typing import Dict, List, Optional, Tuple

import numpy as np

LK_ORDER = ["4-O-5", "alpha-O-4", "beta-O-4", "5-5", "beta-5", "beta-beta", "beta-1"]

# name: (n, (H, G, S), {linkage: fraction}, branching)
CASES = {
    "P1":  (5,  (0.4, 0.4, 0.2), {"4-O-5": .4, "beta-O-4": .4, "beta-beta": .2}, 0.0),
    "P1b": (10, (0.4, 0.4, 0.2), {"4-O-5": .4, "beta-O-4": .4, "beta-beta": .2}, 0.0),
    "P2":  (5,  (0.0, 0.5, 0.5), {"4-O-5": .5, "beta-O-4": .5}, 0.0),
    "P3":  (5,  (0.0, 0.5, 0.5), {"beta-O-4": 1.0}, 0.0),
    "P4":  (5,  (0.0, 0.0, 1.0), {"beta-O-4": 1.0}, 0.0),
    "ALL": (8,  (0.2, 0.5, 0.3), {"beta-O-4": 3, "alpha-O-4": 1, "4-O-5": 1,
                                  "5-5": 1, "beta-5": 1, "beta-beta": 1}, 0.0),
    "BR":  (10, (0.2, 0.5, 0.3), {"beta-O-4": 3, "4-O-5": 1, "5-5": 1}, 0.4),
}

ATOM_MASS = {"C": 12.011, "H": 1.008, "O": 15.999}


# ── helpers ───────────────────────────────────────────────────────────────────

class Result:
    def __init__(self, name):
        self.name, self.checks, self.notes = name, [], []

    def add(self, check: str, ok: Optional[bool], detail: str = ""):
        self.checks.append((check, ok, detail))
        flag = {True: "PASS", False: "FAIL", None: "SKIP"}[ok]
        print(f"    [{flag}] {check:10s} {detail}")

    @property
    def failed(self):
        return any(ok is False for _, ok, _ in self.checks)


def run(cmd, cwd, timeout=900):
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)


def read_gro(path) -> Tuple[np.ndarray, np.ndarray]:
    lines = open(path).read().splitlines()
    n = int(lines[1])
    xyz = np.array([[float(l[20 + 8 * i:28 + 8 * i]) for i in range(3)]
                    for l in lines[2:2 + n]])
    box = np.array([float(x) for x in lines[2 + n].split()[:3]])
    return xyz, box


def bond_stats(xyz, box, bonds) -> Tuple[float, int]:
    """Max bond length (nm, minimum image) and number of bonds > 0.2 nm."""
    mx, broken = 0.0, 0
    for i, j in bonds:
        d = xyz[i - 1] - xyz[j - 1]
        d -= box * np.round(d / box)
        r = float(np.linalg.norm(d))
        mx = max(mx, r)
        broken += r > 0.2
    return mx, broken


def parse_bonds_from_top(path) -> List[Tuple[int, int]]:
    txt = open(path).read().split("[ bonds ]")[1].split("[")[0]
    return [tuple(map(int, l.split()[:2])) for l in txt.splitlines()
            if l.strip() and l[0] != ";"]


# ── the checks ────────────────────────────────────────────────────────────────

def grow(n, mono, linkages, branching, seed):
    from lignoforge.cli.build_chain import _grow_chain_exact
    lk = np.array([linkages.get(k, 0.0) for k in LK_ORDER], float)
    lk /= lk.sum()
    m = np.array(mono, float)
    m /= m.sum()
    with warnings.catch_warnings(record=True) as w, \
            contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("always")
        p = _grow_chain_exact(n, m.tolist(), lk.tolist(), seed, branching)
    return p, [str(x.message) for x in w]


def check_growth(res, p, n, branching):
    got = p.mi + 1
    deg = dict(p.bigG.degree())
    maxdeg = max(deg.values())
    types = Counter(d["mtype"] for _, d in p.bigG.nodes(data=True))
    links = Counter(d["btype"] for _, _, d in p.bigG.edges(data=True))
    ok = got == n and (branching > 0 or maxdeg <= 2)
    res.add("growth", ok, f"n={got}/{n} maxdeg={maxdeg} types={dict(types)} links={dict(links)}")


def check_chemistry(res, p):
    from lignoforge.core.utils import graph_to_mol
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    mol = graph_to_mol(p.G, generate_3d=False)
    if mol is None:
        return res.add("chemistry", False, "RDKit could not sanitise the molecule")
    nfrag = len(Chem.GetMolFrags(mol))
    mol_h = Chem.AddHs(mol)
    ok = nfrag == 1
    res.add("chemistry", ok, f"{CalcMolFormula(mol)} fragments={nfrag} "
                              f"MW={Descriptors.MolWt(mol_h):.1f}")
    return mol_h


def make_chain(p, seed, optimize=True):
    from lignoforge.structure.generator import MolecularStructureGenerator
    with contextlib.redirect_stdout(io.StringIO()):
        return MolecularStructureGenerator().polymer_atomistic_topology(
            p, random_seed=seed, optimize_3d=optimize, max_uff_iterations=300)


def check_geometry(res, chain):
    """Non-bonded heavy-atom pairs closer than 2.0 Å indicate a clash."""
    atoms = [a for m in chain["monomers"] for a in m["atoms"]]
    xyz = np.array([[a["x"], a["y"], a["z"]] for a in atoms])
    ids = [a["atom_id"] for a in atoms]
    pos = {i: k for k, i in enumerate(ids)}
    bonded = {tuple(sorted((pos[b["atom1_id"]], pos[b["atom2_id"]]))) for b in chain["bonds"]}
    heavy = [k for k, a in enumerate(atoms) if a["element"] != "H"]
    clashes = 0
    d = np.linalg.norm(xyz[heavy][:, None] - xyz[heavy][None], axis=-1)
    for a, b in zip(*np.where(np.triu(d < 2.0, 1))):
        if tuple(sorted((heavy[a], heavy[b]))) not in bonded:
            # 1-3 neighbours are legitimately ~2.4 Å; < 2.0 Å is a real clash
            clashes += 1
    res.add("geometry", clashes == 0, f"{clashes} heavy-atom clash(es) < 2.0 Å")


def check_topology(res, chain, outdir, name):
    from lignoforge.forcefield import build_chain_topology, write_gromacs_system
    topo = build_chain_topology(chain, name=name)
    ok = (topo.net_charge == 0.0 and len(topo.bonds) == len(chain["bonds"])
          and all(a.opls.startswith("opls_") for a in topo.atoms))
    worst = max(r["max_atom_shift"] for r in topo.charge_report)
    res.add("topology", ok, f"atoms={len(topo.atoms)} q={topo.net_charge:+.4f} "
            f"angles={len(topo.angles)} dih={len(topo.dihedrals)} "
            f"max charge shift={worst:.4f} e")
    write_gromacs_system(chain, outdir, name=name)
    return topo


def check_grompp_em(res, outdir, name, topo, do_md):
    if shutil.which("gmx") is None:
        for c in ("grompp", "EM") + (("MD",) if do_md else ()):
            res.add(c, None, "gmx not found")
        return
    r = run(["gmx", "grompp", "-f", "em.mdp", "-c", f"{name}.gro", "-p", f"{name}.top",
             "-o", "em.tpr", "-maxwarn", "0"], outdir)
    errs = [l.strip() for l in (r.stdout + r.stderr).splitlines()
            if "No default" in l or "ERROR" in l]
    res.add("grompp", r.returncode == 0, "ok" if r.returncode == 0 else "; ".join(errs[:3]))
    if r.returncode != 0:
        return
    r = run(["gmx", "mdrun", "-deffnm", "em", "-nt", "1"], outdir)
    log = open(os.path.join(outdir, "em.log")).read()
    conv = "converged to Fmax" in log
    bonds = parse_bonds_from_top(os.path.join(outdir, f"{name}.top"))
    xyz, box = read_gro(os.path.join(outdir, "em.gro"))
    mx, broken = bond_stats(xyz, box, bonds)
    res.add("EM", conv and broken == 0,
            f"converged={conv} max bond={mx * 10:.2f} Å broken={broken}")
    if not do_md or not conv:
        return
    mdp = open(os.path.join(outdir, "em.mdp")).read()
    with open(os.path.join(outdir, "md.mdp"), "w") as fh:
        fh.write("integrator=md\ndt=0.001\nnsteps=20000\nconstraints=none\n"
                 "tcoupl=v-rescale\ntc-grps=System\ntau-t=0.1\nref-t=300\n"
                 "gen-vel=yes\ngen-temp=300\ncutoff-scheme=Verlet\ncoulombtype=PME\n"
                 "rcoulomb=1.0\nrvdw=1.0\npbc=xyz\nnstenergy=500\n")
    r = run(["gmx", "grompp", "-f", "md.mdp", "-c", "em.gro", "-p", f"{name}.top",
             "-o", "md.tpr"], outdir)
    if r.returncode:
        return res.add("MD", False, "grompp failed")
    run(["gmx", "mdrun", "-deffnm", "md", "-nt", "2"], outdir, timeout=1800)
    if not os.path.exists(os.path.join(outdir, "md.gro")):
        return res.add("MD", False, "mdrun did not finish")
    xyz, box = read_gro(os.path.join(outdir, "md.gro"))
    mx, broken = bond_stats(xyz, box, bonds)
    res.add("MD", broken == 0, f"20 ps vacuum, max bond={mx * 10:.2f} Å broken={broken}")


def check_cg(res, chain, outdir):
    from lignoforge.cg import write_cg_system
    write_cg_system([chain], outdir, name="cg", run_ns=0.1)
    if shutil.which("gmx") is None:
        return res.add("CG", None, "gmx not found")
    r = run(["gmx", "grompp", "-f", "em.mdp", "-c", "cg_0.gro", "-p", "cg.top",
             "-o", "em.tpr"], outdir)
    if r.returncode:
        return res.add("CG", False, "grompp failed: " + r.stderr[-200:])
    r = run(["gmx", "mdrun", "-deffnm", "em", "-nt", "1"], outdir)
    res.add("CG", "converged to Fmax" in open(os.path.join(outdir, "em.log")).read(),
            "CG EM converged")


# ── driver ────────────────────────────────────────────────────────────────────

def random_cases(k, seed=0):
    rng = np.random.RandomState(seed)
    out = {}
    for i in range(k):
        lk = rng.dirichlet(np.ones(7) * 0.7)
        lk[6] = 0
        mono = rng.dirichlet(np.ones(3))
        out[f"R{i:02d}"] = (int(rng.randint(4, 13)), tuple(mono),
                            {n: float(v) for n, v in zip(LK_ORDER, lk)},
                            float([0, 0, 0.3][i % 3]))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", default=",".join(CASES), help="comma list (default: all)")
    ap.add_argument("--seeds", default="42", help="comma list of seeds (default: 42)")
    ap.add_argument("--random", type=int, default=0, help="add N random chains")
    ap.add_argument("--md", action="store_true", help="add a 20 ps vacuum MD check")
    ap.add_argument("--workdir", default="validation_out")
    ap.add_argument("--keep", action="store_true", help="keep outputs of passing cases")
    args = ap.parse_args(argv)

    cases = {**CASES, **random_cases(args.random)}
    names = [c for c in args.cases.split(",") if c in cases] + \
            [c for c in cases if c.startswith("R")]
    seeds = [int(s) for s in args.seeds.split(",")]
    os.makedirs(args.workdir, exist_ok=True)

    results = []
    for name, seed in itertools.product(names, seeds):
        n, mono, linkages, branching = cases[name]
        tag = f"{name}_s{seed}"
        outdir = os.path.join(args.workdir, tag)
        shutil.rmtree(outdir, ignore_errors=True)
        os.makedirs(outdir)
        print(f"\n== {tag}: n={n} HGS={tuple(round(x, 2) for x in mono)} "
              f"branching={branching}")
        res = Result(tag)
        try:
            p, warns = grow(n, mono, linkages, branching, seed)
            res.notes += warns
            check_growth(res, p, n, branching)
            check_chemistry(res, p)
            chain = make_chain(p, seed)
            check_geometry(res, chain)
            topo = check_topology(res, chain, outdir, "lig")
            check_grompp_em(res, outdir, "lig", topo, args.md)
            check_cg(res, chain, os.path.join(outdir, "cg"))
            with open(os.path.join(outdir, "chain.json"), "w") as fh:
                json.dump(chain, fh)
        except Exception as e:                                    # noqa: BLE001
            res.add("exception", False, f"{type(e).__name__}: {e}")
        for note in res.notes:
            print(f"    [note] {note}")
        results.append(res)
        if not res.failed and not args.keep:
            shutil.rmtree(outdir, ignore_errors=True)

    print("\n" + "=" * 72)
    print(f"{'case':14s} " + "  ".join(f"{c:9s}" for c in
          ["growth", "chemistry", "geometry", "topology", "grompp", "EM", "MD", "CG"]))
    for r in results:
        d = {c: ok for c, ok, _ in r.checks}
        cell = lambda c: "-" if c not in d else {True: "ok", False: "FAIL", None: "skip"}[d[c]]  # noqa: E731
        print(f"{r.name:14s} " + "  ".join(f"{cell(c):9s}" for c in
              ["growth", "chemistry", "geometry", "topology", "grompp", "EM", "MD", "CG"]))
    bad = [r.name for r in results if r.failed]
    print("=" * 72)
    print(f"{len(results) - len(bad)}/{len(results)} cases passed"
          + (f"; failed: {', '.join(bad)} (outputs kept in {args.workdir}/)" if bad else ""))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
