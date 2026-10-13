#!/usr/bin/env python
"""
lignoforge-cg-param — derive a CG model of bulk lignin by running atomistic and CG MD
=====================================================================================

    lignoforge-cg-param init   WORKDIR [--quick] [--scheme ring-side|monomer] [--temperature K]
    lignoforge-cg-param run    WORKDIR [--stage build,reference,map,iterate,validate,export]
                                       [--force] [--gmx PATH] [--dry-run]
    lignoforge-cg-param selftest WORKDIR [--iterations 8] [--mixing free|lb]
    lignoforge-cg-param status WORKDIR [WORKDIR ...] [--follow]   progress and time remaining
    lignoforge-cg-param log    WORKDIR [--md] [-n 40] [--follow]  job log / current mdrun log
    lignoforge-cg-param slurm  WORKDIR [--partition P] [--gres gpu:1] [--time 48:00:00] [--cpus 8]

``init`` writes ``WORKDIR/config.json`` (edit it: ensemble, run lengths, cutoff, ...).
``run`` executes the stages with GROMACS and is resumable.  Use ``--quick`` for
a minutes-long smoke test of the whole chain of tools.  Run directories hold
all simulation data; keep them outside the repository.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(prog="lignoforge-cg-param", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("init", help="write a complete config.json for this machine")
    p.add_argument("workdir")
    p.add_argument("--quick", action="store_true", help="tiny, fast settings for tool validation")
    p.add_argument("--size", choices=["small", "medium", "large"], default=None,
                   help="ensemble size: 16 / 24 (default) / 40 chains")
    p.add_argument("--scheme", default=None, help="ring-side (default) or monomer")
    p.add_argument("--temperature", type=float, default=None, help="K (default 500)")
    p.add_argument("--mixing", choices=["free", "lb"], default=None)
    p.add_argument("--monomers", type=int, default=None, help="degree of polymerisation per chain")
    p.add_argument("--chains", type=int, default=None, help="distinct chain sequences")
    p.add_argument("--copies", type=int, default=None, help="copies of each sequence")
    p.add_argument("--H", type=float, default=None); p.add_argument("--G", type=float, default=None)
    p.add_argument("--S", type=float, default=None, help="monomer fractions (give all three)")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--compress-ns", type=float, default=None)
    p.add_argument("--equil-ns", type=float, default=None)
    p.add_argument("--prod-ns", type=float, default=None)
    p.add_argument("--iterations", type=int, default=None)
    p.add_argument("--iter-ns", type=float, default=None, help="CG time per iteration")
    p.add_argument("--cutoff", type=float, default=None, help="CG cut-off (nm)")
    p.add_argument("--gpu", choices=["auto", "yes", "no"], default="auto",
                   help="use GPU in the atomistic runs (auto: if the gmx found has GPU support)")
    p.add_argument("--threads", type=int, default=None, help="OpenMP threads (default: all cores)")
    p.add_argument("--gmx", default=None, help="path of the gmx binary (default: PATH / $GMX)")
    p.add_argument("--force", action="store_true", help="overwrite an existing config.json")
    p = sub.add_parser("run")
    p.add_argument("workdir")
    p.add_argument("--stage", default=None, help="comma-separated subset of stages")
    p.add_argument("--force", action="store_true")
    p.add_argument("--gmx", default=None)
    p.add_argument("--dry-run", action="store_true", help="write the reference commands only")
    p = sub.add_parser("selftest", help="check the CG optimiser on a synthetic reference (minutes)")
    p.add_argument("workdir")
    p.add_argument("--iterations", type=int, default=8)
    p.add_argument("--mixing", default="free", choices=["free", "lb"])
    p.add_argument("--gmx", default=None)
    p = sub.add_parser("status", help="stages, current mdrun, progress and time remaining")
    p.add_argument("workdir", nargs="+", help="one or more work directories")
    p.add_argument("--follow", "-f", action="store_true", help="refresh until interrupted")
    p.add_argument("--interval", type=float, default=10.0, help="refresh period (s) with --follow")
    p = sub.add_parser("log", help="show the job log or the GROMACS log of the current mdrun")
    p.add_argument("workdir")
    p.add_argument("--md", action="store_true", help="the log of the mdrun in progress instead")
    p.add_argument("-n", type=int, default=40, help="lines")
    p.add_argument("--follow", "-f", action="store_true")
    p = sub.add_parser("slurm")
    p.add_argument("workdir")
    p.add_argument("--partition", default=None)
    p.add_argument("--gres", default="gpu:1")
    p.add_argument("--time", default="48:00:00")
    p.add_argument("--cpus", type=int, default=8)
    args = ap.parse_args(argv)

    from lignoforge.cgfit import STAGES, CGParamJob, quick_config

    wd = Path(args.workdir[0] if isinstance(args.workdir, list) else args.workdir)
    if args.cmd == "init":
        from lignoforge.cgfit.setup import make_config
        if (wd / "config.json").exists() and not args.force:
            print(f"{wd / 'config.json'} exists (use --force to overwrite)")
            return 1
        fr = {"H": args.H, "G": args.G, "S": args.S}
        given = {k: v for k, v in fr.items() if v is not None}
        cfg, info = make_config(
            quick=args.quick, scheme=args.scheme, temperature=args.temperature, mixing=args.mixing,
            size=args.size, n_distinct=args.chains, copies=args.copies, monomers=args.monomers,
            fractions=given or None, seed=args.seed, gpu=args.gpu, threads=args.threads, gmx=args.gmx,
            compress_ns=args.compress_ns, equil_ns=args.equil_ns, prod_ns=args.prod_ns,
            iterations=args.iterations, iter_ns=args.iter_ns, cutoff=args.cutoff)
        wd.mkdir(parents=True, exist_ok=True)
        (wd / "config.json").write_text(json.dumps(cfg, indent=2))
        g = info["gmx"]
        print(f"wrote {wd / 'config.json'}")
        print(f"  GROMACS : {g['path'] or 'NOT FOUND'}" + (f" {g['version']}, GPU support: {g['gpu_kind']}" if g["path"] else ""))
        print(f"  GPU     : {'used in atomistic runs' if info['gpu_used'] else 'not used'}; threads {info['threads']}")
        print(f"  system  : {info['chains']} chains, ~{info['atoms_estimate']} atoms, box ~{info['box_nm_estimate']} nm")
        at = cfg["atomistic"]
        print(f"  runs    : atomistic {at['compress_ps']:.0f}+{at['equil_ps']:.0f}+{at['prod_ps']:.0f} ps; "
              f"CG {cfg['cg']['iterations']} x {cfg['cg']['iter_ps']:.0f} ps; T = {cfg['temperature']} K; "
              f"{cfg['scheme']}, mixing {cfg['cg']['mixing']}")
        for w in info["warnings"]:
            print(f"  WARNING : {w}")
        if not g["path"]:
            print("  WARNING : gmx not found; set --gmx or add GROMACS to PATH before `run`")
        print(f"\nnext: lignoforge-cg-param run {wd}")
        return 0
    if args.cmd == "selftest":
        from lignoforge.cgfit.selftest import run_selftest
        r = run_selftest(wd, args.iterations, mixing=args.mixing, gmx=args.gmx, log=print)
        ok = r["rmse_last"] < 0.4 * r["rmse_first"] and r["rmse_last"] < 0.1
        print(f"\nRMSE(g): {r['rmse_first']:.3f} -> {r['rmse_last']:.3f}; max sigma error "
              f"{100 * r['max_rel_sigma_error']:.1f} %; {'PASS' if ok else 'FAIL'}")
        return 0 if ok else 1
    if args.cmd == "run":
        job = CGParamJob(wd, gmx=args.gmx, dry_run=args.dry_run, log=print)
        stages = args.stage.split(",") if args.stage else None
        job.run(stages, force=args.force)
        return 0
    if args.cmd == "status":
        from lignoforge.cgfit.monitor import format_overview
        import time
        while True:
            if args.follow:
                print("\033[2J\033[H", end="")
            print("\n\n".join(format_overview(w) for w in args.workdir))
            if not args.follow:
                return 0
            try:
                time.sleep(args.interval)
            except KeyboardInterrupt:
                return 0
    if args.cmd == "log":
        import subprocess
        from lignoforge.cgfit.monitor import current_mdrun
        path = wd / "log.txt"
        if args.md:
            cur = current_mdrun(wd)
            if not cur:
                print("no mdrun in progress")
                return 1
            path = wd / cur["run"]
        if not path.exists():
            print(f"no log yet: {path}")
            return 1
        return subprocess.call(["tail", "-n", str(args.n)] + (["-f"] if args.follow else []) + [str(path)])
    if args.cmd == "slurm":
        part = f"#SBATCH --partition={args.partition}\n" if args.partition else ""
        script = (f"#!/usr/bin/env bash\n#SBATCH --job-name=cgparam\n{part}"
                  f"#SBATCH --gres={args.gres}\n#SBATCH --cpus-per-task={args.cpus}\n"
                  f"#SBATCH --time={args.time}\n#SBATCH --output=slurm-%j.out\n"
                  f"export OMP_NUM_THREADS={args.cpus}\n"
                  f"lignoforge-cg-param run {wd.resolve()}\n")
        (wd / "submit.sh").write_text(script)
        print(f"wrote {wd / 'submit.sh'}  (sbatch it; edit mdrun_args in config.json for GPU runs)")
        return 0
    return 1


if __name__ == "__main__":
    sys.exit(main())
