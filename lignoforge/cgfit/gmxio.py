"""
GROMACS input writers and a thin command runner for the CG parametriser.

Atomistic bulk: several distinct chains, each with ``copies`` copies, packed
with ``gmx insert-molecules``.  CG system: one moleculetype per chain with the
bonded terms and per-pair Lennard-Jones (``[ nonbond_params ]``) of a
:class:`~lignoforge.cgfit.model.CGParams`.
"""

from __future__ import annotations

import itertools
import math
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

from lignoforge.cgfit.model import CGMolecule, CGParams, all_type_pairs, pair_key
from lignoforge.forcefield.topology import EXTRA_DIHEDRALTYPES, ChainTopology

AVOGADRO = 6.02214076e23


class CommandError(RuntimeError):
    pass


class Gmx:
    """
    Run GROMACS commands in a working directory.

    ``mdrun_args`` (e.g. ``["-ntmpi", "1", "-ntomp", "8", "-nb", "gpu"]``) are
    appended to every ``mdrun``.  With ``dry_run=True`` commands are only
    appended to ``commands.log`` (useful to review or to hand over to a
    scheduler).
    """

    def __init__(self, gmx: Optional[str] = None, mdrun_args: Sequence[str] = (),
                 dry_run: bool = False):
        self.gmx = gmx or os.environ.get("GMX") or shutil.which("gmx") or "gmx"
        self.mdrun_args = list(mdrun_args)
        self.dry_run = dry_run

    def run(self, args: Sequence[str], cwd: Path, stdin: Optional[str] = None,
            log: Optional[str] = None, extra: Sequence[str] = ()) -> str:
        cmd = [self.gmx, *args, *extra]
        cwd = Path(cwd)
        cwd.mkdir(parents=True, exist_ok=True)
        line = " ".join(shlex.quote(c) for c in cmd)
        with open(cwd / "commands.log", "a") as fh:
            fh.write(line + "\n")
        if self.dry_run:
            return ""
        r = subprocess.run(cmd, cwd=cwd, input=stdin, capture_output=True, text=True)
        if log:
            (cwd / log).write_text(r.stdout + "\n" + r.stderr)
        if r.returncode != 0:
            raise CommandError(f"`{line}` failed in {cwd}:\n{r.stderr[-1500:]}")
        return r.stdout + r.stderr

    def mdrun(self, deffnm: str, cwd: Path, extra: Sequence[str] = ()) -> str:
        args = [*self.mdrun_args, *extra]
        try:
            return self.run(["mdrun", "-deffnm", deffnm], cwd, log=f"{deffnm}.mdrun.log", extra=args)
        except CommandError as e:
            # small systems cannot be split over many thread-MPI ranks: retry on one rank
            if "domain decomposition" in str(e) and "-ntmpi" not in args:
                return self.run(["mdrun", "-deffnm", deffnm], cwd, log=f"{deffnm}.mdrun.log",
                                extra=[*args, "-ntmpi", "1"])
            raise

    def grompp(self, mdp: str, gro: str, top: str, out: str, cwd: Path,
               maxwarn: int = 0, extra: Sequence[str] = ()) -> str:
        return self.run(["grompp", "-f", mdp, "-c", gro, "-p", top, "-o", out,
                         "-maxwarn", str(maxwarn)], cwd, log=f"{out}.grompp.log", extra=extra)


def write_text(path, text: str) -> str:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text if text.endswith("\n") else text + "\n")
    return str(path)


def read_gro_box(path) -> List[float]:
    last = Path(path).read_text().strip().splitlines()[-1].split()
    return [float(x) for x in last[:3]]


def read_gro_xyz(path) -> np.ndarray:
    lines = Path(path).read_text().splitlines()
    n = int(lines[1])
    return np.array([[float(l[20 + 8 * d:28 + 8 * d]) for d in range(3)]
                     for l in lines[2:2 + n]])


# ── atomistic bulk ────────────────────────────────────────────────────────────

def write_atomistic_top(topos: Sequence[ChainTopology], copies: Sequence[int], path) -> str:
    L = ["; atomistic lignin bulk (LignoForge cgfit)", '#include "oplsaa.ff/forcefield.itp"', "",
         EXTRA_DIHEDRALTYPES]
    for t in topos:
        L += t.moleculetype_lines() + [""]
    L += ["[ system ]", "lignin bulk", "", "[ molecules ]"]
    L += [f"{t.name}  {c}" for t, c in zip(topos, copies)]
    return write_text(path, "\n".join(L))


def box_edge_for_density(mass_gmol: float, density_g_cm3: float) -> float:
    """Edge (nm) of the cubic box holding *mass_gmol* at the given density."""
    return (mass_gmol / AVOGADRO / density_g_cm3 * 1e21) ** (1 / 3)


# ── atomistic mdp ─────────────────────────────────────────────────────────────

def atomistic_mdp(kind: str, time_ps: float, temperature: float, dt: float = 0.002,
                  pressure: float = 1.0, out_ps: float = 10.0, seed: int = -1,
                  continuation: bool = True) -> str:
    n = max(1, int(round(time_ps / dt)))
    nout = max(1, int(round(out_ps / dt)))
    ncalc = math.gcd(nout, 100)
    common = f"""\
cutoff-scheme = Verlet
nstlist       = 20
rcoulomb      = 1.0
rvdw          = 1.0
coulombtype   = PME
fourierspacing = 0.12
DispCorr      = EnerPres
pbc           = xyz
"""
    if kind == "em":
        return f"integrator = steep\nemtol = 500\nemstep = 0.001\nnsteps = 20000\n{common}"
    pc = ("pcoupl = no\n" if kind == "nvt" else
          f"pcoupl = C-rescale\npcoupltype = isotropic\ntau-p = 5.0\nref-p = {pressure}\n"
          "compressibility = 4.5e-5\nrefcoord-scaling = com\n")
    return f"""\
integrator    = md
dt            = {dt}
nsteps        = {n}
nstcalcenergy = {ncalc}
nstenergy     = {nout}
nstlog        = {nout}
nstxout-compressed = {nout}
continuation  = {"yes" if continuation else "no"}
constraints   = h-bonds
constraint-algorithm = lincs
tcoupl        = v-rescale
tc-grps       = System
tau-t         = 0.5
ref-t         = {temperature}
{pc}gen-vel       = {"no" if continuation else "yes"}
gen-temp      = {temperature}
gen-seed      = {seed}
{common}"""


# ── CG system ─────────────────────────────────────────────────────────────────

def _lj_c6c12(sigma: float, eps: float):
    return 4 * eps * sigma ** 6, 4 * eps * sigma ** 12


def cg_top_text(mols: Sequence[CGMolecule], copies: Sequence[int], params: CGParams,
                title: str = "cg") -> str:
    types = sorted({t for m in mols for t in m.types})
    L = [f"; LignoForge cgfit CG topology, scheme={params.scheme}, T={params.temperature} K",
         "", "[ defaults ]", "; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ",
         "  1  1  no  1.0  1.0", "", "[ atomtypes ]",
         "; name  at.num  mass  charge  ptype  c6  c12"]
    for t in types:
        L.append(f"{t}  0  0.000  0.000  A  0.0  0.0")
    L += ["", "[ nonbond_params ]", "; i  j  funct  c6  c12"]
    for pk in all_type_pairs(types):
        a, b = pk.split("-")
        q = params.nonbonded.get(pk)
        if q is None:
            raise KeyError(f"No nonbonded parameters for bead-type pair {pk}")
        c6, c12 = _lj_c6c12(q["sigma"], q["epsilon"])
        L.append(f"{a}  {b}  1  {c6:.8e}  {c12:.8e}")
    L.append("")
    for m in mols:
        L += ["[ moleculetype ]", f"{m.name}  1", "", "[ atoms ]",
              ";  nr  type  resnr  residue  atom  cgnr  charge  mass"]
        for b in m.beads:
            L.append(f"{b.index:5d}  {b.type:<4s} {b.resid:5d}  LIG  {b.type:<4s} {b.index:5d}  0.0  {b.mass:9.3f}")
        L += ["", "[ bonds ]", ";  ai  aj  funct  b0  kb"]
        for i, j, key in m.bonds:
            q = params.bond(key)
            L.append(f"{i:5d}{j:5d}  1  {q['r0']:8.4f}  {q['k']:9.1f}   ; {key}")
        if m.angles:
            L += ["", "[ angles ]", ";  ai  aj  ak  funct  theta0  k"]
            for i, j, k, spec, gen in m.angles:
                q = params.angle(spec, gen)
                L.append(f"{i:5d}{j:5d}{k:5d}  1  {q['theta0']:8.2f}  {q['k']:8.2f}   ; {spec}")
        L += ["", "[ exclusions ]"]
        L += [f"{i:5d}{j:5d}" for i, j in m.exclusions]
        L.append("")
    L += ["[ system ]", title, "", "[ molecules ]"]
    L += [f"{m.name}  {c}" for m, c in zip(mols, copies)]
    return "\n".join(L)


def cg_gro_text(types: Sequence[str], xyz_nm: np.ndarray, box: Sequence[float],
                title: str = "cg") -> str:
    L = [title, f"{len(types):5d}"]
    for n, (t, x) in enumerate(zip(types, xyz_nm), start=1):
        L.append(f"{n % 100000:5d}{'LIG':<5s}{t:>5s}{n % 100000:5d}"
                 f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}")
    L.append(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}")
    return "\n".join(L)


def cg_mdp(kind: str, time_ps: float, temperature: float, cutoff: float, dt: float = 0.005,
           out_ps: float = 5.0, tau_t: float = 1.0, pressure: float = 1.0,
           seed: int = -1, continuation: bool = False) -> str:
    common = f"""\
cutoff-scheme = Verlet
nstlist       = 20
verlet-buffer-tolerance = 0.005
coulombtype   = cut-off
rcoulomb      = {cutoff}
vdwtype       = cut-off
vdw-modifier  = potential-shift
rvdw          = {cutoff}
DispCorr      = no
pbc           = xyz
"""
    if kind == "em":
        return f"integrator = steep\nemtol = 50\nemstep = 0.001\nnsteps = 10000\n{common}"
    n = max(1, int(round(time_ps / dt)))
    nout = max(1, int(round(out_ps / dt)))
    pc = ("pcoupl = no\n" if kind == "nvt" else
          f"pcoupl = C-rescale\npcoupltype = isotropic\ntau-p = 5.0\nref-p = {pressure}\n"
          "compressibility = 4.5e-5\n")
    return f"""\
integrator    = sd
dt            = {dt}
nsteps        = {n}
nstcalcenergy = {math.gcd(nout, 100)}
nstenergy     = {nout}
nstlog        = {nout}
nstxout-compressed = {nout}
tc-grps       = System
tau-t         = {tau_t}
ref-t         = {temperature}
continuation  = {"yes" if continuation else "no"}
gen-vel       = {"no" if continuation else "yes"}
gen-temp      = {temperature}
gen-seed      = {seed}
ld-seed       = {seed}
{pc}{common}"""


_PRESSURE = re.compile(r"^Pressure\s+([-\d.eE+]+)\s+([-\d.eE+]+)", re.M)
_DENSITY = re.compile(r"^Density\s+([-\d.eE+]+)\s+([-\d.eE+]+)", re.M)


def energy_average(gmx: Gmx, edr: str, term: str, cwd: Path, begin_ps: float = 0.0) -> Optional[float]:
    """Time average of one ``gmx energy`` term after *begin_ps* (None if dry run / absent)."""
    if gmx.dry_run:
        return None
    out = gmx.run(["energy", "-f", edr, "-b", str(begin_ps)], cwd, stdin=f"{term}\n0\n")
    rx = {"Pressure": _PRESSURE, "Density": _DENSITY}.get(term)
    m = rx.search(out) if rx else None
    if not m:
        return None
    v = float(m.group(1))
    return v / 1000.0 if term == "Density" else v          # kg/m3 → g/cm3
