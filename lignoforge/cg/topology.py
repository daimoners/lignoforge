"""
GROMACS topology and workflow for the monomer-bead CG model.

* one moleculetype per chain, one bead per monomer (mass = monomer mass);
* harmonic bonds per linkage and harmonic angles per linkage pair, from a
  :class:`~lignoforge.cg.parameters.CGParameters` set;
* Lennard-Jones beads (combination rule 2), uncharged, with 1-2 and 1-3
  exclusions;
* stochastic dynamics (``integrator = sd``), i.e. implicit solvent.

The workflow script can place several copies of one or more chains in a cubic
box (``gmx insert-molecules``) at a target density, minimise, equilibrate and
run.  Because the beads attract each other and there is no explicit solvent,
a single chain in vacuum collapses unless ``epsilon`` is lowered; choose the
parameters for the solvent quality of interest.
"""

from __future__ import annotations

import json
import os
import stat
from typing import Dict, List, Optional, Sequence

import numpy as np

from lignoforge.cg.mapping import CGChain, map_topology
from lignoforge.cg.parameters import CGParameters, default_parameters
from lignoforge.forcefield.topology import build_chain_topology

__all__ = ["write_cg_system", "cg_chain_from_atomistic"]

_AVOGADRO = 6.02214076e23


def cg_chain_from_atomistic(chain: dict, name: str) -> CGChain:
    """Map an atomistic chain topology dict onto its CG chain."""
    return map_topology(build_chain_topology(chain, name=name))


def _moleculetype(cg: CGChain, p: CGParameters) -> List[str]:
    types = {b.index: b.mtype for b in cg.beads}
    L = ["[ moleculetype ]", f"{cg.name}  2", "", "[ atoms ]",
         ";  nr  type  resnr residue  atom  cgnr  charge   mass"]
    for b in cg.beads:
        L.append(f"{b.index:5d}  CG{b.mtype}  {b.index:5d}  {b.resname:<5s}  "
                 f"BB{b.mtype}  {b.index:5d}  0.000  {b.mass:9.3f}")
    L += ["", "[ bonds ]", ";  ai   aj  funct  b0(nm)  kb(kJ/mol/nm2)"]
    for i, j, link in cg.bonds:
        q = p.bond(link, types[i], types[j])
        L.append(f"{i:5d}{j:5d}  1  {q['r0']:8.4f}  {q['k']:9.1f}   ; {link}")
    L += ["", "[ angles ]", ";  ai   aj   ak  funct  theta0(deg)  k(kJ/mol/rad2)"]
    for i, j, k, l1, l2 in cg.angles:
        q = p.angle(l1, l2)
        L.append(f"{i:5d}{j:5d}{k:5d}  1  {q['theta0']:8.2f}  {q['k']:8.2f}   ; {l1}|{l2}")
    L.append("")
    return L


def _header(p: CGParameters, used_types: Sequence[str]) -> List[str]:
    L = [f"; LignoForge coarse-grained topology (1 bead / monomer)",
         f"; {p.provenance}", "", "[ defaults ]",
         "; nbfunc  comb-rule  gen-pairs  fudgeLJ  fudgeQQ", "  1  2  no  1.0  1.0", "",
         "[ atomtypes ]", "; name  at.num  mass  charge  ptype  sigma(nm)  epsilon(kJ/mol)"]
    for t in sorted(used_types):
        b = p.beads[t]
        L.append(f"CG{t}  0  0.000  0.000  A  {b['sigma']:.4f}  {b['epsilon']:.4f}")
    L.append("")
    return L


_EM = """\
; CG energy minimisation (LignoForge)
integrator   = steep
emtol        = 100.0
emstep       = 0.005
nsteps       = 20000
cutoff-scheme = Verlet
coulombtype  = cut-off
rcoulomb     = {rc}
rvdw         = {rc}
vdw-modifier = potential-shift
pbc          = xyz
"""

_SD = """\
; CG stochastic dynamics, implicit solvent (LignoForge)
integrator   = sd
dt           = {dt}
nsteps       = {nsteps}
tc-grps      = System
tau-t        = {tau_t}
ref-t        = {T}
nstxout-compressed = {nout}
nstenergy    = {nout}
nstlog       = {nout}
cutoff-scheme = Verlet
nstlist      = 20
coulombtype  = cut-off
rcoulomb     = {rc}
rvdw         = {rc}
vdw-modifier = potential-shift
pbc          = xyz
gen-vel      = yes
gen-temp     = {T}
"""


def _write(path: str, text: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")
    return path


def _gro(cg: CGChain, xyz_nm: np.ndarray, box: Sequence[float]) -> str:
    L = [cg.name, f"{len(cg.beads):5d}"]
    for b, x in zip(cg.beads, xyz_nm):
        L.append(f"{b.index:5d}{b.resname:<5s}{('BB' + b.mtype):>5s}{b.index:5d}"
                 f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}")
    L.append(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}")
    return "\n".join(L)


def write_cg_system(
    chains: Sequence[dict],
    output_dir: str,
    name: str = "cg",
    parameters: Optional[CGParameters] = None,
    copies: int = 1,
    density_g_cm3: float = 0.3,
    temperature: float = 300.0,
    run_ns: float = 100.0,
    dt_ps: float = 0.01,
    cutoff_nm: float = 1.4,
) -> Dict[str, str]:
    """
    Write a CG system: ``<name>.top``, one ``.gro`` per chain, ``em.mdp``,
    ``md.mdp``, ``run.sh``, ``cg_parameters.json``.

    Parameters
    ----------
    chains
        Atomistic chain topology dicts; each becomes its own moleculetype.
    copies
        Number of copies of *each* chain in the box.
    density_g_cm3
        Target bead mass density used to size the cubic box.
    """
    p = parameters or default_parameters()
    cgs = [cg_chain_from_atomistic(c, f"{name}_{i}") for i, c in enumerate(chains)]
    used = sorted({b.mtype for cg in cgs for b in cg.beads})

    top = _header(p, used)
    for cg in cgs:
        top += _moleculetype(cg, p)
    top += ["[ system ]", name, "", "[ molecules ]"]
    top += [f"{cg.name}  {copies}" for cg in cgs]
    paths = {"top": _write(os.path.join(output_dir, f"{name}.top"), "\n".join(top))}

    total_mass = copies * sum(cg.mass_array().sum() for cg in cgs)        # g/mol
    volume_nm3 = total_mass / _AVOGADRO / density_g_cm3 * 1e21           # cm3 → nm3
    longest = max(float(np.ptp(cg.ref_xyz, axis=0).max()) for cg in cgs)
    box = max(volume_nm3 ** (1 / 3), longest + 2 * cutoff_nm, 2 * cutoff_nm + 0.5)

    for cg in cgs:
        xyz = cg.ref_xyz - cg.ref_xyz.mean(axis=0) + box / 2
        paths[f"gro_{cg.name}"] = _write(
            os.path.join(output_dir, f"{cg.name}.gro"), _gro(cg, xyz, [box] * 3))

    nsteps = int(round(run_ns * 1000.0 / dt_ps))
    nout = max(1, nsteps // 1000)
    paths["em"] = _write(os.path.join(output_dir, "em.mdp"), _EM.format(rc=cutoff_nm))
    paths["md"] = _write(os.path.join(output_dir, "md.mdp"), _SD.format(
        dt=dt_ps, nsteps=nsteps, tau_t=1.0, T=temperature, nout=nout, rc=cutoff_nm))

    if len(cgs) == 1 and copies == 1:
        build = f"cp {cgs[0].name}.gro system.gro\n"
    else:
        build = "".join(
            ("$GMX insert-molecules " +
             (f"-f system.gro " if i else f"-box {box:.4f} {box:.4f} {box:.4f} ") +
             f"-ci {cg.name}.gro -nmol {copies} -o system.gro -rot xyz\n")
            for i, cg in enumerate(cgs)
        )
    run = (
        "#!/usr/bin/env bash\nset -euo pipefail\nGMX=${GMX:-gmx}\nNT=${NT:-0}\n\n"
        + build +
        f"$GMX grompp -f em.mdp -c system.gro -p {name}.top -o em.tpr\n"
        '$GMX mdrun -nt "$NT" -deffnm em\n'
        f"$GMX grompp -f md.mdp -c em.gro -p {name}.top -o md.tpr\n"
        '$GMX mdrun -nt "$NT" -deffnm md\n'
    )
    paths["run"] = _write(os.path.join(output_dir, "run.sh"), run)
    os.chmod(paths["run"], os.stat(paths["run"]).st_mode | stat.S_IXUSR)

    p.save(os.path.join(output_dir, "cg_parameters.json"))
    paths["parameters"] = os.path.join(output_dir, "cg_parameters.json")
    _write(os.path.join(output_dir, "workflow.json"), json.dumps({
        "model": "monomer-bead CG (1 bead/monomer)", "chains": [c.name for c in cgs],
        "copies": copies, "box_nm": box, "density_g_cm3": density_g_cm3,
        "temperature_K": temperature, "run_ns": run_ns, "dt_ps": dt_ps,
        "cutoff_nm": cutoff_nm, "parameters_provenance": p.provenance,
    }, indent=2))
    return paths
