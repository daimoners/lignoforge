"""
Complete atomistic GROMACS workflow (OPLS-AA) for one lignin chain.

``write_atomistic_workflow`` writes a self-contained directory::

    <name>.top  <name>.gro      topology and coordinates (lignoforge.forcefield)
    em.mdp nvt.mdp npt.mdp prod.mdp
    run.sh                       box → solvate → EM → NVT → NPT → production
    workflow.json                all parameters used

The script only calls standard GROMACS tools and stops on the first error.
"""

from __future__ import annotations

import json
import os
import stat
from typing import Dict

from lignoforge.forcefield.topology import build_chain_topology
from lignoforge.md import mdp

__all__ = ["write_atomistic_workflow"]

_WATER_MODELS = {"tip3p": "oplsaa.ff/tip3p.itp", "spce": "oplsaa.ff/spce.itp"}


def _write(path: str, text: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def _run_script(name: str, solvent: str, box_d: float) -> str:
    head = "#!/usr/bin/env bash\nset -euo pipefail\nGMX=${GMX:-gmx}\nNT=${NT:-0}\n\n"
    mdrun = '$GMX mdrun -nt "$NT" -deffnm'
    if solvent == "none":
        return head + f"""\
# Vacuum / implicit-environment run: chain in a large periodic box.
$GMX editconf -f {name}.gro -o box.gro -bt cubic -d {box_d} -c
$GMX grompp -f em.mdp -c box.gro -p {name}.top -o em.tpr
{mdrun} em
$GMX grompp -f nvt.mdp -c em.gro -p {name}.top -o nvt.tpr
{mdrun} nvt
$GMX grompp -f prod.mdp -c nvt.gro -t nvt.cpt -p {name}.top -o prod.tpr
{mdrun} prod
"""
    return head + f"""\
$GMX editconf -f {name}.gro -o box.gro -bt cubic -d {box_d} -c
$GMX solvate -cp box.gro -cs spc216.gro -o solv.gro -p {name}.top
$GMX grompp -f em.mdp -c solv.gro -p {name}.top -o em.tpr
{mdrun} em
$GMX grompp -f nvt.mdp -c em.gro -p {name}.top -o nvt.tpr
{mdrun} nvt
$GMX grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p {name}.top -o npt.tpr
{mdrun} npt
$GMX grompp -f prod.mdp -c npt.gro -t npt.cpt -p {name}.top -o prod.tpr
{mdrun} prod
"""


def write_atomistic_workflow(
    chain: dict,
    output_dir: str,
    name: str = "lignin",
    temperature: float = 300.0,
    solvent: str = "tip3p",
    box_distance_nm: float = 1.2,
    nvt_ps: float = 100.0,
    npt_ps: float = 200.0,
    prod_ns: float = 10.0,
) -> Dict[str, str]:
    """
    Write the full atomistic workflow for one chain topology dict.

    Parameters
    ----------
    solvent
        ``"tip3p"``, ``"spce"`` or ``"none"`` (periodic vacuum
        box; NVT only).
    """
    if solvent != "none" and solvent not in _WATER_MODELS:
        raise ValueError(
            f"solvent must be 'none' or one of {sorted(_WATER_MODELS)}; got {solvent!r}"
        )
    topo = build_chain_topology(chain, name=name)
    includes = [] if solvent == "none" else [_WATER_MODELS[solvent]]

    paths = {
        "top": topo.write_top(os.path.join(output_dir, f"{name}.top"),
                              extra_includes=includes),
        "gro": topo.write_gro(os.path.join(output_dir, f"{name}.gro")),
        "em": _write(os.path.join(output_dir, "em.mdp"), mdp.em_mdp()),
        "nvt": _write(os.path.join(output_dir, "nvt.mdp"),
                      mdp.nvt_mdp(nvt_ps, temperature)),
        "prod": _write(os.path.join(output_dir, "prod.mdp"),
                       mdp.prod_mdp(prod_ns * 1000.0, temperature,
                                    npt=solvent != "none")),
    }
    if solvent != "none":
        paths["npt"] = _write(os.path.join(output_dir, "npt.mdp"),
                              mdp.npt_mdp(npt_ps, temperature))
    run = _write(os.path.join(output_dir, "run.sh"),
                 _run_script(name, solvent, box_distance_nm))
    os.chmod(run, os.stat(run).st_mode | stat.S_IXUSR)
    paths["run"] = run
    paths["config"] = _write(
        os.path.join(output_dir, "workflow.json"),
        json.dumps({
            "name": name, "forcefield": "OPLS-AA", "solvent": solvent,
            "temperature_K": temperature, "box_distance_nm": box_distance_nm,
            "nvt_ps": nvt_ps, "npt_ps": npt_ps if solvent != "none" else None,
            "prod_ns": prod_ns, "net_charge": topo.net_charge,
            "n_atoms": len(topo.atoms), "linkages": topo.linkages,
            "residues": topo.charge_report,
        }, indent=2),
    )
    return paths
