"""
LignoForge
==========

A top-down framework for automated generation of lignin structural models
from experimental and industrial input data.

Pipeline overview:
    1. Load experimental/industrial JSON input
    2. Estimate structural priors (S/G/H, linkage distribution, MW, branching)
    3. Translate priors to kMC simulation parameters
    4. Run kMC/MCMC polymer growth simulation
    5. Export population of atomistic polymer structures (SMILES, MOL, SDF, ...)

Example usage::

    from lignoforge.pipeline import LigninPipeline

    pipeline = LigninPipeline.from_json("examples/lignin_input_example.json")
    results  = pipeline.run()
    results.export("output/")
"""

__version__ = "0.1.0"
__author__  = "LignoForge contributors"

from lignoforge.pipeline.pipeline import LigninPipeline
from lignoforge.structure import MolecularStructureGenerator, PDBStructureWriter

__all__ = ["LigninPipeline", "MolecularStructureGenerator", "PDBStructureWriter"]
