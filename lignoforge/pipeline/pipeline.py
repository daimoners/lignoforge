"""
Top-down lignin modelling pipeline.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional

from lignoforge.io.schema import InputSchemaValidator
from lignoforge.priors.estimator import LigninPriorEstimator
from lignoforge.pipeline.translator import ParameterTranslator
from lignoforge.simulation.population import Simulation
from lignoforge.io.exporters import LigninExporter


@dataclass
class LigninResults:
    input_data: dict
    priors: dict
    simulation_kwargs: dict
    polymers: list
    output_dir: str
    artifacts: dict


class LigninPipeline:
    """
    End-to-end pipeline:
    input JSON -> schema validation -> prior estimation -> parameter translation
    -> kMC population simulation -> export.
    """

    def __init__(self, input_data: dict, output_dir: str | None = None):
        self.input_data = input_data
        self.output_dir = output_dir or os.path.join(os.getcwd(), "lignoforge_results")
        os.makedirs(self.output_dir, exist_ok=True)

        self.priors: Optional[dict] = None
        self.simulation_kwargs: Optional[dict] = None
        self.polymers: Optional[list] = None
        self.artifacts: dict = {}

    @classmethod
    def from_json(cls, json_path: str, schema_path: str | None = None, output_dir: str | None = None):
        validator = InputSchemaValidator(schema_path=schema_path)
        data = validator.validate_file(json_path)
        return cls(data, output_dir=output_dir)

    @classmethod
    def from_dict(cls, data: dict, schema_path: str | None = None, output_dir: str | None = None, validate: bool = True):
        """
        Create pipeline from a dictionary, optionally validating against schema.

        Parameters
        ----------
        data : dict
            Input specification dictionary
        schema_path : str, optional
            Path to JSON schema. If None, uses lignin_info_schema.json
        output_dir : str, optional
            Output directory. If None, uses lignoforge_results in current directory
        validate : bool, default=True
            Validate data against schema before accepting

        Returns
        -------
        LigninPipeline
            Initialized pipeline
        """
        if validate:
            validator = InputSchemaValidator(schema_path=schema_path)
            validator.validate_dict(data)
        return cls(data, output_dir=output_dir)

    def estimate_priors(self, random_seed: int | None = None) -> dict:
        estimator = LigninPriorEstimator(self.input_data, random_seed=random_seed)
        self.priors = estimator.run()
        return self.priors

    def translate_parameters(self) -> dict:
        if self.priors is None:
            self.estimate_priors()
        translator = ParameterTranslator(self.priors, self.input_data)
        self.simulation_kwargs = translator.to_simulation_kwargs()
        return self.simulation_kwargs

    def run(
        self,
        random_seed: int = 1,
        library_name: str = "lignoforge_library",
        simulation_overrides: Optional[dict] = None,
    ) -> LigninResults:
        if self.priors is None:
            self.estimate_priors(random_seed=random_seed)
        if self.simulation_kwargs is None:
            self.translate_parameters()

        sim_kwargs = dict(self.simulation_kwargs)
        if simulation_overrides:
            sim_kwargs.update(simulation_overrides)

        export_option_keys = {
            "generate_3d",
            "include_hydrogens",
            "optimize_3d",
            "max_uff_iterations",
            "n_workers",
        }
        export_options = {}
        for key in list(sim_kwargs.keys()):
            if key in export_option_keys:
                export_options[key] = sim_kwargs.pop(key)

        self.simulation_kwargs = sim_kwargs
        sim_kwargs.update({
            "seed_init": random_seed,
            "library_name": library_name,
            "results_name": "population",
            "save_path": self.output_dir,
        })

        sim = Simulation(**sim_kwargs)
        self.polymers = sim.run()

        exporter = LigninExporter(self.output_dir)
        self.artifacts = exporter.export_full_pipeline_bundle(
            input_data=self.input_data,
            priors=self.priors,
            sim_kwargs=self.simulation_kwargs,
            polymers=self.polymers,
            export_options=export_options,
        )

        return LigninResults(
            input_data=self.input_data,
            priors=self.priors,
            simulation_kwargs=self.simulation_kwargs,
            polymers=self.polymers,
            output_dir=self.output_dir,
            artifacts=self.artifacts,
        )

    def run_single_chain(
        self,
        random_seed: int = 1,
        library_name: str = "lignoforge_single_chain",
        simulation_overrides: Optional[dict] = None,
    ) -> LigninResults:
        """
        Convenience method to generate one model lignin chain.

        It forces n_population=1 and keeps all other parameters inferred from
        the priors unless explicitly overridden.
        """
        overrides = {
            "n_population": 1,
            "i_max": 40,
            "i_max_out": 40,
            "i_max_ring": 0,
        }
        if simulation_overrides:
            overrides.update(simulation_overrides)

        return self.run(
            random_seed=random_seed,
            library_name=library_name,
            simulation_overrides=overrides,
        )

    def export_single_pdb(
        self,
        filename: str = "single_chain.pdb",
        random_seed: int = 42,
        explicit_hydrogens: bool = True,
        optimize_3d: bool = True,
        max_uff_iterations: int = 200,
    ) -> str:
        if not self.polymers:
            raise ValueError("No polymers available. Run `run()` or `run_single_chain()` first.")
        exporter = LigninExporter(self.output_dir)
        return exporter.export_single_pdb(
            self.polymers[0],
            filename=filename,
            random_seed=random_seed,
            explicit_hydrogens=explicit_hydrogens,
            optimize_3d=optimize_3d,
            max_uff_iterations=max_uff_iterations,
        )

    def export_population_pdbs(
        self,
        folder: str = "pdb",
        basename: str = "polymer",
        random_seed: int = 42,
        explicit_hydrogens: bool = True,
        optimize_3d: bool = True,
        max_uff_iterations: int = 200,
    ) -> list[str]:
        if not self.polymers:
            raise ValueError("No polymers available. Run `run()` or `run_single_chain()` first.")
        exporter = LigninExporter(self.output_dir)
        return exporter.export_population_pdbs(
            self.polymers,
            folder=folder,
            basename=basename,
            random_seed=random_seed,
            explicit_hydrogens=explicit_hydrogens,
            optimize_3d=optimize_3d,
            max_uff_iterations=max_uff_iterations,
        )
