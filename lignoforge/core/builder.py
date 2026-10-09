"""
Library-level chain construction.

``resolve_spec`` turns a user specification (explicit fractions and/or an
experimental-input JSON) into fully explicit monomer and linkage distributions;
``grow_chain_exact`` grows a chain of exactly *n* monomers.  Both are shared by
the command line, the project API and the web interface.
"""

from __future__ import annotations

import warnings
from typing import Dict, Optional

import numpy as np

from lignoforge.core.rules import linkage_names, monomer_types

__all__ = ["resolve_spec", "grow_chain_exact", "normalise_fractions"]

_DEFAULT_INPUT = {
    "material_origin": {"biomass_type": "hardwood"},
    "extraction_process": {"process_type": "kraft"},
}


def normalise_fractions(values: Dict[str, float], keys) -> Dict[str, float]:
    """Return ``{key: fraction}`` over *keys*, normalised to sum 1."""
    v = np.array([max(0.0, float(values.get(k, 0.0))) for k in keys])
    total = v.sum()
    if total <= 0:
        raise ValueError(f"All fractions are zero for {list(keys)}")
    return {k: float(x) for k, x in zip(keys, v / total)}


def resolve_spec(spec: dict) -> dict:
    """
    Resolve a chain specification.

    Recognised keys (all optional): ``n_monomers``, ``monomers`` (``{"H","G","S"}``
    fractions), ``linkages`` (``{linkage: fraction}``), ``branching``, ``seed``,
    ``experimental_input`` (LignoForge input JSON used to estimate priors for
    anything not given explicitly).

    Returns a dict with explicit, normalised ``monomers``/``linkages`` and
    ``warnings`` (list of strings).
    """
    from lignoforge.core.polymer import ENABLE_BETA_1_LINKAGE

    seed = int(spec.get("seed", 42))
    notes = []

    base_m = {"H": 0.2, "G": 0.5, "S": 0.3}
    base_l = {k: 1.0 for k in linkage_names}
    expected_dp: Optional[int] = None

    exp = spec.get("experimental_input")
    if exp is not None:
        from lignoforge.io.schema import InputSchemaValidator
        from lignoforge.pipeline.translator import ParameterTranslator
        from lignoforge.priors.estimator import LigninPriorEstimator

        InputSchemaValidator().validate_dict(exp)     # raises on invalid input
        data = exp
        priors = LigninPriorEstimator(data, random_seed=seed).run()
        kw = ParameterTranslator(priors, data).to_simulation_kwargs()
        base_m = dict(zip(monomer_types, kw["monomer_distribution_input"]))
        base_l = dict(zip(linkage_names, kw["linkage_distribution_input"]))
        expected_dp = int(kw.get("expected_size", 20))

    monomers = dict(base_m)
    monomers.update({k: float(v) for k, v in (spec.get("monomers") or {}).items()
                     if v is not None})
    unknown = set(monomers) - set(monomer_types)
    if unknown:
        raise ValueError(f"Unknown monomer type(s): {sorted(unknown)}")

    explicit_l = {k: float(v) for k, v in (spec.get("linkages") or {}).items()
                  if v is not None}
    unknown = set(explicit_l) - set(linkage_names)
    if unknown:
        raise ValueError(f"Unknown linkage(s): {sorted(unknown)}")
    # Explicit linkage fractions replace the estimated ones wholesale; those
    # not mentioned are 0 (a partial dict is a full specification).
    linkages = {k: explicit_l.get(k, 0.0) for k in linkage_names} if explicit_l \
        else dict(base_l)

    if not ENABLE_BETA_1_LINKAGE and linkages.get("beta-1", 0.0) > 0:
        if explicit_l.get("beta-1", 0.0) > 0:
            notes.append("beta-1 is not supported by the structure builder; "
                         "its fraction was set to 0.")
        linkages["beta-1"] = 0.0

    n = spec.get("n_monomers", expected_dp or 20)
    return {
        "n_monomers": int(n),
        "monomers": normalise_fractions(monomers, monomer_types),
        "linkages": normalise_fractions(linkages, linkage_names),
        "branching": float(spec.get("branching") or 0.0),
        "seed": seed,
        "warnings": notes,
    }


def grow_chain_exact(
    n: int,
    monomer_dist: list,
    linkage_dist: list,
    seed: int,
    branching: Optional[float],
    verbose: bool = False,
):
    """
    Grow a chain of exactly *n* monomers.

    ``branching`` of 0 / ``None`` means strictly linear growth (only terminal
    monomers accept a new unit).  If no compatible linkage is found after 20
    attempts the chain is returned shorter and a ``RuntimeWarning`` is raised.
    """
    from lignoforge.core.monomer import Monomer
    from lignoforge.core.polymer import Polymer
    from lignoforge.core.utils import (
        generate_random_branching_state,
        generate_random_monomer,
        set_random_state,
    )

    rstate = set_random_state(seed)
    m_dist = np.asarray(monomer_dist)
    l_dist = np.asarray(linkage_dist)

    polymer = Polymer(Monomer(generate_random_monomer(m_dist, rstate)), verbose=verbose)
    p_branch = branching or 0.0
    for step in range(n - 1):
        for _attempt in range(20):
            # Re-drawn on every attempt so a branch request that cannot be
            # satisfied (no interior monomer with a free site) does not
            # exhaust the attempts.
            b_state = (generate_random_branching_state(p_branch, rstate)
                       if p_branch > 0.0 else False)
            if polymer.add_random_monomer(
                monomer_distribution=m_dist,
                linkage_distribution=l_dist,
                branching_state=b_state,
                random_state=rstate,
            ):
                break
        else:
            warnings.warn(
                f"Chain growth stalled at {step + 1} of {n} monomers: no "
                f"compatible linkage found after 20 attempts "
                f"(check linkage/monomer fractions).",
                RuntimeWarning,
            )
            break
    return polymer
