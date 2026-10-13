"""
CG parametriser: derives a coarse-grained model of bulk lignin from atomistic
OPLS-AA simulations that it runs itself (see :mod:`lignoforge.cgfit.workflow`).
"""

from lignoforge.cgfit.model import SCHEMES, CGMolecule, CGParams, map_chain
from lignoforge.cgfit.workflow import DEFAULT_CONFIG, STAGES, CGParamJob, quick_config

__all__ = ["SCHEMES", "CGMolecule", "CGParams", "map_chain", "CGParamJob", "DEFAULT_CONFIG",
           "STAGES", "quick_config"]
