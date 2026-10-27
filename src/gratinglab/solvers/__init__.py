"""Solver backends.

Importing this package registers the built-in solvers. Third-party backends
register through the ``gratinglab.solvers`` entry point.
"""

from .base import (
    Capabilities,
    Progress,
    SolveCancelled,
    Solver,
    UnsupportedConfiguration,
    available_solvers,
    get_solver,
    register,
)
from .cmethod import CMethodSolver, cmethod
from .integral import IntegralSolver, integral
from .scalar import ScalarSolver, interference_factor, scalar

__all__ = [
    "Capabilities",
    "Progress",
    "SolveCancelled",
    "Solver",
    "UnsupportedConfiguration",
    "available_solvers",
    "get_solver",
    "register",
    "CMethodSolver",
    "cmethod",
    "ScalarSolver",
    "IntegralSolver",
    "integral",
    "interference_factor",
    "scalar",
]
