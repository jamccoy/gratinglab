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
from .rcwa import RcwaSolver, rcwa
from .scalar import ScalarSolver, interference_factor, scalar


def solve(problem, illumination, wavelengths, method: str = "scalar", **options):
    """Efficiency over ``wavelengths`` by the solver named ``method``.

    A shorthand for ``get_solver(method).solve(...)``, for a script that
    wants one call. ``options`` are that solver's own settings, such as
    ``quadrature_points`` for ``"scalar"`` or ``truncation`` for ``"rcwa"``.
    Returns an :class:`~gratinglab.result.EfficiencyScan`.
    """
    return get_solver(method).solve(problem, illumination, wavelengths, **options)

__all__ = [
    "Capabilities",
    "Progress",
    "SolveCancelled",
    "Solver",
    "UnsupportedConfiguration",
    "available_solvers",
    "get_solver",
    "register",
    "solve",
    "CMethodSolver",
    "cmethod",
    "ScalarSolver",
    "IntegralSolver",
    "integral",
    "RcwaSolver",
    "rcwa",
    "interference_factor",
    "scalar",
]
