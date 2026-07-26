"""The solver interface, and the registry that lets a backend live anywhere.

Two rules make the comparison harness trustworthy:

1. **A solver declares what it can do.** :class:`Capabilities` says whether the
   backend handles conical mounts, which polarizations it supports, and what
   its accuracy knob is called. The convergence harness reads that last field,
   which is why it exists from the very first solver even though scalar theory
   has no knob to turn.

2. **A solver refuses what it cannot do.** Silently approximating an
   unsupported configuration is the failure mode that makes cross-method plots
   lie -- a C-method quietly smoothing a vertical facet returns a plausible,
   wrong number. Raise :class:`UnsupportedConfiguration` instead.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib.metadata import entry_points
from typing import Protocol, runtime_checkable

from numpy.typing import ArrayLike

from ..illumination import Illumination, Polarization
from ..problem import Problem
from ..result import EfficiencyScan

__all__ = [
    "Capabilities",
    "Solver",
    "UnsupportedConfiguration",
    "register",
    "get_solver",
    "available_solvers",
]


class UnsupportedConfiguration(NotImplementedError):
    """A solver was asked for something outside its declared capabilities."""


@dataclass(frozen=True, slots=True)
class Capabilities:
    """What a backend can actually do, declared rather than discovered.

    Attributes
    ----------
    name
        Registry key, e.g. ``"scalar"``, ``"rcwa"``, ``"pcgrate:file"``.
    conical
        Handles off-plane mounts (``gamma != 90``). A backend that only does
        in-plane must say so; the off-plane X-ray case is the primary
        application here, and a solver silently treating it as in-plane would
        be wrong by a factor of ``sin(gamma)`` in the phase.
    polarizations
        Which polarizations the backend resolves. Scalar theory neglects
        polarization entirely, so it advertises only ``"unpolarized"``.
    accuracy_knob
        Name of the keyword controlling numerical accuracy -- ``"truncation"``
        for RCWA, ``"boundary_points"`` for the integral method. ``None`` for a
        closed-form method. The convergence harness sweeps this parameter.
    rigorous
        False for approximate theories. Marks results that should never be
        presented as reference values.
    handles_undercut
        True only for methods that parametrise the boundary rather than
        assuming a height function -- in practice, the integral method.
    """

    name: str
    conical: bool = True
    polarizations: tuple[Polarization, ...] = ("TE", "TM", "unpolarized")
    accuracy_knob: str | None = None
    rigorous: bool = True
    handles_undercut: bool = False

    def check(self, problem: Problem, illumination: Illumination) -> None:
        """Raise :class:`UnsupportedConfiguration` if this case is out of scope."""
        if not self.conical and not illumination.is_in_plane:
            raise UnsupportedConfiguration(
                f"{self.name} handles in-plane mounts only, but gamma is "
                f"{illumination.gamma_deg}deg. Treating an off-plane mount as "
                "in-plane would drop a factor of sin(gamma) from the phase."
            )
        if illumination.polarization not in self.polarizations:
            raise UnsupportedConfiguration(
                f"{self.name} does not resolve {illumination.polarization} "
                f"polarization; it supports {', '.join(self.polarizations)}"
            )
        if not self.handles_undercut and not problem.profile.is_single_valued():
            raise UnsupportedConfiguration(
                f"{self.name} cannot represent an undercut profile. Only a "
                "method that parametrises the boundary (integral) can."
            )


@runtime_checkable
class Solver(Protocol):
    """What every backend provides."""

    capabilities: Capabilities

    def solve(
        self,
        problem: Problem,
        illumination: Illumination,
        wavelengths: ArrayLike,
        **options,
    ) -> EfficiencyScan:
        """Efficiency over a wavelength scan at fixed geometry."""


_REGISTRY: dict[str, Solver] = {}


def register(solver: Solver) -> Solver:
    """Register a solver under its declared name. Usable as a decorator."""
    name = solver.capabilities.name
    if name in _REGISTRY and _REGISTRY[name] is not solver:
        raise ValueError(f"a different solver is already registered as {name!r}")
    _REGISTRY[name] = solver
    return solver


def _load_entry_points() -> None:
    """Pull in backends published by other distributions.

    A backend can live in its own package and still appear here, which is the
    point of keeping the problem spec free of solver fields.
    """
    for entry in entry_points(group="gratinglab.solvers"):
        if entry.name not in _REGISTRY:
            register(entry.load())


def get_solver(name: str) -> Solver:
    """Look up a registered solver by name."""
    if name not in _REGISTRY:
        _load_entry_points()
    if name not in _REGISTRY:
        known = ", ".join(sorted(_REGISTRY)) or "none"
        raise KeyError(f"no solver registered as {name!r}; available: {known}")
    return _REGISTRY[name]


def available_solvers() -> tuple[str, ...]:
    """Every registered solver name, sorted."""
    _load_entry_points()
    return tuple(sorted(_REGISTRY))
