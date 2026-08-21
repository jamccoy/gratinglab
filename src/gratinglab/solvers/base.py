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

from dataclasses import dataclass
from importlib.metadata import entry_points
from typing import Callable, Protocol, runtime_checkable

from numpy.typing import ArrayLike

from ..illumination import Illumination, Polarization
from ..problem import Problem
from ..result import EfficiencyScan

__all__ = [
    "Capabilities",
    "Progress",
    "Solver",
    "SolveCancelled",
    "UnsupportedConfiguration",
    "register",
    "get_solver",
    "available_solvers",
]

#: ``(done, total)``, in units the solver chooses -- wavelengths, for every
#: backend so far. See :class:`Solver` for the full contract.
Progress = Callable[[int, int], None]


class UnsupportedConfiguration(NotImplementedError):
    """A solver was asked for something outside its declared capabilities."""


class SolveCancelled(Exception):
    """Raised **by a caller's progress callback** to abandon a solve.

    A Python thread cannot be killed and a NumPy-bound call has no check
    point, so a caller that wants to interrupt a solve has to be given a
    moment in which to do it. The progress callback is that moment: raise this
    from inside it and the solver unwinds.

    A solver must let it through. Wrapping the callback in ``except
    Exception`` would swallow it and quietly restore the old behaviour, in
    which cancelling only stops the *waiting*.

    Deriving from ``Exception`` rather than ``BaseException`` is deliberate,
    and it is the weaker of the two against accidental swallowing. The
    alternative slips past every ordinary handler -- including the broad one
    in ``gui/qt/worker.py`` that exists to keep a solver failure from aborting
    the process -- and an unhandled ``BaseException`` on a worker thread is a
    worse failure than the one being guarded against. The contract is stated
    and tested instead.
    """


@dataclass(frozen=True, slots=True)
class Capabilities:
    """What a backend can actually do, declared rather than discovered.

    Attributes
    ----------
    name
        Registry key, e.g. ``"scalar"``, ``"rcwa"``, ``"integral"``.
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
    reports_progress
        Accepts a ``progress`` keyword and honours the contract in
        :class:`Solver`. Declared rather than assumed because a backend that
        has not been updated -- including a third-party one arriving through
        the ``gratinglab.solvers`` entry point -- has an explicit signature
        that would raise ``TypeError`` on an unexpected keyword. A caller
        checks this before passing one.

        A solver declaring ``False`` is simply not cancellable. That is a fact
        about the backend, and a UI should say so rather than offer a Cancel
        that does nothing.
    """

    name: str
    conical: bool = True
    polarizations: tuple[Polarization, ...] = ("TE", "TM", "unpolarized")
    accuracy_knob: str | None = None
    rigorous: bool = True
    handles_undercut: bool = False
    reports_progress: bool = False

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
    """What every backend provides.

    **The progress contract**, for a backend declaring
    ``Capabilities.reports_progress``. ``solve`` takes a keyword-only
    ``progress: Progress | None = None`` and:

    1. calls it once with ``(0, total)`` **before doing any work**, then once
       after each unit completes, ending exactly at ``(total, total)``;
    2. keeps ``done`` monotone non-decreasing and ``total`` constant;
    3. lets any exception the callback raises **propagate unchanged** -- this
       is what makes :class:`SolveCancelled` work, and a solver that catches
       it silently removes the caller's only way out;
    4. costs nothing measurable when ``progress`` is ``None``.

    The leading ``(0, total)`` is not decoration. For RCWA a single wavelength
    can take a minute, so a cancellation point *before* the first one is the
    difference between stopping now and stopping in a minute.
    """

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
