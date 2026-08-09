r"""Demonstrating convergence, or refusing to claim it.

Every :class:`~gratinglab.result.Provenance` this project has ever produced
carries ``converged=None`` -- "not checked". That is honest but useless: it
means no number here has been shown correct, only computed. This module is
what turns that field into evidence.

The idea is the ordinary one. A solver declares its accuracy knob
(:attr:`~gratinglab.solvers.base.Capabilities.accuracy_knob` -- ``truncation``
for RCWA, ``boundary_points`` for the integral method, ``quadrature_points``
for scalar quadrature). Sweep it upward, watch the answer stop moving, and
report where it stopped.

Why a *plateau*, not the first pair that agrees
===============================================

The tempting rule -- "converged once two consecutive refinements agree within
tolerance" -- is not safe, and this project has the measurement to prove it
rather than the principle. Scalar quadrature on the reference blazed grating,
maximum change in efficiency per doubling:

======  ==========  ===================
``n``   change      vs. previous change
======  ==========  ===================
512     2.265e-04   0.44x
1024    9.645e-06   **23x smaller**
2048    2.426e-05   **2.5x LARGER**
4096    4.835e-06   5.0x smaller
8192    5.964e-07   8.1x smaller
======  ==========  ===================

At a tolerance of 1e-5 the naive rule stops at ``n=512`` on the strength of
that 9.6e-6 -- and the very next refinement then moves the answer by 2.4e-5,
more than the tolerance just certified. The cause is geometric, not numerical
noise: the blazed profile has a slope kink at :math:`t = 0.833`, and doubling
``n`` changes how near a quadrature node falls to it, so the error is not
monotone in ``n``. Any method with a corner, a resonance, or an anomaly can do
this.

So agreement must persist. :data:`DEFAULT_PLATEAU` consecutive differences
must all fall below tolerance -- three knob values in a row that agree -- and
the accidental single agreement above is rejected by exactly one more solve.
Raise ``plateau`` for a method known to be worse behaved; there is no way to
lower it below 1, because a single agreement is the thing this exists to
distrust.

What is compared
================

Absolute difference in efficiency, maximised over every wavelength and order.
Absolute rather than relative because efficiencies live in :math:`[0, 1]` and
a relative error on an order carrying 1e-12 of the power is noise being
promoted to a headline. The order grid and the propagating mask must match
between rungs; if they do not, this raises rather than comparing misaligned
arrays, which is the failure mode the whole project is arranged against.

What it does *not* do
=====================

Prove correctness. A converged wrong model is still wrong -- scalar theory
converges beautifully to an answer that neglects polarization entirely. This
measures numerical self-consistency, and :mod:`gratinglab.checks` measures the
physics. They are different questions and both are needed.
"""

from __future__ import annotations

import dataclasses
import inspect
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from .illumination import Illumination
from .problem import Problem
from .result import EfficiencyScan
from .solvers.base import Progress, UnsupportedConfiguration

__all__ = [
    "ConvergenceReport",
    "DEFAULT_PLATEAU",
    "DEFAULT_TOLERANCE",
    "check_convergence",
    "doubling_ladder",
]

#: Consecutive agreements required. Two differences below tolerance means
#: three knob values in a row agree. See the module docstring for the measured
#: case that makes one insufficient.
DEFAULT_PLATEAU = 2

#: Comfortably below any experimental precision this would be compared against,
#: and reachable: the reference blazed case gets there by ``n = 4096``.
DEFAULT_TOLERANCE = 1e-6

#: How many rungs below the solver's own default :func:`doubling_ladder` starts,
#: when no explicit ladder is given. Starting *below* the default is the point
#: -- the actionable output is how coarse a setting would have done, which for
#: RCWA at N ~ 200-400 is the difference between a run and an afternoon.
_RUNGS_BELOW_DEFAULT = 3

#: Total rungs in the default ladder. Generous because the sweep stops as soon
#: as it finds a plateau, so a well-behaved case never pays for the tail.
_DEFAULT_RUNGS = 10


@dataclass(frozen=True, slots=True)
class ConvergenceReport:
    """What a sweep found, including when it found nothing.

    Attributes
    ----------
    knob
        Name of the swept keyword, from the solver's declared capabilities.
    values
        Knob values actually solved, ascending.
    differences
        ``differences[i]`` is the largest change in efficiency between
        ``values[i]`` and ``values[i + 1]``. One shorter than ``values``.
    tolerance, plateau
        The criterion applied.
    converged_at
        Coarsest knob value demonstrated adequate, or ``None`` if the sweep
        ended without a plateau. **This is the actionable number**: it is what
        a production run should use, not the finest value tried.
    scan
        The finest scan computed, with its provenance stamped -- ``converged``
        set to a real boolean at last, and ``notes["convergence"]`` carrying
        this evidence so a result can defend itself away from this object.
    skipped
        ``(value, reason)`` for rungs the solver refused, e.g. a quadrature
        too coarse to resolve the highest order. Recorded, not silently
        dropped: a ladder that mostly did not run is a different situation
        from one that ran and disagreed.
    """

    knob: str
    values: tuple[int, ...]
    differences: tuple[float, ...]
    tolerance: float
    plateau: int
    converged_at: int | None
    scan: EfficiencyScan
    skipped: tuple[tuple[int, str], ...] = ()

    @property
    def passed(self) -> bool:
        """True only if a plateau was actually found."""
        return self.converged_at is not None

    @property
    def final_difference(self) -> float | None:
        """The last change measured -- the tightest bound the sweep supports."""
        return self.differences[-1] if self.differences else None

    def as_notes(self) -> dict[str, object]:
        """The evidence, in a form that survives in ``Provenance.notes``."""
        return {
            "knob": self.knob,
            "values": list(self.values),
            "differences": [float(d) for d in self.differences],
            "tolerance": self.tolerance,
            "plateau": self.plateau,
            "converged_at": self.converged_at,
        }

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        if self.passed:
            return (
                f"convergence: demonstrated at {self.knob}={self.converged_at} "
                f"({self.plateau} consecutive changes below {self.tolerance:g}; "
                f"swept to {self.values[-1]})"
            )
        last = self.final_difference
        detail = f"best change {last:.3g}" if last is not None else "nothing to compare"
        return (
            f"convergence: NOT demonstrated -- swept {self.knob} to "
            f"{self.values[-1] if self.values else '?'}, {detail}, "
            f"tolerance {self.tolerance:g}"
        )


def doubling_ladder(start: int, rungs: int) -> tuple[int, ...]:
    """``start``, ``2*start``, ``4*start``, ... -- ``rungs`` values.

    Doubling because discretisation error typically falls geometrically, so a
    linear ladder spends most of its solves learning nothing while a geometric
    one brackets the plateau in a handful.
    """
    if start < 1:
        raise ValueError(f"start must be positive, got {start}")
    if rungs < 2:
        raise ValueError(f"a sweep needs at least two rungs to compare, got {rungs}")
    return tuple(start * 2**i for i in range(rungs))


def default_ladder(solver) -> tuple[int, ...]:
    """A ladder straddling the solver's own default for its accuracy knob.

    Read from the signature rather than hardcoded, so a backend that changes
    its default does not silently leave this sweeping the wrong range.
    """
    knob = _knob_of(solver)
    default = inspect.signature(solver.solve).parameters[knob].default
    if not isinstance(default, int) or default < 1:
        raise ValueError(
            f"{solver.capabilities.name}.solve has no usable default for "
            f"{knob!r} ({default!r}), so a ladder cannot be inferred -- pass "
            "`values` explicitly"
        )
    return doubling_ladder(max(1, default // 2**_RUNGS_BELOW_DEFAULT), _DEFAULT_RUNGS)


def check_convergence(
    solver,
    problem: Problem,
    illumination: Illumination,
    wavelengths: ArrayLike,
    *,
    values: "tuple[int, ...] | None" = None,
    tolerance: float = DEFAULT_TOLERANCE,
    plateau: int = DEFAULT_PLATEAU,
    progress: "Progress | None" = None,
    **options,
) -> ConvergenceReport:
    """Sweep a solver's accuracy knob until the answer stops moving.

    Stops at the first plateau, so a well-behaved case costs only the rungs it
    needed. A case that never plateaus runs the whole ladder and says so --
    :attr:`ConvergenceReport.passed` is False and
    ``scan.provenance.converged`` is ``False``, which is a stronger and more
    useful statement than the ``None`` it started as.

    Parameters
    ----------
    values
        Explicit ladder, ascending. Defaults to :func:`default_ladder`.
    tolerance
        Largest efficiency change treated as "stopped moving".
    plateau
        Consecutive differences that must fall below ``tolerance``. See the
        module docstring; below 1 is refused.
    progress
        Called ``(rungs_done, len(ladder))``.

        Two granularities, deliberately different. The callback is **invoked
        once per wavelength**, so a
        :class:`~gratinglab.solvers.base.SolveCancelled` raised from it lands
        within milliseconds rather than at the end of a rung -- and a rung is
        the expensive unit here, up to a quarter of the whole sweep. But the
        numbers it *reports* are per rung, because they are the only ones that
        are both monotone and meaningful: a per-wavelength counter would
        restart on every rung and drive a bar backwards, and composing one
        across rungs would misrepresent progress anyway, since rungs double in
        cost and the last is half the total work.

        So a bar sits still within a rung and steps once per rung. That is the
        honest shape of this computation.

        Unlike the :class:`~gratinglab.solvers.base.Solver` contract, this
        does **not** end at ``(total, total)``: the sweep stops at the first
        plateau, so a well-behaved case genuinely never runs the tail of the
        ladder. Finishing early is the good outcome, and running the counter
        up to the end to make a bar look tidy would claim work that was never
        done. A caller learns the sweep is over by it returning.
    **options
        Passed to every solve, so a sweep runs at the same settings as the
        production run it is certifying. The knob itself may not appear here.

    Raises
    ------
    SolveCancelled
        Propagated untouched from ``progress``. **No partial report is
        returned**: a sweep that was stopped has not demonstrated anything,
        and handing back a report that implied otherwise would be the mistake
        this module exists to prevent.
    """
    knob = _knob_of(solver)
    if plateau < 1:
        raise ValueError(
            f"plateau must be at least 1, got {plateau}; a sweep that requires "
            "no agreement at all demonstrates nothing"
        )
    if knob in options:
        raise ValueError(
            f"{knob!r} is the knob being swept and cannot also be pinned in "
            "options -- pass `values` to control the ladder"
        )

    ladder = tuple(values) if values is not None else default_ladder(solver)
    if len(ladder) < 2:
        raise ValueError("a sweep needs at least two knob values to compare")
    if list(ladder) != sorted(set(ladder)):
        raise ValueError(f"values must be ascending and distinct, got {ladder}")

    solved: list[int] = []
    scans: list[EfficiencyScan] = []
    differences: list[float] = []
    skipped: list[tuple[int, str]] = []
    converged_at: int | None = None

    for rung, value in enumerate(ladder):
        forwarded = _rung_progress(solver, progress, rung, len(ladder))
        try:
            scan = solver.solve(problem, illumination, wavelengths, **{knob: value},
                                **forwarded, **options)
        except (ValueError, UnsupportedConfiguration) as exc:
            # Narrow on purpose. A coarse rung below a solver's own floor
            # raises ValueError (scalar's Nyquist guard) or refuses outright;
            # a TypeError or AttributeError is a bug in the call and must not
            # be quietly logged as a skipped rung.
            if solved:
                # A refusal *after* a success is not a coarse-grid floor; the
                # ladder walked into something real and hiding it would leave
                # a sweep that looks complete and is not.
                raise
            skipped.append((value, str(exc)))
            continue

        if scans:
            _require_same_grid(scans[-1], scan, solved[-1], value, knob)
            differences.append(_largest_change(scans[-1], scan))
        solved.append(value)
        scans.append(scan)

        if len(differences) >= plateau and all(
            d <= tolerance for d in differences[-plateau:]
        ):
            # The coarsest value whose successors all agree -- what a
            # production run should use, and cheaper than what we just solved.
            converged_at = solved[-(plateau + 1)]
            break

    if not scans:
        raise RuntimeError(
            f"no rung of the ladder {ladder} could be solved; refusals: "
            + "; ".join(f"{v}: {why}" for v, why in skipped)
        )
    if len(scans) < 2:
        raise RuntimeError(
            f"only one rung of the ladder {ladder} could be solved "
            f"({solved[0]}), so nothing could be compared"
        )

    report = ConvergenceReport(
        knob=knob,
        values=tuple(solved),
        differences=tuple(differences),
        tolerance=tolerance,
        plateau=plateau,
        converged_at=converged_at,
        scan=scans[-1],
        skipped=tuple(skipped),
    )
    return dataclasses.replace(report, scan=_stamp(scans[-1], report))


def _rung_progress(solver, progress: "Progress | None", rung: int, rungs: int) -> dict:
    """``{"progress": ...}`` reporting rungs while being called per wavelength.

    The wavelength counter is discarded on purpose -- see
    :func:`check_convergence`. What it buys is the call *frequency*: every
    wavelength is a chance for the caller to raise, so cancelling a sweep does
    not have to wait out the rung it is in.
    """
    if progress is None or not solver.capabilities.reports_progress:
        return {}
    return {"progress": lambda _done, _total: progress(rung, rungs)}


def _knob_of(solver) -> str:
    knob = solver.capabilities.accuracy_knob
    if knob is None:
        raise ValueError(
            f"{solver.capabilities.name} declares no accuracy_knob, so there is "
            "nothing to sweep. A closed-form method may well be exact, but that "
            "claim belongs to whoever can defend it, not to a sweep of a "
            "parameter that does not exist."
        )
    return knob


def _require_same_grid(
    before: EfficiencyScan, after: EfficiencyScan, coarse: int, fine: int, knob: str
) -> None:
    """Refuse to compare scans that are not about the same orders.

    The knob controls accuracy, not geometry, so this should never fire -- the
    order grid comes from the grating equation. It exists because a silent
    broadcast between differently-shaped arrays would produce a number that
    looks like a convergence measure and is not.
    """
    if not np.array_equal(before.orders, after.orders):
        raise ValueError(
            f"{knob}={coarse} and {knob}={fine} returned different order grids "
            f"({before.orders} vs {after.orders}); the accuracy knob must not "
            "change which orders exist"
        )
    if not np.array_equal(before.wavelengths, after.wavelengths):
        raise ValueError(
            f"{knob}={coarse} and {knob}={fine} returned different wavelength "
            "grids"
        )
    if not np.array_equal(before.propagating, after.propagating):
        raise ValueError(
            f"{knob}={coarse} and {knob}={fine} disagree about which orders "
            "propagate, which is geometry and cannot depend on accuracy"
        )


def _largest_change(before: EfficiencyScan, after: EfficiencyScan) -> float:
    """Largest absolute efficiency change, over every wavelength and order."""
    return float(np.abs(np.asarray(after.efficiency) - np.asarray(before.efficiency)).max())


def _stamp(scan: EfficiencyScan, report: ConvergenceReport) -> EfficiencyScan:
    """A copy of ``scan`` whose provenance carries the verdict and its evidence.

    ``converged`` becomes a real boolean either way. False is a result, not a
    failure to produce one: it says the sweep ran and the answer was still
    moving, which is exactly what a reader needs in order not to quote the
    number.
    """
    notes = dict(scan.provenance.notes)
    notes["convergence"] = report.as_notes()
    provenance = dataclasses.replace(
        scan.provenance,
        converged=report.passed,
        truncation=report.values[-1],
        notes=notes,
    )
    if not report.passed:
        provenance = provenance.with_warning(
            f"convergence not demonstrated: {report.knob} swept to "
            f"{report.values[-1]} and the largest change was still "
            f"{report.final_difference:.3g}, above the {report.tolerance:g} "
            "tolerance. Treat this efficiency as indicative, not defensible."
        )
    return dataclasses.replace(scan, provenance=provenance)


def converged_scan(
    solver,
    problem: Problem,
    illumination: Illumination,
    wavelengths: ArrayLike,
    **kwargs,
) -> EfficiencyScan:
    """:func:`check_convergence`, keeping only the stamped scan.

    For callers that want a defensible result rather than the study behind it.
    The study is still on ``scan.provenance.notes["convergence"]``.
    """
    return check_convergence(solver, problem, illumination, wavelengths, **kwargs).scan
