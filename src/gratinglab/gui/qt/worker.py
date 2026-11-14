"""The work the efficiency window runs off the UI thread.

The thread is endstation's (`endstation.qt.jobs.JobRunner`, owned by the
window); this module is only the work: a solve, a convergence sweep and a
comparison, each a function the runner calls with a ``progress`` callback.

The scalar solver takes about 70 ms. RCWA at N of a few hundred and the
integral method take seconds to minutes, which is why solving happens here
and not on the window's thread.

Cancel stops the work itself. The solvers call ``progress`` once per
wavelength (`Capabilities.reports_progress`), and the runner's callback
raises :class:`~endstation.progress.Cancelled`, the base class of
:class:`~gratinglab.solvers.base.SolveCancelled`, once Cancel is clicked, so
the solve unwinds at the next wavelength. A result that was already on its
way when Cancel landed is dropped by the runner, which is what the window's
token used to do.
"""

from __future__ import annotations

import inspect
import os
from typing import TYPE_CHECKING, NamedTuple

from endstation.progress import Cancelled

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ...checks import EnergyReport
    from ...result import EfficiencyScan
    from ..state import Parsed

__all__ = ["SolveResult", "compare", "converge", "solve"]


class SolveResult(NamedTuple):
    """A finished solve and its energy verdict, computed together.

    Both cross the thread boundary by reference and neither is copied.
    ``EfficiencyScan`` is a frozen dataclass whose arrays are marked
    ``write=False``, so there is nothing here that another thread could
    mutate; the immutability work in `result.py` makes this safe without a
    lock.
    """

    scan: "EfficiencyScan"
    energy: "EnergyReport"


def _energy(scan) -> "EnergyReport":
    """Check energy balance at a tolerance the conserved quantity can meet.

    The 1e-6 default is appropriate for a one-sided quantity. A scalar sum, or
    a perfect conductor's, approaches 1 from below, so any excess at all is a
    modelling error and should be caught at machine precision. This check
    found a reference run summing to 3.6.

    A scan carrying an absorbed fraction is different in kind. ``R + A = 1`` is
    a theorem constrained from both sides, so the discretisation error falls
    above unity as often as below. A converged run measured +2.7e-6 and was
    rendered as "EXCEEDS UNITY" in red, which is the "nothing correct looks
    wrong" failure that `gui/provenance.py` is written to prevent. Such a scan
    is therefore checked at the same deviation that the solver itself reports,
    so that the panel's red text and the solver's own warning agree instead of
    contradicting each other.
    """
    from ...checks import check_energy_balance
    from ...solvers.integral import _ENERGY_DEVIATION_WARN

    if scan.absorption is None:
        return check_energy_balance(scan)
    return check_energy_balance(scan, tolerance=_ENERGY_DEVIATION_WARN)


#: Threads a scan's wavelengths are solved on, for the solvers that take
#: ``workers``. The library defaults to 1, and nothing in the window passed
#: anything else, so a scan used one core of however many there were. A
#: 32-wavelength integral scan on MC1 (16 cores) took 5.05 s at 1, 1.69 s at
#: 4, 1.34 s at 8 and 1.22 s at 12, identical to the serial result each time,
#: so 8 is where it stops paying. One core is left for the window itself.
_SCAN_WORKERS = max(1, min(8, (os.cpu_count() or 2) - 1))


def _with_workers(method: str, options: dict) -> dict:
    """``options`` with :data:`_SCAN_WORKERS` added, if ``method`` takes it.

    A ``workers`` already present is kept. The solver's own signature is what
    decides, so a solver that gains the option later is covered without an
    edit here, and one without it (scalar, the C-method) is never sent a
    keyword it would refuse.
    """
    from ...solvers.base import get_solver

    if "workers" in options:
        return options
    if "workers" not in inspect.signature(get_solver(method).solve).parameters:
        return options
    return {**options, "workers": _SCAN_WORKERS}


def solve(method: str, geometry: "Parsed", options: dict, *, progress) -> SolveResult:
    """One solver over the scan, with its energy verdict."""
    from ...compare import sweep

    scan = sweep(
        geometry.problem,
        geometry.illumination,
        geometry.wavelengths,
        [method],
        options={method: _with_workers(method, options)},
        progress=progress,
    )[0]
    return SolveResult(scan, _energy(scan))


def converge(method: str, geometry: "Parsed", options: dict, *, progress):
    """Sweep the solver's accuracy knob instead of solving once.

    Returns ``(SolveResult, ConvergenceReport)``: the result separately from
    the study, so a tab can reuse its ordinary result path and show the
    evidence through the provenance panel it already has.

    A sweep runs up to ten solves, which makes it long enough to merit a
    progress bar and the first operation a user is likely to cancel. The knob
    is dropped from ``options``, because the sweep chooses it;
    `check_convergence` would refuse both, and refusing here with a clear
    reason reads better than a solver error.
    """
    from ...convergence import check_convergence
    from ...solvers.base import get_solver

    solver = get_solver(method)
    knob = solver.capabilities.accuracy_knob
    report = check_convergence(
        solver,
        geometry.problem,
        geometry.illumination,
        geometry.wavelengths,
        progress=progress,
        **{k: v for k, v in _with_workers(method, options).items() if k != knob},
    )
    return SolveResult(report.scan, _energy(report.scan)), report


def compare(geometry: "Parsed", requested: tuple, options: dict, refused: dict,
            *, progress):
    """Solve one geometry with each solver in ``options``, then line them up.

    ``requested`` is every ticked solver in run order. ``options`` maps the
    ones still to run to their keyword arguments; ``refused`` holds any the
    window already turned away (a form error on that solver's own tab) with
    the reason, and those are never started.

    Each solver goes through :func:`~gratinglab.compare.sweep` on its own
    rather than all of them in one call, because one solver declining the
    problem must not cost the others their answer: an
    ``UnsupportedConfiguration`` from the C-method on a vertical facet is
    recorded against the C-method and the next solver starts. The scans that
    came back are put on one order set by :func:`~gratinglab.compare.align`.

    Progress is one bar over the whole run, in wavelengths, and a solver that
    refuses steps the bar past its block so it never stalls there. Cancel
    stops the run and nothing partial is returned: half a comparison would
    read as the solvers that happened to come first being the only ones that
    apply.
    """
    from ...compare import align, sweep
    from ...solvers.base import UnsupportedConfiguration
    from ..comparison import ComparisonOutcome
    from ..failures import describe_failure

    refused = dict(refused)
    per_method = len(geometry.wavelengths)
    total = max(1, len(options) * per_method)
    scans, energies, failed = [], {}, {}
    for index, (method, method_options) in enumerate(options.items()):
        offset = index * per_method

        def step(done: int, _total: int, offset=offset) -> None:
            progress(offset + done, total)

        try:
            scan = sweep(
                geometry.problem,
                geometry.illumination,
                geometry.wavelengths,
                [method],
                options={method: _with_workers(method, method_options)},
                progress=step,
            )[0]
        except Cancelled:
            raise
        except UnsupportedConfiguration as exc:
            refused[method] = str(exc)
        except Exception as exc:  # noqa: BLE001 - one solver's failure is recorded
            failed[method] = describe_failure(exc)
        else:
            scans.append(scan)
            energies[method] = _energy(scan)
            continue
        progress(offset + per_method, total)
    return ComparisonOutcome(
        requested=tuple(requested),
        comparison=align(scans) if scans else None,
        energies=energies,
        refused=refused,
        failed=failed,
    )
