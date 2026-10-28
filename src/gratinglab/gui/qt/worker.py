"""Solve off the UI thread.

The scalar solver takes about 70 ms, so at present this gives no benefit a
user would notice. RCWA at N ≈ 200–400 and the integral method will take
seconds to minutes, and retrofitting an event flow onto a window that already
freezes is harder than designing for it from the start.

Cancel now stops the work itself, and not only the waiting. Previously it did
only the latter, because ``sweep`` reaches a NumPy-bound call with no check
point and a Python thread cannot be killed, so there was nothing to interrupt
from here. `Capabilities.reports_progress` changed that: the solver calls back
once per wavelength, and a callback that raises
:class:`~endstation.progress.Cancelled` (or its subclass
:class:`~gratinglab.solvers.base.SolveCancelled`) unwinds it.

The token is still needed, because real cancellation does not remove the race
it was covering. The worker may already have finished when Cancel is clicked,
so a result can still be in flight and must still be recognised as stale when
it arrives. The two mechanisms therefore cover two different situations.

A fresh event is used for each solve, created by the window and carried in
the request. An event owned by this object and cleared at the top of `run`
would introduce a real race, in which a Cancel clicked between the request
and the clear is silently lost.
"""

from __future__ import annotations

import threading
from typing import TYPE_CHECKING, NamedTuple

from endstation.progress import Cancelled, cancelling, throttled
from PySide6.QtCore import QObject, Signal, Slot

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ...checks import EnergyReport
    from ...result import EfficiencyScan
    from ..state import Parsed

__all__ = ["SolveResult", "SolveWorker"]


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


class SolveWorker(QObject):
    """Run one solve at a time on a dedicated thread, for whichever solver asked.

    The worker lives on a long-lived ``QThread`` owned by the window and is
    shared by every solver tab, since only one solve is ever in flight in the
    whole window (see `main_window.py`). ``token`` is echoed back unchanged so
    that the window can distinguish a result it still wants from one it
    cancelled, and ``method`` is echoed in the same way so that the window
    knows which tab the result belongs to, since several tabs share this one
    worker.
    """

    finished = Signal(int, str, object)  # token, method, SolveResult
    failed = Signal(int, str, str)  # token, method, message
    #: token, method -- a stop the user asked for, which is not an error and
    #: must not be rendered as one.
    cancelled = Signal(int, str)
    progress = Signal(int, int, int)  # token, done, total
    #: token, method, SolveResult, ConvergenceReport. Carries the result
    #: separately from the study so a tab can reuse its ordinary result path
    #: and show the evidence through the provenance panel it already has.
    converged = Signal(int, str, object, object)

    @Slot(int, str, object, dict, object)
    def run(
        self,
        token: int,
        method: str,
        geometry: "Parsed",
        options: dict,
        cancel: threading.Event,
    ) -> None:
        from ...compare import sweep

        try:
            scan = sweep(
                geometry.problem,
                geometry.illumination,
                geometry.wavelengths,
                [method],
                options={method: options},
                progress=self._reporter(token, cancel),
            )[0]
            self.finished.emit(token, method, SolveResult(scan, _energy(scan)))
        except Cancelled:
            # Before the broad handler: the user asked for this, so
            # it is not a failure and rendering it as one would be the same
            # category mistake M8 fixed for the no-coating default.
            self.cancelled.emit(token, method)
        except Exception as exc:  # noqa: BLE001 - see below
            # Broad, because a Python exception escaping a slot under
            # PySide6 does not surface as a friendly dialog: depending on
            # version and excepthook state it prints a traceback and aborts
            # the process. Turning any failure into a signal keeps the window
            # alive and able to say what went wrong.
            self.failed.emit(token, method, f"{type(exc).__name__}: {exc}")

    @Slot(int, str, object, dict, object)
    def run_convergence(
        self,
        token: int,
        method: str,
        geometry: "Parsed",
        options: dict,
        cancel: threading.Event,
    ) -> None:
        """Sweep the solver's accuracy knob instead of solving once.

        A sweep runs up to ten solves, which makes it the first operation in
        the project long enough to merit a progress bar, and the first that a
        user is likely to want to cancel.

        The knob is dropped from ``options``, because the sweep chooses it.
        `check_convergence` would refuse both being passed in any case, and
        refusing here with a clear reason is preferable to surfacing the
        refusal as a solver error.
        """
        from ...convergence import check_convergence
        from ...solvers.base import get_solver

        try:
            solver = get_solver(method)
            knob = solver.capabilities.accuracy_knob
            report = check_convergence(
                solver,
                geometry.problem,
                geometry.illumination,
                geometry.wavelengths,
                progress=self._reporter(token, cancel),
                **{k: v for k, v in options.items() if k != knob},
            )
            self.converged.emit(
                token,
                method,
                SolveResult(report.scan, _energy(report.scan)),
                report,
            )
        except Cancelled:
            self.cancelled.emit(token, method)
        except Exception as exc:  # noqa: BLE001 - see `run`
            self.failed.emit(token, method, f"{type(exc).__name__}: {exc}")

    def _reporter(self, token: int, cancel: threading.Event):
        """The callback handed to the solver, which emits progress and stops if asked.

        Progress is throttled to whole percent. Otherwise a 200-wavelength
        scan would emit 200 queued cross-thread signals for a bar that has
        ~100 distinguishable states, and the signals that are dropped are ones
        that could not have been seen. The cancellation check is not
        throttled: it runs on every call, so a Cancel click takes effect at
        the next wavelength instead of the next percent.

        Both halves are endstation's: :func:`~endstation.progress.throttled`
        for the bar and :func:`~endstation.progress.cancelling` for the stop,
        which raises :class:`~endstation.progress.Cancelled`, the base class
        of :class:`~gratinglab.solvers.base.SolveCancelled`.
        """

        def emit(done: int, total: int) -> None:
            self.progress.emit(token, done, total)

        return cancelling(throttled(emit), cancel.is_set)
