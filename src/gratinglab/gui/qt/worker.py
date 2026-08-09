"""Solving off the UI thread.

The scalar solver takes about 70 ms, so today this buys nothing a user would
notice. RCWA at N ≈ 200–400 and the integral method will take seconds to
minutes, and retrofitting an event flow around a window that has already
frozen is worse than building for it once.

**Cancel now stops the work**, not merely the waiting. It used to do only the
latter -- ``sweep`` reaches a NumPy-bound call with no check point and a
Python thread cannot be killed, so there was nothing to interrupt from out
here. `Capabilities.reports_progress` changed that: the solver calls back once
per wavelength, and a callback that raises
:class:`~gratinglab.solvers.base.SolveCancelled` unwinds it.

**The token is still needed.** Real cancellation does not remove the race it
was covering for: the worker may already have finished when Cancel is clicked,
so a result can still be in flight and must still be recognised as stale on
arrival. Belt and braces, for two different situations.

**A fresh event per solve, created by the window and carried in the request.**
Not one owned by this object and cleared at the top of `run` -- that has a real
race, in which a Cancel clicked between the request and the clear is silently
lost.
"""

from __future__ import annotations

import threading
from typing import TYPE_CHECKING, NamedTuple

from PySide6.QtCore import QObject, Signal, Slot

from ...solvers.base import SolveCancelled

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ...checks import EnergyReport
    from ...result import EfficiencyScan
    from ..state import Parsed

__all__ = ["SolveResult", "SolveWorker"]


class SolveResult(NamedTuple):
    """A finished solve and its energy verdict, computed together.

    Both cross the thread boundary by reference and neither is copied.
    ``EfficiencyScan`` is a frozen dataclass whose arrays are marked
    ``write=False``, so there is nothing here another thread could mutate --
    the immutability work in `result.py` is what makes this safe without a
    lock.
    """

    scan: "EfficiencyScan"
    energy: "EnergyReport"


class SolveWorker(QObject):
    """Runs one solve at a time on a thread of its own, for whichever solver
    asked.

    Lives on a long-lived ``QThread`` owned by the window, shared by every
    solver tab -- there is only ever one solve in flight window-wide (see
    `main_window.py`). ``token`` is echoed back untouched so the window can
    tell a result it still wants from one it cancelled; ``method`` is echoed
    the same way so the window knows which tab's result this is, since
    several tabs share this one worker.
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
        from ...checks import check_energy_balance
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
            self.finished.emit(
                token, method, SolveResult(scan, check_energy_balance(scan))
            )
        except SolveCancelled:
            # Before the broad handler on purpose. The user asked for this, so
            # it is not a failure and rendering it as one would be the same
            # category mistake M8 fixed for the no-coating default.
            self.cancelled.emit(token, method)
        except Exception as exc:  # noqa: BLE001 - see below
            # Deliberately broad. A Python exception escaping a slot under
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

        Up to ten solves, which is what makes this the first operation in the
        project long enough for a progress bar to be worth having -- and the
        first a user is likely to want out of.

        The knob is dropped from ``options``: the sweep chooses it. Passing
        both would be refused by `check_convergence` anyway, and refusing here
        with a clear reason beats surfacing that as a solver error.
        """
        from ...checks import check_energy_balance
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
                SolveResult(report.scan, check_energy_balance(report.scan)),
                report,
            )
        except SolveCancelled:
            self.cancelled.emit(token, method)
        except Exception as exc:  # noqa: BLE001 - see `run`
            self.failed.emit(token, method, f"{type(exc).__name__}: {exc}")

    def _reporter(self, token: int, cancel: threading.Event):
        """The callback handed to the solver: emit, then bail if asked to.

        Throttled to whole percent. A 200-wavelength scan otherwise emits 200
        queued cross-thread signals for a bar that has ~100 distinguishable
        states, and the ones it drops are the ones nobody could have seen.
        The *cancellation* check is not throttled -- it runs on every call, so
        a Cancel click lands at the next wavelength rather than the next
        percent.
        """
        last = -1

        def report(done: int, total: int) -> None:
            nonlocal last
            percent = done * 100 // total if total else 100
            if percent != last:
                last = percent
                self.progress.emit(token, done, total)
            if cancel.is_set():
                raise SolveCancelled(f"cancelled at {done} of {total}")

        return report
