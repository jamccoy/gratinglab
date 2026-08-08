"""Solving off the UI thread.

The scalar solver takes about 70 ms, so today this buys nothing a user would
notice. RCWA at N ≈ 200–400 and the integral method will take seconds to
minutes, and retrofitting an event flow around a window that has already
frozen is worse than building for it once.

**Cancel abandons a result; it does not interrupt one.** ``sweep`` reaches a
single NumPy-bound call with no check point, and a Python thread cannot be
killed, so there is nothing to interrupt from out here. Cancelling therefore
bumps a token: the in-flight result still arrives, is recognised as stale, and
is dropped. The window stops waiting; the CPU does not stop working. The panel
says exactly that rather than implying more.

The honest fix is upstream, and is worth doing when a slow solver lands: an
optional ``progress`` callback on the `Solver` protocol, invoked per
wavelength, gives real progress *and* real cancellation (raise a sentinel from
inside the callback) in one change. See `docs/roadmap.md`.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

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

    @Slot(int, str, object, dict)
    def run(
        self, token: int, method: str, geometry: "Parsed", options: dict
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
            )[0]
            self.finished.emit(
                token, method, SolveResult(scan, check_energy_balance(scan))
            )
        except Exception as exc:  # noqa: BLE001 - see below
            # Deliberately broad. A Python exception escaping a slot under
            # PySide6 does not surface as a friendly dialog: depending on
            # version and excepthook state it prints a traceback and aborts
            # the process. Turning any failure into a signal keeps the window
            # alive and able to say what went wrong.
            self.failed.emit(token, method, f"{type(exc).__name__}: {exc}")
