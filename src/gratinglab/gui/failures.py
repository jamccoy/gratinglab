"""How a solve that raised is described on a solver tab.

Called by the worker (:mod:`gratinglab.gui.qt.worker`) for every exception
other than a cancel, and its text is what the tab shows after "solve failed:".

Two kinds are told apart. A **refusal** is a solver saying this problem is
outside what it handles: :class:`~gratinglab.solvers.base.UnsupportedConfiguration`,
a profile it cannot represent, or a ``ValueError`` from checking the inputs.
Its message already says what to change, so it is shown as written, without
the Python class name. Anything else is **unexpected**, so the class name is
kept and the text asks for a report.
"""

from __future__ import annotations

from ..profiles._types import ProfileRepresentationError
from ..solvers.base import UnsupportedConfiguration
from .docs import REPOSITORY

__all__ = ["REFUSALS", "is_refusal", "describe_failure"]

#: Exceptions a solver raises on purpose when a problem is outside its scope.
REFUSALS = (UnsupportedConfiguration, ProfileRepresentationError, ValueError)


def is_refusal(exc: BaseException) -> bool:
    """True when the solver declined the problem rather than broke on it."""
    return isinstance(exc, REFUSALS)


def describe_failure(exc: BaseException) -> str:
    """One line or paragraph for the tab, saying what went wrong."""
    message = str(exc).strip() or "(no message)"
    if is_refusal(exc):
        return f"this solver cannot take this problem. {message}"
    return (
        f"unexpected {type(exc).__name__}: {message}. This is probably a bug; "
        f"please report it at {REPOSITORY}/issues"
    )
