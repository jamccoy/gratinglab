r"""Dense solves that report how well conditioned they were.

The rigorous solvers each spend their time in a few dense linear systems: the
integral solver's Nystrom matrix, RCWA's mode matrices and S-matrix
denominators. A system that is close to singular returns a confident answer
with fewer correct digits than the efficiency it feeds, and nothing
downstream can see that. Energy balance cannot either: RCWA's holds to
rounding by construction, so it says nothing about lost digits.

So every such solve goes through here. :func:`solve` and :func:`inverse`
factor the matrix once with LAPACK's ``getrf`` (what ``numpy.linalg.solve``
calls), solve with ``getrs`` from that factorisation, and ask ``gecon`` for
the reciprocal 1-norm condition number from the same factors. ``gecon`` is
Hager and Higham's estimator: ``O(n^2)`` after the ``O(n^3)`` factorisation
the solve needed anyway, so the estimate is close to free, and it is reliable
to within a small factor of the true ``kappa_1``, which is all a warning
threshold needs.

``kappa * eps`` bounds the relative error a backward-stable solve can leave,
and it is a worst case: measured at a Rayleigh anomaly, both rigorous solvers
stay accurate with ``kappa`` far past ``1/sqrt(eps)``, because the
near-singular direction is the grazing order's own amplitude and the
efficiencies do not lean on it (``docs/theory/rcwa.md`` section 9,
``docs/theory/integral.md``). So each backend sets its own threshold, placed
by measurement well above what it reaches on ordinary problems: crossing it
means the problem is unusual, not that the answer is wrong.

:class:`ConditionLog` keeps the worst estimate seen per named group, so a
backend can record one number per group in ``Provenance.notes``. An exactly
singular matrix raises :class:`numpy.linalg.LinAlgError`, as
``numpy.linalg.solve`` would.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import get_lapack_funcs

__all__ = ["ConditionLog", "solve", "inverse", "ESTIMATOR"]

#: What the recorded numbers are, for ``Provenance.notes``.
ESTIMATOR = "LAPACK gecon: 1-norm condition estimate from the solve's own LU"


def _factor(a: NDArray) -> tuple:
    """``(lu, piv, getrs, kappa)`` for a square ``a``."""
    a = np.asarray(a)
    if a.ndim != 2 or a.shape[0] != a.shape[1]:
        raise ValueError(f"expected a square matrix, got shape {a.shape}")
    getrf, getrs, gecon, lange = get_lapack_funcs(
        ("getrf", "getrs", "gecon", "lange"), (a,)
    )
    norm = float(lange("1", a))
    lu, piv, info = getrf(a)
    if info > 0:
        raise np.linalg.LinAlgError("Singular matrix")
    if info < 0:  # pragma: no cover - an argument error is a bug here
        raise ValueError(f"getrf: illegal argument {-info}")
    if not np.isfinite(norm):
        kappa = np.inf
    else:
        rcond, info = gecon(lu, norm, norm="1")
        kappa = np.inf if info != 0 or rcond <= 0.0 else 1.0 / float(rcond)
    return lu, piv, getrs, kappa


def solve(a: NDArray, b: NDArray) -> tuple[NDArray, float]:
    """``(x, kappa_1)``: the solution of ``a x = b`` and ``a``'s condition."""
    lu, piv, getrs, kappa = _factor(a)
    b = np.asarray(b)
    x, info = getrs(lu, piv, b.astype(lu.dtype, copy=False))
    if info != 0:  # pragma: no cover - an argument error is a bug here
        raise ValueError(f"getrs: illegal argument {-info}")
    return x, kappa


def inverse(a: NDArray) -> tuple[NDArray, float]:
    """``(a^-1, kappa_1)``, the inverse from the same factorisation."""
    a = np.asarray(a)
    return solve(a, np.eye(len(a), dtype=a.dtype))


class ConditionLog:
    """The worst condition estimate seen, per named group of systems."""

    __slots__ = ("worst",)

    def __init__(self) -> None:
        self.worst: dict[str, float] = {}

    def record(self, group: str, kappa: float) -> None:
        if kappa > self.worst.get(group, 0.0):
            self.worst[group] = float(kappa)

    def merge(self, other: "ConditionLog | dict[str, float]") -> None:
        items = other.worst if isinstance(other, ConditionLog) else other
        for group, kappa in items.items():
            self.record(group, kappa)

    def warning(
        self, method: str, flagged: list[float], threshold: float
    ) -> str | None:
        """The provenance warning, or ``None`` when nothing crossed.

        ``flagged`` are the wavelengths at which some system crossed
        ``threshold``; the groups over it are named from ``worst``.
        """
        over = {g: k for g, k in self.worst.items() if k > threshold}
        if not over or not flagged:
            return None
        groups = ", ".join(f"{g} {k:.2g}" for g, k in sorted(over.items()))
        span = f"{min(flagged):.6g}-{max(flagged):.6g} nm"
        return (
            f"{method}: a dense system's condition estimate exceeds "
            f"{threshold:.0e} at {len(flagged)} wavelength(s) in {span} "
            f"(worst: {groups}), far above anything this solver reaches on "
            f"an ordinary problem; the solve's error bound kappa*eps is "
            f"{threshold * 2.2e-16:.0e} or worse there, so check those rows "
            "against their neighbours before relying on them"
        )

    def solve(self, group: str, a: NDArray, b: NDArray) -> NDArray:
        x, kappa = solve(a, b)
        self.record(group, kappa)
        return x

    def inverse(self, group: str, a: NDArray) -> NDArray:
        x, kappa = inverse(a)
        self.record(group, kappa)
        return x
