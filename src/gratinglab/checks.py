r"""Physics self-checks: validate the *model*, not just the arithmetic.

The closed-form tests in ``tests/test_scalar.py`` confirm that a Fourier
integral is evaluated correctly. They cannot confirm that the integrand is the
right one, because the expected values are derived from the same formula the
solver implements. Two checks here close that gap, and neither needs reference
data or a closed form:

**Reciprocity** (:func:`check_reciprocity`) follows from Lorentz reciprocity and
constrains the *structure* of the phase function. It is sharp: on the scalar
solver it holds to ~1e-17, and breaking the :math:`\alpha \leftrightarrow \beta`
symmetry of :math:`\Phi` pushes the violation to ~4e-1 -- sixteen orders of
magnitude of discrimination.

**Energy balance** (:func:`check_energy_balance`) catches the other class of
error. Summed efficiency exceeding unity is unphysical for any passive
structure, whatever method produced it. This is the check that found the
unphysical finite-conductivity run in the reference corpus.

Both take a solver rather than a scan where possible, so every backend added
later inherits them without new code.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .geometry import beta as beta_from_sin
from .geometry import is_propagating, sin_beta
from .illumination import Illumination
from .problem import Problem
from .result import EfficiencyScan

__all__ = [
    "ReciprocityReport",
    "EnergyReport",
    "check_reciprocity",
    "check_energy_balance",
]

#: Illumination is undefined at exactly grazing, so reciprocal geometries
#: within this margin of +/-90 degrees are skipped rather than clamped.
_GRAZING_MARGIN_DEG = 0.5


@dataclass(frozen=True, slots=True)
class ReciprocityReport:
    """Outcome of a reciprocity check."""

    max_violation: float
    worst_order: int | None
    worst_wavelength: float | None
    pairs_tested: int
    tolerance: float

    @property
    def passed(self) -> bool:
        return self.pairs_tested > 0 and self.max_violation <= self.tolerance

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        if not self.pairs_tested:
            return "reciprocity: no testable pairs"
        verdict = "pass" if self.passed else "FAIL"
        return (
            f"reciprocity: {verdict} -- worst {self.max_violation:.3e} "
            f"at order {self.worst_order}, lambda={self.worst_wavelength:g} nm "
            f"({self.pairs_tested} pairs, tol {self.tolerance:.1e})"
        )


@dataclass(frozen=True, slots=True)
class EnergyReport:
    """Outcome of an energy-balance check."""

    total: NDArray[np.float64]
    wavelengths: NDArray[np.float64]
    tolerance: float
    lossless: bool

    @property
    def max_excess(self) -> float:
        """How far above unity the summed efficiency ever goes.

        Positive means unphysical: a passive grating cannot return more power
        than it receives, by any method.
        """
        return float(self.total.max() - 1.0)

    @property
    def max_deficit(self) -> float:
        """How far below unity. Expected for absorbing or non-rigorous methods."""
        return float(1.0 - self.total.min())

    @property
    def unphysical(self) -> NDArray[np.bool_]:
        """Wavelengths where efficiency exceeds unity beyond tolerance."""
        return self.total > 1.0 + self.tolerance

    @property
    def passed(self) -> bool:
        if self.lossless:
            return bool(np.abs(self.total - 1.0).max() <= self.tolerance)
        return not bool(self.unphysical.any())

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        verdict = "pass" if self.passed else "FAIL"
        return (
            f"energy balance: {verdict} -- total in "
            f"[{self.total.min():.4f}, {self.total.max():.4f}], "
            f"{int(self.unphysical.sum())} of {len(self.total)} points above unity"
        )


def check_reciprocity(
    solver,
    problem: Problem,
    illumination: Illumination,
    wavelengths: ArrayLike,
    *,
    tolerance: float = 1e-9,
    max_orders: int | None = 12,
    **options,
) -> ReciprocityReport:
    r"""Verify Lorentz reciprocity for a solver.

    Reversing the diffracted ray of order :math:`m` gives an incidence azimuth
    :math:`\alpha' = \beta_m`, and order :math:`m` then exits along
    :math:`\alpha`. Reciprocity requires

    .. math::
        \mathscr{E}_m(\alpha) = \mathscr{E}_m(\beta_m)

    This constrains the phase function to be symmetric under
    :math:`\alpha \leftrightarrow \beta_m`, which is exactly the structure
    :math:`\cos\alpha + \cos\beta_m` in the scalar model. A solver that used
    :math:`2\cos\alpha`, or subtracted rather than added, fails immediately.

    Parameters
    ----------
    max_orders
        Cap on how many orders to test per wavelength, since each one costs an
        extra solve. ``None`` tests all of them.

    Notes
    -----
    Orders diffracting within half a degree of grazing are skipped: the
    reciprocal illumination is undefined there, and clamping it would compare
    two different geometries.
    """
    wavelengths = np.atleast_1d(np.asarray(wavelengths, dtype=np.float64))

    # Collected rather than tracked with a running maximum. A running max
    # seeded at 0.0 never updates when a solver is *exactly* reciprocal, which
    # left the report with worst_order=None on platforms whose BLAS returns a
    # bitwise-zero difference -- the better the solver, the worse the report.
    violations: list[tuple[float, int, float]] = []

    for wavelength in wavelengths:
        forward = solver.solve(problem, illumination, [wavelength], **options).at(
            wavelength
        )
        orders = forward.propagating_orders
        if max_orders is not None and len(orders) > max_orders:
            # Keep the strongest orders: a violation there matters most, and
            # near-zero orders make the comparison meaningless.
            strength = np.array([forward[int(m)] for m in orders])
            orders = orders[np.argsort(strength)[::-1][:max_orders]]

        for order in orders:
            sine = sin_beta(
                order,
                wavelength,
                problem.period,
                illumination.sin_alpha,
                illumination.sin_gamma,
            )
            if not is_propagating(sine):
                continue
            beta_deg = float(np.degrees(beta_from_sin(sine)))
            if abs(beta_deg) > 90.0 - _GRAZING_MARGIN_DEG:
                continue

            reversed_illumination = Illumination(
                alpha_deg=beta_deg,
                gamma_deg=illumination.gamma_deg,
                polarization=illumination.polarization,
            )
            backward = solver.solve(
                problem, reversed_illumination, [wavelength], **options
            ).at(wavelength)

            violations.append(
                (
                    abs(forward[int(order)] - backward[int(order)]),
                    int(order),
                    float(wavelength),
                )
            )

    if not violations:
        return ReciprocityReport(
            max_violation=0.0,
            worst_order=None,
            worst_wavelength=None,
            pairs_tested=0,
            tolerance=tolerance,
        )

    worst, worst_order, worst_wavelength = max(violations, key=lambda item: item[0])
    return ReciprocityReport(
        max_violation=worst,
        worst_order=worst_order,
        worst_wavelength=worst_wavelength,
        pairs_tested=len(violations),
        tolerance=tolerance,
    )


def check_energy_balance(
    scan: EfficiencyScan,
    *,
    tolerance: float = 1e-6,
    lossless: bool = False,
) -> EnergyReport:
    """Verify that summed efficiency is physically admissible.

    Parameters
    ----------
    lossless
        Require ``sum == 1`` exactly (within tolerance). Use for a perfectly
        conducting or otherwise lossless structure, where the sum rule is a
        theorem. Default ``False`` only requires ``sum <= 1``, which holds for
        **any** passive grating by **any** method -- a weaker claim, but one
        that never has a legitimate exception.

    Notes
    -----
    The default direction matters. A deficit is ordinary: power goes into
    absorption, into evanescent orders, or is missed by an approximate method.
    An *excess* is not, and this is the check that identified the unphysical
    finite-conductivity run in the reference corpus, where summed efficiency
    reached 3.6.
    """
    return EnergyReport(
        total=np.asarray(scan.total, dtype=np.float64),
        wavelengths=scan.wavelengths,
        tolerance=tolerance,
        lossless=lossless,
    )
