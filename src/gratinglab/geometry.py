"""Core grating geometry: the generalized grating equation and order bookkeeping.

All angles in this module are in **radians**. Degrees are converted once, at the
public API boundary (see :mod:`gratinglab.illumination`).

Symbols follow ``docs/conventions.md`` -- in particular ``period`` and ``depth``
are always spelled out, because the two primary source references use ``d`` to
mean opposite things.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "sin_beta",
    "beta",
    "cos_beta",
    "order_range",
    "is_propagating",
    "facet_graze",
    "blaze_direction",
    "blaze_wavelength",
]


def sin_beta(
    order: ArrayLike,
    wavelength: float,
    period: float,
    sin_alpha: float,
    sin_gamma: float = 1.0,
) -> NDArray[np.float64]:
    r"""Return :math:`\sin\beta_m` from the generalized grating equation.

    .. math::
        \sin\alpha + \sin\beta_m = \frac{m\lambda}{p \sin\gamma}

    ``sin_gamma = 1`` is the in-plane (classical) case. Values with
    ``|sin_beta| > 1`` correspond to evanescent orders; they are returned as-is
    rather than clipped, so callers can test them with :func:`is_propagating`.

    Parameters
    ----------
    order
        Diffraction order index :math:`m`. Scalar or array.
    wavelength, period
        Both in nm (any consistent unit works -- only the ratio enters).
    sin_alpha
        Sine of the azimuthal incidence angle.
    sin_gamma
        Sine of the half-cone angle. Must be > 0.
    """
    if period <= 0:
        raise ValueError(f"period must be positive, got {period}")
    if wavelength <= 0:
        raise ValueError(f"wavelength must be positive, got {wavelength}")
    if not 0 < sin_gamma <= 1:
        raise ValueError(f"sin_gamma must lie in (0, 1], got {sin_gamma}")

    m = np.asarray(order, dtype=np.float64)
    return m * wavelength / (period * sin_gamma) - sin_alpha


def is_propagating(sin_beta_m: ArrayLike) -> NDArray[np.bool_]:
    """True where an order propagates, i.e. :math:`|\\sin\\beta_m| \\le 1`."""
    return np.abs(np.asarray(sin_beta_m, dtype=np.float64)) <= 1.0


def cos_beta(sin_beta_m: ArrayLike, transmitted: bool = False) -> NDArray[np.float64]:
    r"""Return :math:`\cos\beta_m`, signed by branch.

    Reflected orders have :math:`\cos\beta_m > 0`; transmitted orders have
    :math:`\cos\beta_m < 0` (``docs/conventions.md`` §4 -- we do not flip the
    sign of :math:`\sin\beta_m` for transmission).

    Evanescent orders yield ``nan``; they carry no propagation direction.
    """
    s = np.asarray(sin_beta_m, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        c = np.sqrt(1.0 - s**2)
    return -c if transmitted else c


def beta(sin_beta_m: ArrayLike, transmitted: bool = False) -> NDArray[np.float64]:
    r"""Return :math:`\beta_m` in radians, on the correct branch.

    Reflected orders land in :math:`[-\pi/2, \pi/2]`. Transmitted orders are
    reflected through :math:`\pi/2` so that :math:`\cos\beta_m < 0`.
    """
    s = np.asarray(sin_beta_m, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        b = np.arcsin(np.where(np.abs(s) <= 1.0, s, np.nan))
    return np.pi - b if transmitted else b


def order_range(
    wavelength: float,
    period: float,
    sin_alpha: float,
    sin_gamma: float = 1.0,
) -> NDArray[np.int64]:
    r"""Return every propagating order index, ascending.

    Solves :math:`|\sin\alpha + m\lambda/(p\sin\gamma)| \le 1` for integer
    :math:`m`. Always includes ``m = 0``.

    The result is contiguous, so it is safe to use as a column index basis --
    but do not: order mapping should come from an explicit index array, never
    from a positional offset (``docs/conventions.md`` §4).
    """
    scale = period * sin_gamma / wavelength
    lo = int(np.ceil((sin_alpha - 1.0) * scale))
    hi = int(np.floor((sin_alpha + 1.0) * scale))

    # Guard the closed-interval endpoints against floating-point drift: an order
    # exactly at grazing (|sin_beta| == 1) is a passing-off point and belongs in
    # the set, but rounding can push it either way.
    while lo <= hi and not is_propagating(
        sin_beta(lo, wavelength, period, sin_alpha, sin_gamma)
    ):
        lo += 1
    while hi >= lo and not is_propagating(
        sin_beta(hi, wavelength, period, sin_alpha, sin_gamma)
    ):
        hi -= 1

    return np.arange(lo, hi + 1, dtype=np.int64)


def facet_graze(gamma: float, blaze_angle: float, alpha: float) -> float:
    r"""Graze angle onto a sawtooth groove facet, in radians.

    .. math::
        \sin\zeta = \sin\gamma \, \cos(\delta - \alpha)

    This is the angle that must stay below the critical angle for total
    external reflection, :math:`\theta_c \approx \sqrt{2\delta_n}`.
    """
    return float(np.arcsin(np.sin(gamma) * np.cos(blaze_angle - alpha)))


def blaze_direction(blaze_angle: float, alpha: float) -> float:
    r"""Azimuthal direction of peak blaze response, :math:`\beta_b = 2\delta - \alpha`.

    Equivalent to specular reflection from a mirror tilted at the facet angle.
    """
    return 2.0 * blaze_angle - alpha


def blaze_wavelength(
    order: ArrayLike,
    period: float,
    blaze_angle: float,
    alpha: float,
    gamma: float,
) -> NDArray[np.float64]:
    r"""Blaze wavelength for the given order(s), in the units of ``period``.

    .. math::
        m \lambda_b = 2 p \sin\zeta \, \sin\delta,
        \qquad \sin\zeta = \sin\gamma \cos(\delta - \alpha)

    This is ISSI eq. following (15). It agrees exactly with the thesis form
    :math:`\lambda_b = 2 p \sin\gamma \sin\delta \cos(\delta-\alpha)/m`
    (Appendix-D.tex:666) -- the two were cross-checked and match.

    ``order = 0`` yields ``inf``: zeroth order has no blaze condition.
    """
    m = np.asarray(order, dtype=np.float64)
    zeta = facet_graze(gamma, blaze_angle, alpha)
    numerator = 2.0 * period * np.sin(zeta) * np.sin(blaze_angle)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(m != 0, numerator / m, np.inf)
