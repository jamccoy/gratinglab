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
    "flux_obliquity",
    "facet_graze",
    "sin_facet_graze",
    "horizon_visible",
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


def flux_obliquity(cos_alpha: ArrayLike, cos_beta_m: ArrayLike) -> NDArray[np.float64]:
    r"""Flux projection factor for a diffracted order.

    .. math::
        O_m = \frac{4\cos\alpha\,\cos\beta_m}{(\cos\alpha + \cos\beta_m)^2}

    Multiplying :math:`|G_m|^2` by this is what makes scalar theory agree with
    **first-order Rayleigh perturbation theory** in the shallow-groove limit,
    where both are valid and must therefore give the same answer. Perturbation
    theory says :math:`\eta_m = 4k_{z,0}k_{z,m}|\hat{g}_m|^2`; the bare Fourier
    coefficient gives :math:`k_{z,0}` and :math:`k_{z,m}` only through their
    *sum*, and this is exactly the ratio between the two. See
    ``tests/test_perturbation.py``, which derives it from two closed forms with
    no solver in the loop.

    Three properties, all load-bearing:

    - **Symmetric** under :math:`\alpha \leftrightarrow \beta_m`, so Lorentz
      reciprocity survives it untouched. That is what makes it admissible at
      all: the obliquity factor of thesis Appendix D,
      :math:`\cos\beta_m/\cos\alpha`, is the *asymmetric* one and breaks
      reciprocity, which is why ``docs/theory/scalar.md`` rejected it. The
      rejection was right; the conclusion that no flux factor belongs was not.
    - **At most 1**, by AM-GM, with equality exactly when
      :math:`\cos\alpha = \cos\beta_m` -- Littrow, specular, and every
      :math:`m = 0` (where :math:`\beta_0 = -\alpha`). So it can only reduce an
      efficiency, and it leaves the shallow-groove energy identity alone.
    - **Vanishes at grazing exit**, :math:`\cos\beta_m \to 0`. An order about to
      pass off carries no flux through a plane parallel to the surface, and the
      unfactored :math:`|G_m|^2` claimed it did.

    Both arguments are cosines, not angles -- the callers already have them from
    :func:`cos_beta` and :attr:`Illumination.cos_alpha`, and converting to
    angles and back would be a chance to lose a sign.
    """
    cos_alpha = np.asarray(cos_alpha, dtype=np.float64)
    cos_beta_m = np.asarray(cos_beta_m, dtype=np.float64)
    return 4.0 * cos_alpha * cos_beta_m / (cos_alpha + cos_beta_m) ** 2


def sin_facet_graze(
    gamma: float, tilt: ArrayLike, angle: ArrayLike
) -> NDArray[np.float64]:
    r"""":math:`\sin\zeta = \sin\gamma\,\cos(\text{tilt} - \text{angle})`, unclipped.

    The vector form of :func:`facet_graze`, and deliberately **not** wrapped in
    ``arcsin``. Two reasons:

    - A **negative** value is the meaningful signal that the facet is turned
      away from the direction in question -- back-facing, and contributing
      nothing. ``arcsin`` of a negative number is a perfectly good angle and
      would hide that; clipping to zero would hide it too. A caller resolving
      reflectivity across a groove tests this sign to build its visibility mask.
    - ``tilt`` varies point by point along a groove, so this is called on whole
      arrays where :func:`facet_graze` takes scalars.

    ``angle`` is :math:`\alpha` for the incident direction and :math:`\beta_m`
    for a diffracted one. It is spelled generically because the two uses are the
    same formula: reflection off a tilted facet does not care which way the
    light is going, which is the symmetry that keeps a groove-resolved
    reflectivity model reciprocal.
    """
    tilt = np.asarray(tilt, dtype=np.float64)
    angle = np.asarray(angle, dtype=np.float64)
    return np.sin(gamma) * np.cos(tilt - angle)


def horizon_visible(
    height_nm: ArrayLike, period: float, angle: float
) -> NDArray[np.bool_]:
    r"""Which points on the groove a ray at ``angle`` can actually reach.

    The mask for **cast** shadows -- the ones the groove apex throws across
    the trough onto surface that faces the ray perfectly well. The local
    orientation test (:func:`sin_facet_graze` ``> 0``) sees only
    *self*-shadowing, a facet turned away from the ray; on the reference
    sawtooth at its working angles the two differ by a fraction of a percent
    of the period on the incident side and by 10-50% on the exit side, which
    is why this function exists (``docs/findings.md``).

    **The half-cone angle is absent, deliberately.** The surface is invariant
    along the grooves, so a 3D ray occludes exactly as its projection into the
    transverse plane does, and ``angle`` is the transverse azimuth --
    :math:`\alpha` for the incident direction, :math:`\beta_m` for a
    diffracted one, measured from the normal. :math:`\gamma` scales the ray's
    groove-parallel component, which slides *along* the invariant direction
    and can neither create nor remove an occlusion.

    **The travel direction is the trap.** The profile parameter runs against
    the dispersion direction (``docs/findings.md``, "The profile parameter
    runs backwards"; conventions section 3), so a ray with ``angle > 0``
    travels toward **+t** and its occluders sit at *smaller* ``t``. Getting
    this backwards does not fail loudly -- it reports shadows on the wrong
    facet, which on a sawtooth mostly coincides with the self-shadowing mask
    and looks plausible. The anchor that pins it: on an ideal sawtooth with
    blaze slope :math:`s_b`, anti-blaze slope :math:`s_a` and apex at
    :math:`t_a`, a ray at ``angle`` :math:`> \operatorname{arccot} s_a`
    (19.5 deg for the 29.5/70.5 reference groove) throws a shadow of width

    .. math::
        \Delta = (1 - t_a)\,\frac{s_a - s_r}{s_b + s_r}, \qquad
        s_r = \cot(\text{angle})

    past the trough onto the blaze facet, and ``tests/test_scalar.py``
    derives that independently and checks it here.

    The scan itself is the classic running-horizon argument, exact and
    :math:`O(n)`: with :math:`u(t) = g(t) + p\,\cot|\theta|\,t` (sign of the
    linear term following the travel direction), a point is lit iff
    :math:`u` at that point matches the running maximum of :math:`u` over
    everything upstream of it. Upstream whole periods enter in closed form --
    each sits exactly :math:`p\cot|\theta|` lower than the last, so only the
    nearest one can ever compete -- which is what removes any need to tile.

    Parameters
    ----------
    height_nm
        Groove height in nm on a **uniform** grid over one period,
        ``t = arange(n)/n``. Uniformity is load-bearing: the linear term is
        built from the index.
    period
        Groove period in nm -- the same unit as ``height_nm``, because the
        ray's run-to-drop ratio compares the two directly.
    angle
        Transverse azimuth of the ray in radians, from the grating normal.
        Must be a propagating direction, :math:`|\theta| \le \pi/2`; exactly
        0 is straight down and shadows nothing.
    """
    height_nm = np.asarray(height_nm, dtype=np.float64)
    n = len(height_nm)
    if angle == 0.0:
        return np.ones(n, dtype=bool)

    t = np.arange(n) / n
    cot = 1.0 / np.tan(abs(angle))
    if angle > 0.0:
        u = height_nm + cot * t * period
        upstream = np.maximum.accumulate(u)
    else:
        # The mirror case: travel toward -t, occluders at larger t, and the
        # linear term flips sign with the run direction.
        u = height_nm - cot * t * period
        upstream = np.maximum.accumulate(u[::-1])[::-1]
    # Every complete upstream period is a copy of u shifted down by one
    # period's drop, so the nearest bounds them all.
    horizon = np.maximum(upstream, u.max() - cot * period)
    return u >= horizon


def facet_graze(gamma: float, blaze_angle: float, alpha: float) -> float:
    r"""Graze angle onto a sawtooth groove facet, in radians.

    .. math::
        \sin\zeta = \sin\gamma \, \cos(\delta - \alpha)

    This is the angle that must stay below the critical angle for total
    external reflection, :math:`\theta_c \approx \sqrt{2\delta_n}`.

    The scalar case of :func:`sin_facet_graze`, and defined in terms of it so
    there is one formula rather than two that have to be kept agreeing.
    """
    return float(np.arcsin(sin_facet_graze(gamma, blaze_angle, alpha)))


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
