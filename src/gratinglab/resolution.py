r"""Resolving power and the finite-:math:`N` line profile.

The infinite grating diffracts each wavelength into delta-function orders; a
real grating illuminates :math:`N` grooves and each order becomes a line with
the interference profile of ISSI eq. (8),

.. math::
    \left[\frac{\sin(Ns)}{N\sin(s)}\right]^2, \qquad
    s = \left[\sin\alpha + \sin\beta\right]\sin\gamma\,\frac{p\pi}{\lambda},

whose first zeros sit at :math:`s = m\pi \pm \pi/N`. The Rayleigh criterion --
one wavelength's peak on the other's first zero -- then gives

.. math::
    R \equiv \frac{\lambda}{\Delta\lambda} = |m|\,N,

so the closed form and the numeric profile are one implementation checking the
other, exactly the pattern the scalar solver uses for its sawtooth closed forms.

**Scope, deliberately staged.** This module answers the ideal-grating question:
:math:`N` identical grooves at exact spacing. The successor -- resolving power
degraded by groove placement and period errors, computed from a measured groove
ensemble -- is roadmap work and is *not* approximated here. An
:class:`~gratinglab.problem.Problem` whose ``n_grooves`` is ``None`` is refused
rather than silently treated as infinite, because "R is undefined" and
"R is very large" are different answers.

``n_grooves`` is the **illuminated** groove count, not the total count on the
part -- :meth:`BoundaryProfile.to_problem` fills it with the number of grooves
the scan actually measured, which is what makes a measured surface flow through
to a spectrograph number in-process.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .geometry import beta, interference_factor, is_propagating, sin_beta
from .illumination import Illumination
from .problem import Problem

__all__ = ["ResolvingPower", "resolving_power", "line_profile", "describe"]


@dataclass(frozen=True, slots=True)
class ResolvingPower:
    """The resolving power of one order at one wavelength.

    Attributes
    ----------
    order, wavelength, n_grooves
        The question asked: order index ``m``, wavelength in nm, illuminated
        groove count ``N``.
    r
        :math:`\\lambda/\\Delta\\lambda` at the stated criterion.
    delta_lambda
        The just-resolved wavelength difference, nm.
    criterion
        How "resolved" is judged. Only ``"rayleigh"`` exists; the field is
        carried so a future criterion (e.g. FWHM-based) is a value, not an
        incompatible reinterpretation of ``r``.
    """

    order: int
    wavelength: float
    n_grooves: int
    r: float
    delta_lambda: float
    criterion: str = "rayleigh"


def _check(problem: Problem, illumination: Illumination, wavelength: float, order: int) -> None:
    """The three refusals shared by every entry point in this module."""
    if problem.n_grooves is None:
        raise ValueError(
            "n_grooves is None: the infinite-grating limit has delta-function "
            "orders and no finite resolving power. Set Problem.n_grooves to "
            "the illuminated groove count -- a measured BoundaryProfile "
            "carries it via to_problem()."
        )
    if order == 0:
        raise ValueError(
            "order 0 does not disperse: every wavelength lands at beta_0 = "
            "-alpha, so its resolving power is undefined."
        )
    s = sin_beta(order, wavelength, problem.period, illumination.sin_alpha,
                 illumination.sin_gamma)
    if not is_propagating(s):
        raise ValueError(
            f"order {order} is evanescent at {wavelength:g} nm "
            f"(sin beta = {float(s):.4g}): there is no propagating line to "
            "resolve."
        )


def resolving_power(
    problem: Problem,
    illumination: Illumination,
    wavelength: float,
    order: int,
) -> ResolvingPower:
    r"""Rayleigh resolving power :math:`R = |m|N` of one propagating order.

    Derivation, so the closed form is checkable against :func:`line_profile`:
    the line of wavelength :math:`\lambda` peaks where
    :math:`\sin\alpha + \sin\beta = m\lambda/(p\sin\gamma)` and first vanishes
    where the right-hand side becomes :math:`(m \pm 1/N)\lambda/(p\sin\gamma)`.
    Placing the peak of :math:`\lambda + \Delta\lambda` on that zero gives
    :math:`m\,\Delta\lambda = \lambda/N`, i.e. :math:`R = |m|N` -- independent
    of the mount, because :math:`\sin\gamma` scales both spacings alike.

    Refuses (``ValueError``) rather than approximates: ``n_grooves`` unset,
    ``order == 0``, or an order that is evanescent at this wavelength.
    """
    _check(problem, illumination, wavelength, order)
    assert problem.n_grooves is not None  # narrowed by _check
    r = abs(order) * problem.n_grooves
    return ResolvingPower(
        order=order,
        wavelength=wavelength,
        n_grooves=problem.n_grooves,
        r=float(r),
        delta_lambda=wavelength / r,
    )


def line_profile(
    problem: Problem,
    illumination: Illumination,
    wavelength: float,
    order: int,
    *,
    n_points: int = 1001,
    half_widths: float = 3.0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    r"""The angular line shape of one order: ``(beta_deg, intensity)``.

    Samples :math:`\beta` around :math:`\beta_m` and evaluates the
    interference function, normalised to 1 at the peak. ``half_widths`` is the
    half-window in units of the first-zero spacing :math:`\pi/N` (in ``s``),
    so the default window always contains the line and its first
    ``half_widths`` zeros regardless of :math:`N`.

    The window is clipped where it would leave :math:`|\sin\beta| \le 1`: a
    line sitting near grazing exit is returned truncated rather than padded
    with angles that do not exist. Same refusals as :func:`resolving_power`.
    """
    _check(problem, illumination, wavelength, order)
    assert problem.n_grooves is not None  # narrowed by _check
    n = problem.n_grooves

    center = float(sin_beta(order, wavelength, problem.period,
                            illumination.sin_alpha, illumination.sin_gamma))
    # ds = pi/N in s corresponds to this spacing in sin(beta).
    zero_spacing = wavelength / (problem.period * illumination.sin_gamma * n)
    lo = max(center - half_widths * zero_spacing, -1.0)
    hi = min(center + half_widths * zero_spacing, 1.0)
    sin_b = np.linspace(lo, hi, n_points)

    s = (illumination.sin_alpha + sin_b) * illumination.sin_gamma * (
        problem.period * np.pi / wavelength
    )
    intensity = interference_factor(s, n)
    beta_deg = np.degrees(beta(sin_b))
    return beta_deg, intensity


def describe(
    problem: Problem,
    illumination: Illumination,
    wavelength: float,
    orders: NDArray[np.int64] | list[int],
) -> str:
    """One human-readable line per order with a defined resolving power.

    Presentation helper for a GUI readout: orders without a defined R (the
    zeroth, and any evanescent at this wavelength) are skipped rather than
    refused, because a readout enumerating orders is a different contract from
    a computation asked about one. Returns ``""`` when nothing qualifies --
    including when ``n_grooves`` is unset, which is the common case for an
    idealised problem.
    """
    if problem.n_grooves is None:
        return ""
    lines = []
    for m in orders:
        try:
            rp = resolving_power(problem, illumination, wavelength, int(m))
        except ValueError:
            continue
        lines.append(
            f"m = {rp.order:+d}: R = {rp.r:,.0f}  "
            f"(N = {rp.n_grooves}, Δλ = {rp.delta_lambda:.3g} nm)"
        )
    return "\n".join(lines)
