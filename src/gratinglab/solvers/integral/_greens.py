r"""Quasi-periodic Green's function and its normal-derivative kernel.

Everything here lives in the transverse plane of the conical reduction:
``x`` along the dispersion direction, ``y`` along the grating normal, both in
nm, with the transverse wavenumber ``k = 2 pi sin(gamma) / wavelength``. The
Green's function for a grating of period ``d`` illuminated with incident
x-wavenumber ``alpha0`` is the spectral (Rayleigh) sum

.. math::
    G(X, Z) = \frac{1}{2 i d} \sum_m \frac{1}{\gamma_m}
        e^{i \alpha_m X + i \gamma_m |Z|},
    \qquad \alpha_m = \alpha_0 + m K,\quad
    \gamma_m = \sqrt{k^2 - \alpha_m^2}

(Maystre & Popov eq. 4.76, in pseudo-periodic form), with ``K = 2 pi / d``
and the branch ``Im(gamma_m) >= 0`` so evanescent terms decay upward.

The raw sum converges like ``1/m`` -- uselessly slowly near ``X = Z = 0``.
Following M&P section 4.6.2 we subtract the large-``|m|`` asymptote, whose sum
is a closed-form logarithm (the Kummer trick, eqs. 4.82-4.85):

.. math::
    G_\infty = \frac{1}{4\pi}\left[
        e^{\alpha_0 (iX - |Z|)} \ln(1 - e^{K(iX - |Z|)})
      + e^{\alpha_0 (iX + |Z|)} \ln(1 - e^{-K(iX + |Z|)}) \right]

leaving a remainder series whose +m/-m pairs decay like ``m^{-3}`` on the
diagonal (eq. 4.87) and at worst ``m^{-2}`` elsewhere. The remaining
``ln`` singularity of ``G_\infty`` at coincident points is integrable and is
handled by the quadrature layer (``_nystrom``), not here.

The Neumann kernel is the normal derivative at the source point,

.. math::
    N(s, s') = \frac{\partial G}{\partial n'}
    = -\frac{1}{2d} \sum_m e^{i\alpha_m X + i\gamma_m |Z|}
      \left[\frac{\alpha_m}{\gamma_m} n'_x + \mathrm{sgn}(Z)\, n'_y\right],

derived directly from the spectral form (M&P eq. 4.77 states the same kernel;
the derivation here is independent so the signs rest on this module's tests,
not on transcription). Its asymptote sums to closed-form geometric terms
(eq. 4.91) and the paired remainder again decays fast.

``k`` may be complex: the same code serves the metal-side Green's function of
the planned finite-conductivity milestone (Goray & Schmidt 2010), where the
transverse wavenumber carries the complex refractive index.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "wavenumbers",
    "greens_function",
    "greens_remainder_diagonal",
    "neumann_function",
]

#: Terms are accumulated in chunks of this many spectral orders so the
#: broadcast ``(points, points, orders)`` intermediates stay small.
_CHUNK = 32


def wavenumbers(
    k: complex, alpha0: float, period: float, orders: NDArray[np.int64]
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    r"""``(alpha_m, gamma_m)`` for the given spectral orders.

    ``gamma_m = sqrt(k^2 - alpha_m^2)`` on the branch with ``Im >= 0``:
    positive real for propagating orders, positive imaginary for evanescent
    ones, so every term decays or radiates upward.
    """
    alpha_m = alpha0 + 2.0 * np.pi * np.asarray(orders, dtype=np.float64) / period
    gamma_m = np.sqrt(k**2 - alpha_m.astype(np.complex128) ** 2)
    # numpy's sqrt of a negative real returns +i tau, but a complex k (the
    # finite-conductivity seam) can land on the other sheet; pin the branch.
    flip = gamma_m.imag < 0.0
    if flip.any():
        gamma_m = np.where(flip, -gamma_m, gamma_m)
    return alpha_m, gamma_m


def _asymptote_factors(
    x: NDArray[np.float64], z: NDArray[np.float64], alpha0: float, period: float
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """Prefactors and ratios of the large-``|m|`` asymptotic terms.

    Returns ``(pre_plus, xi_plus, pre_minus, xi_minus)`` such that the
    asymptotic Green's-function term for order ``+p`` is
    ``-(pre_plus / 4 pi p) xi_plus**p`` and for ``-p`` is
    ``-(pre_minus / 4 pi p) xi_minus**p``.
    """
    big_k = 2.0 * np.pi / period
    absz = np.abs(z)
    up = 1j * x - absz  # exponent direction for the +m branch
    down = 1j * x + absz
    return (
        np.exp(alpha0 * up),
        np.exp(big_k * up),
        np.exp(alpha0 * down),
        np.exp(-big_k * down),
    )


def greens_function(
    x: NDArray[np.float64],
    z: NDArray[np.float64],
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
) -> NDArray[np.complex128]:
    """Accelerated ``G(X, Z)`` at separations ``x = X``, ``z = Z`` (nm).

    ``terms`` is the spectral truncation ``M``: the remainder series keeps
    orders ``-M..M``; the asymptotic logs carry the rest to infinity in
    closed form. Coincident points (``X = Z = 0`` modulo the period) are the
    kernel's genuine logarithmic singularity -- the quadrature layer supplies
    the diagonal, so this function is only called off-diagonal.
    """
    x = np.asarray(x, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    absz = np.abs(z)
    pre_p, xi_p, pre_m, xi_m = _asymptote_factors(x, z, alpha0, period)

    # m = 0 exactly.
    a0, g0 = wavenumbers(k, alpha0, period, np.array([0]))
    total = np.exp(1j * a0[0] * x + 1j * g0[0] * absz) / (2j * period * g0[0])

    for start in range(1, terms + 1, _CHUNK):
        p = np.arange(start, min(start + _CHUNK, terms + 1))
        ap, gp = wavenumbers(k, alpha0, period, p)
        am, gm = wavenumbers(k, alpha0, period, -p)
        shape = x[..., None]
        exact = np.exp(1j * ap * shape + 1j * gp * absz[..., None]) / (
            2j * period * gp
        ) + np.exp(1j * am * shape + 1j * gm * absz[..., None]) / (2j * period * gm)
        asym = (
            pre_p[..., None] * xi_p[..., None] ** p
            + pre_m[..., None] * xi_m[..., None] ** p
        ) / (-4.0 * np.pi * p)
        total = total + (exact - asym).sum(axis=-1)

    g_inf = (np.log1p(-xi_p) * pre_p + np.log1p(-xi_m) * pre_m) / (4.0 * np.pi)
    return total + g_inf


def greens_remainder_diagonal(
    *, k: complex, alpha0: float, period: float, terms: int
) -> complex:
    """``lim (G - G_inf)`` at coincident points (M&P eq. 4.86).

    The paired terms decay like ``m^-3`` (eq. 4.87), so modest ``terms``
    suffice. The quadrature layer adds the integrated log singularity.
    """
    p = np.arange(1, terms + 1)
    _, g0 = wavenumbers(k, alpha0, period, np.array([0]))
    _, gp = wavenumbers(k, alpha0, period, p)
    _, gm = wavenumbers(k, alpha0, period, -p)
    paired = 1.0 / (2j * period * gp) + 1.0 / (2j * period * gm) + 1.0 / (
        2.0 * np.pi * p
    )
    return complex(1.0 / (2j * period * g0[0]) + paired.sum())


def neumann_function(
    x: NDArray[np.float64],
    z: NDArray[np.float64],
    source_nx: NDArray[np.float64],
    source_ny: NDArray[np.float64],
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
) -> NDArray[np.complex128]:
    """Accelerated ``N(s, s') = dG/dn'`` off-diagonal.

    ``source_nx, source_ny`` is the outward unit normal at the *source* point
    ``s'``; ``x = x(s) - x(s')``, ``z = y(s) - y(s')``. The kernel is
    continuous through coincident points on a smooth curve (the limit picks up
    the local curvature, M&P eq. 4.92), but this function is only valid off
    the diagonal -- the quadrature layer supplies the diagonal as the
    numerical principal-value limit, by averaging the kernel at the two
    physically adjacent nodes. That estimate was measured to converge as
    O(spacing^2) against the exact limit, which keeps it below the spectral
    truncation error at any realistic node count.
    """
    x = np.asarray(x, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    source_nx = np.broadcast_to(np.asarray(source_nx, dtype=np.float64), x.shape)
    source_ny = np.broadcast_to(np.asarray(source_ny, dtype=np.float64), x.shape)
    absz = np.abs(z)
    sgn = np.sign(z)
    pre_p, xi_p, pre_m, xi_m = _asymptote_factors(x, z, alpha0, period)

    a0, g0 = wavenumbers(k, alpha0, period, np.array([0]))
    total = (
        -np.exp(1j * a0[0] * x + 1j * g0[0] * absz)
        * ((a0[0] / g0[0]) * source_nx + sgn * source_ny)
        / (2.0 * period)
    )

    for start in range(1, terms + 1, _CHUNK):
        p = np.arange(start, min(start + _CHUNK, terms + 1))
        ap, gp = wavenumbers(k, alpha0, period, p)
        am, gm = wavenumbers(k, alpha0, period, -p)
        bracket_p = (ap / gp) * source_nx[..., None] + (sgn * source_ny)[..., None]
        bracket_m = (am / gm) * source_nx[..., None] + (sgn * source_ny)[..., None]
        exact = (
            np.exp(1j * ap * x[..., None] + 1j * gp * absz[..., None]) * bracket_p
            + np.exp(1j * am * x[..., None] + 1j * gm * absz[..., None]) * bracket_m
        ) / (-2.0 * period)
        # Large-|m| limits alpha_m/gamma_m -> -i (m>0) and +i (m<0).
        asym_bracket_p = -1j * source_nx[..., None] + (sgn * source_ny)[..., None]
        asym_bracket_m = 1j * source_nx[..., None] + (sgn * source_ny)[..., None]
        asym = (
            pre_p[..., None] * xi_p[..., None] ** p * asym_bracket_p
            + pre_m[..., None] * xi_m[..., None] ** p * asym_bracket_m
        ) / (-2.0 * period)
        total = total + (exact - asym).sum(axis=-1)

    # Geometric closed forms xi/(1-xi) for the asymptotic branches
    # (M&P eq. 4.91).
    geo_p = xi_p / (1.0 - xi_p)
    geo_m = xi_m / (1.0 - xi_m)
    n_inf = (
        pre_p * geo_p * (-1j * source_nx + sgn * source_ny)
        + pre_m * geo_m * (1j * source_nx + sgn * source_ny)
    ) / (-2.0 * period)
    return total + n_inf
