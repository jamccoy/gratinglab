r"""Separable (GEMM) assembly of the accelerated kernels.

``_greens`` evaluates the Kummer-accelerated kernels pointwise, which for a
full ``(P, P)`` operator means broadcast ``(P, P, orders)`` intermediates hit
with ``exp`` and complex powers -- measured at 97-99% of the whole solver's
runtime, dwarfing the dense LU. This module assembles the *same* sums as a
handful of BLAS-3 products, exploiting that every spectral term separates
into node factors:

.. math::
    e^{i \alpha_m (x_j - x_l)} = e^{i \alpha_m x_j}\, e^{-i \alpha_m x_l},
    \qquad
    e^{i \gamma_m |y_j - y_l|} =
      \begin{cases}
        e^{i \gamma_m (y_j - c)}\, e^{-i \gamma_m (y_l - c)} & y_j \ge y_l \\
        e^{-i \gamma_m (y_j - c)}\, e^{+i \gamma_m (y_l - c)} & y_j < y_l
      \end{cases}

so the truncated sum over orders ``-M..M`` is two ``(P, 2M+1) @ (2M+1, P)``
matrix products -- one per sign branch of ``Z = y_j - y_l``, selected
entrywise (the height is not monotonic in node index, so the branches form
no triangle; computing both and selecting doubles the flops but keeps two
clean GEMMs, still three orders of magnitude below the broadcast form).

Three safeguards make the split exact rather than merely fast:

- **The reference shift** ``c = (y_min + y_max) / 2`` bounds every factor's
  exponent by ``Im(gamma_m) h / 2`` with ``h`` the groove height: the split
  factors individually can overflow where their product cannot, and the
  midpoint shift is the best any single reference can do.
- **The hard-order fallback**: orders with ``Im(gamma_m) h / 2 >
  _EXP_GUARD`` (safely under ``log`` of the double-precision max) are
  evaluated by the original broadcast form restricted to those orders,
  which exponentiates the pair *difference* and never overflows. The set is
  empty for realistic meshes (it needs roughly ``pi * terms * h / period >
  1380``) and deep gratings degrade gracefully instead of returning
  ``0 * inf``.
- **Geometry-cached Kummer tails**: the asymptote ratios ``xi_pm =
  exp(pm K(iX mp |Z|))`` carry no wavelength (only the prefactors
  ``pre_pm = exp(alpha0(iX mp |Z|))`` do), so the accelerated series'
  entire asymptotic content collapses to two per-boundary arrays each:
  ``sum_{p>M} xi^p / p`` for the Green's function (via ``log1p`` minus a
  Horner partial sum -- no powers, ``|xi| <= 1`` throughout) and the
  geometric ``sum_{p>M} xi^p = xi^{M+1} / (1 - xi)`` for the Neumann
  kernel. Built once per (boundary, terms), reused at every wavelength.

The Neumann kernel is served as its two node-factorable pieces,

.. math::
    N(s, s') = A\, v_x(s') + \mathrm{sgn}(Z)\, B\, v_y(s'), \qquad
    A = -\frac{1}{2d} \sum_m \frac{\alpha_m}{\gamma_m}
        e^{i\alpha_m X + i\gamma_m |Z|},\quad
    B = -\frac{1}{2d} \sum_m e^{i\alpha_m X + i\gamma_m |Z|},

so one evaluation serves the double layer (vector = source normal, columns),
the adjoint layer (observation normal, rows, sign flipped) and the
tangential layer (observation tangent, rows, sign flipped).

Diagonals are garbage here by design -- the quadrature layer (``_nystrom``)
overwrites them with the analytic and principal-value limits, exactly as it
does for the pointwise kernels, and the cached tails are zeroed on the
diagonal so no ``inf`` leaks into the elementwise combines.

Shared formulation envelope (not introduced here): the Kummer acceleration
-- pointwise and separable alike -- cancels asymptote terms of magnitude
``exp(|alpha0| |Z|)`` down to O(1), so past ``|alpha0| h ~ 30`` double
precision retains nothing. Measured against a cancellation-free reference,
this assembly tracks the pointwise kernels to ``~exp(|alpha0| h) * 1e-16``
and is, if anything, the slightly *more* accurate of the two throughout the
degradation window. The regime needs a short transverse wavelength driven
through a deep groove (an in-plane mount at hard-X-ray wavelengths); the
grazing-conical mounts the solver targets keep ``alpha0`` small, and the
energy-balance theorem check fails loudly if the envelope is ever left.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ._boundary import PhysicalBoundary
from ._greens import wavenumbers

__all__ = [
    "KernelGeometry",
    "SpectralFactors",
    "kernel_geometry",
    "spectral_factors",
    "greens_matrix",
    "neumann_core",
]

#: Largest ``Im(gamma_m) * height / 2`` a split factor may carry -- safely
#: under ``log(DBL_MAX) ~ 709.78``. Orders beyond it take the fallback.
_EXP_GUARD = 690.0

#: Fallback orders are accumulated in chunks of this many, matching
#: ``_greens._CHUNK`` so the broadcast intermediates stay small.
_CHUNK = 32


@dataclass(frozen=True, slots=True)
class KernelGeometry:
    """Everything about one (boundary, terms) pair that no wavelength touches.

    Roughly eleven ``(P, P)`` arrays -- ~20 MB at ``P = 400``, ~500 MB at
    ``P = 2000`` -- built once per scan and amortised over every wavelength
    and both media.
    """

    period: float
    terms: int
    x: NDArray[np.float64]
    y_shifted: NDArray[np.float64]  #: ``y - y_center``
    height: float
    big_x: NDArray[np.float64]  #: true separations, no diagonal dummy
    absz: NDArray[np.float64]
    upper: NDArray[np.bool_]  #: ``y_j >= y_l``
    sgn: NDArray[np.float64]
    up_dir: NDArray[np.complex128]  #: ``iX - |Z|``
    down_dir: NDArray[np.complex128]  #: ``iX + |Z|``
    greens_tail_p: NDArray[np.complex128]  #: ``sum_{p>M} xi_+^p / p``, diag 0
    greens_tail_m: NDArray[np.complex128]
    neumann_tail_p: NDArray[np.complex128]  #: ``sum_{p>M} xi_+^p``, diag 0
    neumann_tail_m: NDArray[np.complex128]


@dataclass(frozen=True, slots=True)
class SpectralFactors:
    """The per-(k, alpha0) node factors of the spectral sum.

    ``row/col`` matrices cover the *safe* orders only; the ``hard`` arrays
    hold the wavenumbers of orders whose split factors would overflow, for
    the broadcast fallback in :func:`_spectral_sum`.
    """

    k: complex
    alpha0: float
    alpha_safe: NDArray[np.float64]
    gamma_safe: NDArray[np.complex128]
    alpha_hard: NDArray[np.float64]
    gamma_hard: NDArray[np.complex128]
    row_up: NDArray[np.complex128]  #: ``(P, S) exp(i a x_j + i g (y_j - c))``
    col_up: NDArray[np.complex128]  #: ``(S, P) exp(-i a x_l - i g (y_l - c))``
    row_dn: NDArray[np.complex128]
    col_dn: NDArray[np.complex128]
    pre_p: NDArray[np.complex128]  #: ``exp(alpha0 (iX - |Z|))``
    pre_m: NDArray[np.complex128]


def kernel_geometry(
    boundary: PhysicalBoundary, *, period: float, terms: int
) -> KernelGeometry:
    """Build the wavelength-independent half of the kernel assembly."""
    x, y = boundary.x, boundary.y
    big_x = x[:, None] - x[None, :]
    big_z = y[:, None] - y[None, :]
    absz = np.abs(big_z)
    up_dir = 1j * big_x - absz
    down_dir = 1j * big_x + absz
    diag = np.eye(len(x), dtype=bool)
    big_k = 2.0 * np.pi / period

    def tails(xi: NDArray[np.complex128], direction: NDArray[np.complex128]):
        # sum_{p<=M} xi^p / p by Horner: no powers, |xi| <= 1 throughout.
        partial = np.zeros_like(xi)
        for p in range(terms, 0, -1):
            partial = xi * (1.0 / p + partial)
        with np.errstate(divide="ignore", invalid="ignore"):
            log_tail = -np.log1p(-xi) - partial
            geometric_tail = np.exp((terms + 1) * direction) / (1.0 - xi)
        # xi = 1 exactly on the diagonal (and nowhere else within one
        # period); zero it so no inf reaches the elementwise combines --
        # the quadrature layer overwrites the diagonal regardless.
        log_tail[diag] = 0.0
        geometric_tail[diag] = 0.0
        return log_tail, geometric_tail

    greens_p, neumann_p = tails(np.exp(big_k * up_dir), big_k * up_dir)
    greens_m, neumann_m = tails(np.exp(-big_k * down_dir), -big_k * down_dir)

    y_min, y_max = float(y.min()), float(y.max())
    return KernelGeometry(
        period=period,
        terms=terms,
        x=x,
        y_shifted=y - 0.5 * (y_min + y_max),
        height=y_max - y_min,
        big_x=big_x,
        absz=absz,
        upper=big_z >= 0.0,
        sgn=np.sign(big_z),
        up_dir=up_dir,
        down_dir=down_dir,
        greens_tail_p=greens_p,
        greens_tail_m=greens_m,
        neumann_tail_p=neumann_p,
        neumann_tail_m=neumann_m,
    )


def spectral_factors(
    geometry: KernelGeometry, *, k: complex, alpha0: float
) -> SpectralFactors:
    """Node factors for one medium at one wavelength.

    One :func:`~._greens.wavenumbers` call covers all orders ``-M..M`` --
    the branch pin lives there and only there.
    """
    orders = np.arange(-geometry.terms, geometry.terms + 1)
    alpha_m, gamma_m = wavenumbers(k, alpha0, geometry.period, orders)
    safe = gamma_m.imag * (0.5 * geometry.height) <= _EXP_GUARD

    a = alpha_m[safe]
    g = gamma_m[safe]
    x = geometry.x
    ys = geometry.y_shifted
    return SpectralFactors(
        k=complex(k),
        alpha0=float(alpha0),
        alpha_safe=a,
        gamma_safe=g,
        alpha_hard=alpha_m[~safe],
        gamma_hard=gamma_m[~safe],
        row_up=np.exp(1j * a[None, :] * x[:, None] + 1j * g[None, :] * ys[:, None]),
        col_up=np.exp(-1j * a[:, None] * x[None, :] - 1j * g[:, None] * ys[None, :]),
        row_dn=np.exp(1j * a[None, :] * x[:, None] - 1j * g[None, :] * ys[:, None]),
        col_dn=np.exp(-1j * a[:, None] * x[None, :] + 1j * g[:, None] * ys[None, :]),
        pre_p=np.exp(alpha0 * geometry.up_dir),
        pre_m=np.exp(alpha0 * geometry.down_dir),
    )


def _spectral_sum(geometry, factors, weight) -> NDArray[np.complex128]:
    """``sum_m weight(alpha_m, gamma_m) exp(i alpha_m X + i gamma_m |Z|)``.

    Safe orders as two branch GEMMs selected on ``sign(Z)``; the discarded
    branch may overflow harmlessly (each selected entry's own accumulation
    only ever decays), hence the suppressed warnings. Hard orders take the
    chunked broadcast form on the pair differences, which cannot overflow.
    """
    w = weight(factors.alpha_safe, factors.gamma_safe)
    with np.errstate(over="ignore", invalid="ignore"):
        up = factors.row_up @ (w[:, None] * factors.col_up)
        down = factors.row_dn @ (w[:, None] * factors.col_dn)
        total = np.where(geometry.upper, up, down)
    for start in range(0, len(factors.alpha_hard), _CHUNK):
        a = factors.alpha_hard[start : start + _CHUNK]
        g = factors.gamma_hard[start : start + _CHUNK]
        phases = np.exp(
            1j * a * geometry.big_x[..., None]
            + 1j * g * geometry.absz[..., None]
        )
        total = total + (weight(a, g) * phases).sum(axis=-1)
    return total


def greens_matrix(
    geometry: KernelGeometry, factors: SpectralFactors
) -> NDArray[np.complex128]:
    """The full ``(P, P)`` accelerated Green's function, diagonal garbage.

    Algebraically identical to :func:`~._greens.greens_function` on the true
    separations: the truncated spectral sum minus the tail of the asymptote
    series (the per-order pairing ``exact - asym`` plus the closed-form logs
    telescopes to exactly this).
    """
    period = geometry.period
    total = _spectral_sum(
        geometry, factors, lambda a, g: 1.0 / (2j * period * g)
    )
    return total - (
        factors.pre_p * geometry.greens_tail_p
        + factors.pre_m * geometry.greens_tail_m
    ) / (4.0 * np.pi)


def neumann_core(
    geometry: KernelGeometry, factors: SpectralFactors
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """The Neumann kernel's two vector components ``(A, B)``.

    ``N = A v_x + sgn(Z) B v_y`` for any unit vector ``v`` at source or
    observation point -- the caller broadcasts ``v`` along columns (double
    layer) or rows (adjoint and tangential layers, with the sign flip).
    The tail terms carry the large-``|m|`` bracket limits ``alpha_m /
    gamma_m -> mp i`` exactly as the pointwise kernel's asymptote does.
    """
    period = geometry.period
    a_matrix = _spectral_sum(
        geometry, factors, lambda a, g: -a / (2.0 * period * g)
    )
    b_matrix = _spectral_sum(
        geometry,
        factors,
        lambda a, g: np.full(g.shape, -1.0 / (2.0 * period), dtype=np.complex128),
    )
    a_matrix = a_matrix + (
        -1j * factors.pre_p * geometry.neumann_tail_p
        + 1j * factors.pre_m * geometry.neumann_tail_m
    ) / (-2.0 * period)
    b_matrix = b_matrix + (
        factors.pre_p * geometry.neumann_tail_p
        + factors.pre_m * geometry.neumann_tail_m
    ) / (-2.0 * period)
    return a_matrix, b_matrix
