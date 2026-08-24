r"""Nystrom discretisation of the boundary integral operators.

The rectangular rule on equally spaced arc-length nodes is spectrally
accurate for periodic integrands (M&P section 4.6.1), which is why
``Profile.boundary`` resamples at equal arc length. Four operators are
built -- the shared vocabulary of the perfectly conducting solve and the
finite-conductivity coupled system (Goray & Schmidt 2010):

**Single layer** ``V`` -- ``\int G(s,s') sigma(s') ds'`` as a dense matrix.
``G`` has an integrable log singularity at ``s' = s``, handled per M&P
section 4.6.3: subtract the periodic comparison function

.. math::
    \tilde G_\infty(u) = \frac{1}{2\pi}
        \left[\ln\frac{u}{L} + \ln\left(1 - \frac{u}{L}\right)\right],
    \qquad u = |s - s'|,

whose bracket with the unknown is continuous (its diagonal limit is
``\ln(2 \pi L / d) / 2\pi`` plus the remainder diagonal of eq. 4.86) and whose
own integral is exactly ``-L / \pi`` (eq. 4.100). Off the diagonal the matrix
entry is simply ``spacing * G``; everything singular lands on the diagonal:

.. math::
    A_{jj} = \Delta s \left[R(0,0) + \frac{\ln(2\pi L/d)}{2\pi}
        - \tilde S\right] - \frac{L}{\pi},
    \qquad
    \tilde S = \frac{1}{\pi}\left[\ln(P-1)! - (P-1)\ln P\right],

with ``\tilde S`` the (node-independent) rectangular-rule sum of
``\tilde G_\infty`` over the off-diagonal nodes.

**Double layer** ``K`` -- ``\int N(s,s') psi(s') ds'`` with ``N = dG/dn'``,
the normal at the *source* point (the sign follows from the double-layer
jump relation for our kernel orientation; the flat-mirror and Table 4.1
tests pin it). The kernel is continuous on a smooth curve, so the plain
rectangular rule applies; the diagonal is the numerical principal-value
limit, taken as the average of the kernel evaluated at the two physically
adjacent nodes (through the period wrap), which converges as O(spacing^2).

**Adjoint double layer** ``L`` -- the same kernel machinery with the normal
moved from source to observation point and the sign flipped
(``dG/dn = -N``-form with ``n`` taken at ``s``), same principal-value
diagonal. It is the principal value of the normal derivative of the
single-layer potential: the jump relation reads
``d(V phi)/dn |_pm = (L pm I/2) phi``, the ``+`` side being the one the
outward normal points into. The finite-conductivity system needs it to
match normal derivatives across the interface.

**Tangential derivative** ``D_t`` -- ``d/ds`` along the boundary of an
``alpha0``-quasi-periodic trace, as a dense differentiation matrix: strip
the Bloch phase ``exp(i alpha0 x(s))``, differentiate the now-periodic
factor spectrally on the uniform arc-length grid, restore the phase (the
product rule contributes ``i alpha0 x'(s)``, and ``x'(s) = n_y`` is the
tangent's x-component for our orientation: nodes ordered by increasing
``x``, outward normal up). Geometry-only -- built once per boundary and
reused at every wavelength. The conical finite-conductivity coupling is
``V (D_t V)`` products, so no new singular quadrature is needed.

The perfectly conducting solve keeps its milestone-1 names as thin
compositions: ``dirichlet_matrix`` (TE) is ``V`` itself and
``neumann_matrix`` (TM) is ``I/2 + K``.
"""

from __future__ import annotations

from math import lgamma

import numpy as np
from numpy.typing import NDArray

from ._boundary import PhysicalBoundary
from ._greens import greens_function, greens_remainder_diagonal, neumann_function

__all__ = [
    "single_layer_matrix",
    "double_layer_matrix",
    "adjoint_layer_matrix",
    "tangential_derivative_matrix",
    "dirichlet_matrix",
    "neumann_matrix",
]


def _separations(
    boundary: PhysicalBoundary,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_]]:
    """Pairwise ``(X, Z)`` separations with the diagonal masked off.

    The diagonal is replaced by a far-away dummy separation so the kernel
    functions never see their singular point; callers overwrite it.
    """
    x, y = boundary.x, boundary.y
    big_x = x[:, None] - x[None, :]
    big_z = y[:, None] - y[None, :]
    diag = np.eye(len(x), dtype=bool)
    big_x[diag] = 0.25 * (x.max() - x.min()) + boundary.spacing
    big_z[diag] = 0.0
    return big_x, big_z, diag


def _wrapped_neighbours(
    boundary: PhysicalBoundary, period: float
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """Separations to the two physically adjacent nodes, through the wrap.

    The wrapped neighbours cross the period boundary, so their x separation
    carries the period explicitly rather than picking up a spurious Bloch
    phase from the kernel's quasi-periodicity.
    """
    x, y = boundary.x, boundary.y
    x_next, y_next = np.roll(x, -1), np.roll(y, -1)
    x_next[-1] += period
    x_prev, y_prev = np.roll(x, 1), np.roll(y, 1)
    x_prev[0] -= period
    return x - x_next, y - y_next, x - x_prev, y - y_prev


def single_layer_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
) -> NDArray[np.complex128]:
    """The discretised single-layer operator ``V``."""
    big_x, big_z, diag = _separations(boundary)
    spacing, length = boundary.spacing, boundary.arc_length
    points = len(boundary.x)

    matrix = spacing * greens_function(
        big_x, big_z, k=k, alpha0=alpha0, period=period, terms=terms
    )

    remainder = greens_remainder_diagonal(
        k=k, alpha0=alpha0, period=period, terms=terms
    )
    comparison_sum = (lgamma(points) - (points - 1) * np.log(points)) / np.pi
    matrix[diag] = (
        spacing
        * (
            remainder
            + np.log(2.0 * np.pi * length / period) / (2.0 * np.pi)
            - comparison_sum
        )
        - length / np.pi
    )
    return matrix


def double_layer_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
) -> NDArray[np.complex128]:
    """The discretised double-layer operator ``K`` (no identity part)."""
    big_x, big_z, diag = _separations(boundary)
    spacing = boundary.spacing
    nx, ny = boundary.nx, boundary.ny

    kernel = neumann_function(
        big_x,
        big_z,
        nx[None, :],
        ny[None, :],
        k=k,
        alpha0=alpha0,
        period=period,
        terms=terms,
    )

    # Diagonal: principal-value limit as the average over the two physically
    # adjacent nodes, whose *source* normals are the neighbours' own.
    dx_next, dy_next, dx_prev, dy_prev = _wrapped_neighbours(boundary, period)
    forward = neumann_function(
        dx_next, dy_next, np.roll(nx, -1), np.roll(ny, -1),
        k=k, alpha0=alpha0, period=period, terms=terms,
    )
    backward = neumann_function(
        dx_prev, dy_prev, np.roll(nx, 1), np.roll(ny, 1),
        k=k, alpha0=alpha0, period=period, terms=terms,
    )
    kernel[diag] = 0.5 * (forward + backward)

    return spacing * kernel


def adjoint_layer_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
) -> NDArray[np.complex128]:
    """The discretised adjoint double-layer operator ``L``.

    The kernel is the normal derivative of ``G`` at the *observation* point:
    ``n(s) . grad G = -``\\ [the ``N``-form with the normal taken at ``s``],
    so the same accelerated kernel serves with the normals broadcast along
    rows instead of columns and the sign flipped. The diagonal is the same
    adjacent-node principal-value average, with the node's own normal.
    """
    big_x, big_z, diag = _separations(boundary)
    spacing = boundary.spacing
    nx, ny = boundary.nx, boundary.ny

    kernel = -neumann_function(
        big_x,
        big_z,
        nx[:, None],
        ny[:, None],
        k=k,
        alpha0=alpha0,
        period=period,
        terms=terms,
    )

    dx_next, dy_next, dx_prev, dy_prev = _wrapped_neighbours(boundary, period)
    forward = -neumann_function(
        dx_next, dy_next, nx, ny, k=k, alpha0=alpha0, period=period, terms=terms
    )
    backward = -neumann_function(
        dx_prev, dy_prev, nx, ny, k=k, alpha0=alpha0, period=period, terms=terms
    )
    kernel[diag] = 0.5 * (forward + backward)

    return spacing * kernel


def tangential_derivative_matrix(
    boundary: PhysicalBoundary, *, alpha0: float
) -> NDArray[np.complex128]:
    """The discretised tangential derivative ``D_t`` for quasi-periodic traces.

    Applied to node values of an ``alpha0``-quasi-periodic function of arc
    length (such as the trace of a single-layer potential), returns ``d/ds``:
    the Bloch phase is stripped, the periodic factor differentiated
    spectrally on the uniform arc-length grid, the phase restored, and the
    product-rule term ``i alpha0 x'(s)`` added with ``x'(s) = n_y``.
    Geometry-only: no wavelength enters.
    """
    points = len(boundary.x)
    wavenumber = 2j * np.pi * np.fft.fftfreq(points, d=boundary.spacing)
    if points % 2 == 0:
        # The unpaired Nyquist mode has no odd (derivative) representation.
        wavenumber[points // 2] = 0.0
    spectral = np.fft.ifft(
        wavenumber[:, None] * np.fft.fft(np.eye(points), axis=0), axis=0
    ).real
    phase = np.exp(1j * alpha0 * boundary.x)
    matrix = (phase[:, None] / phase[None, :]) * spectral
    matrix[np.diag_indices(points)] += 1j * alpha0 * boundary.ny
    return matrix


def dirichlet_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
) -> NDArray[np.complex128]:
    """The TE (Dirichlet) system matrix -- the single-layer operator itself."""
    return single_layer_matrix(
        boundary, k=k, alpha0=alpha0, period=period, terms=terms
    )


def neumann_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
) -> NDArray[np.complex128]:
    """``I/2 + K`` for the TM (Neumann) problem."""
    return 0.5 * np.eye(len(boundary.x)) + double_layer_matrix(
        boundary, k=k, alpha0=alpha0, period=period, terms=terms
    )
