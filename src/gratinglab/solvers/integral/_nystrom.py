r"""Nystrom discretisation of the boundary integral operators.

The rectangular rule on equally spaced arc-length nodes is spectrally
accurate for periodic integrands (M&P section 4.6.1), which is why
``Profile.boundary`` resamples at equal arc length. Two operators are built:

**Dirichlet (TE)** -- the single-layer operator ``\int G(s,s') sigma(s') ds'``
as a dense matrix. ``G`` has an integrable log singularity at ``s' = s``,
handled per M&P section 4.6.3: subtract the periodic comparison function

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

**Neumann (TM)** -- ``psi/2 + \int N(s,s') psi(s') ds' = psi_i`` with
``N = dG/dn'`` (the sign follows from the double-layer jump relation for our
kernel orientation; the flat-mirror and Table 4.1 tests pin it). The kernel
is continuous on a smooth curve, so the plain
rectangular rule applies; the diagonal is the numerical principal-value
limit, taken as the average of the kernel evaluated at the two physically
adjacent nodes (through the period wrap), which converges as O(spacing^2).
"""

from __future__ import annotations

from math import lgamma

import numpy as np
from numpy.typing import NDArray

from ._boundary import PhysicalBoundary
from ._greens import greens_function, greens_remainder_diagonal, neumann_function

__all__ = ["dirichlet_matrix", "neumann_matrix"]


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


def dirichlet_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
) -> NDArray[np.complex128]:
    """The discretised single-layer operator for the TE (Dirichlet) problem."""
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


def neumann_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
) -> NDArray[np.complex128]:
    """``I/2 + spacing * N`` for the TM (Neumann) problem."""
    big_x, big_z, diag = _separations(boundary)
    spacing = boundary.spacing
    x, y, nx, ny = boundary.x, boundary.y, boundary.nx, boundary.ny

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
    # adjacent nodes. The wrapped neighbours cross the period boundary, so
    # their x separation carries the period explicitly rather than picking up
    # a spurious Bloch phase from the kernel's quasi-periodicity.
    x_next, y_next = np.roll(x, -1), np.roll(y, -1)
    x_next[-1] += period
    x_prev, y_prev = np.roll(x, 1), np.roll(y, 1)
    x_prev[0] -= period
    forward = neumann_function(
        x - x_next, y - y_next, np.roll(nx, -1), np.roll(ny, -1),
        k=k, alpha0=alpha0, period=period, terms=terms,
    )
    backward = neumann_function(
        x - x_prev, y - y_prev, np.roll(nx, 1), np.roll(ny, 1),
        k=k, alpha0=alpha0, period=period, terms=terms,
    )
    kernel[diag] = 0.5 * (forward + backward)

    return 0.5 * np.eye(len(x)) + spacing * kernel
