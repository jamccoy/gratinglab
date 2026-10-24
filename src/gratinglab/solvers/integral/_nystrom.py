r"""Nystrom discretisation of the boundary integral operators.

The rectangular rule on equally spaced arc-length nodes is spectrally
accurate for periodic integrands (M&P section 4.6.1), which is why
``Profile.boundary`` resamples at equal arc length. Four operators are
built, and together they form the shared vocabulary of the perfectly
conducting solve and the finite-conductivity coupled system (Goray & Schmidt
2010).

Single layer ``V``: ``\int G(s,s') sigma(s') ds'`` as a dense matrix.
``G`` has an integrable log singularity at ``s' = s``, which is handled
following M&P section 4.6.3 by subtracting the periodic comparison function

.. math::
    \tilde G_\infty(u) = \frac{1}{2\pi}
        \left[\ln\frac{u}{L} + \ln\left(1 - \frac{u}{L}\right)\right],
    \qquad u = |s - s'|,

whose bracket with the unknown is continuous (its diagonal limit is
``\ln(2 \pi L / d) / 2\pi`` plus the remainder diagonal of eq. 4.86) and whose
own integral is exactly ``-L / \pi`` (eq. 4.100). Off the diagonal the matrix
entry is ``spacing * G``, and every singular contribution falls on the
diagonal:

.. math::
    A_{jj} = \Delta s \left[R(0,0) + \frac{\ln(2\pi L/d)}{2\pi}
        - \tilde S\right] - \frac{L}{\pi},
    \qquad
    \tilde S = \frac{1}{\pi}\left[\ln(P-1)! - (P-1)\ln P\right],

with ``\tilde S`` the (node-independent) rectangular-rule sum of
``\tilde G_\infty`` over the off-diagonal nodes.

Double layer ``K``: ``\int N(s,s') psi(s') ds'`` with ``N = dG/dn'``, the
normal being taken at the source point. The sign follows from the
double-layer jump relation for our kernel orientation, and the flat-mirror
and Table 4.1 tests pin it. The kernel is continuous on a smooth curve, so
the plain rectangular rule applies. The diagonal is the numerical
principal-value limit, taken as the average of the kernel evaluated at the
two physically adjacent nodes (through the period wrap), which converges as
O(spacing^2).

Adjoint double layer ``L``: the same kernel machinery with the normal moved
from the source to the observation point and the sign flipped
(``dG/dn = -N``-form with ``n`` taken at ``s``), and the same
principal-value diagonal. It is the principal value of the normal
derivative of the single-layer potential. The jump relation reads
``d(V phi)/dn |_pm = (L pm I/2) phi``, where the ``+`` side is the one into
which the outward normal points. The finite-conductivity system needs this
operator to match normal derivatives across the interface.

Tangential layer ``D_t V``: ``d/ds`` of a single-layer trace, as one
operator. The tangential derivative is continuous across the boundary, so
this is a plain principal-value integral of ``t(s) . grad_s G`` with no jump
term. It is built as a single kernel instead of as a discrete ``d/ds``
applied to an assembled ``V``, because that composition multiplies the
quadrature error of ``V`` by the ``O(N)`` norm of the differentiation
matrix, leaving an ``O(1)`` error that refinement does not remove. The
kernel is genuinely Cauchy-singular (unlike ``K`` and ``L``, whose normal
component cancels the singularity), so the rectangular rule is only
``O(h)`` on it. The singularity is therefore removed with a Bloch-phased
periodic Hilbert comparison kernel whose principal value is exact in the
Fourier basis, leaving a continuous remainder for the usual rule.

Tangential derivative ``D_t``: ``d/ds`` along the boundary of an
``alpha0``-quasi-periodic trace, as a dense differentiation matrix. The
Bloch phase ``exp(i alpha0 x(s))`` is stripped, the now-periodic factor is
differentiated spectrally on the uniform arc-length grid, and the phase is
restored. The product rule contributes ``i alpha0 x'(s)``, and
``x'(s) = n_y`` is the x-component of the tangent for our orientation
(nodes ordered by increasing ``x``, outward normal up). The matrix depends
on geometry only. It is valid only for a trace that is smooth as a function
of arc length, and so it is not suitable for differentiating a single-layer
trace; ``tangential_layer_matrix`` serves that purpose.

The perfectly conducting solve keeps its milestone-1 names as thin
compositions: ``dirichlet_matrix`` (TE) is ``V`` itself and
``neumann_matrix`` (TM) is ``I/2 + K``.

Graded meshes. Everything above is stated on the equal-arc mesh. Each
operator is in fact assembled in the mesh parameter ``tau`` of
``_boundary``, with node weights ``w_l = s'(tau_l) / P`` in place of the
constant spacing, and with each singular comparison function moved from arc
length into ``tau``, where the nodes are always equally spaced:

* ``V``: the comparison is ``H(v) = [ln v + ln(1 - v)] / 2 pi`` with
  ``v = |tau - tau'|``, whose integral is ``-1/pi`` and whose off-diagonal
  rectangular sum is the same ``\tilde S``. Since
  ``G ~ ln(2 pi u / d) / 2 pi + R(0,0)`` and ``u ~ s'_j v``, the diagonal is
  ``w_j [R(0,0) + ln(2 pi s'_j / d) / 2 pi - \tilde S] - s'_j / pi``. With
  ``s' = L`` it is the formula above.
* ``D_t V``: the integrand is the kernel times ``f = sigma s'``, which the
  grading makes smooth even where ``sigma`` is singular. The comparison is
  ``e^{i alpha0 (x - x')} cot(pi (tau - tau')) / (2 s'_j)``, the Hilbert
  part is the same exact Fourier multiplier applied to ``f`` in ``tau``, and
  it is divided by ``s'_j`` afterwards. Applying the multiplier to
  ``sigma`` itself instead is exact on smooth data and wrong in a solve: a
  corner density is unbounded in ``tau``, and it sent the
  finite-conductivity absorption to 1e5 (measured).
* ``D_t``: ``d/ds = (1 / s') d/dtau``, differentiated spectrally in ``tau``.

On a profile with corners this is what lets the graded mesh of
``Profile.boundary(n, grading=p)`` cluster nodes where the Meixner edge
behaviour makes the densities singular, without losing the singular
corrections (M&P section 4.6.6; Kress 1990).

One error the grading does not remove, and :func:`corner_subtraction`
does. For a node a distance ``a`` from a corner, the double-layer kernel
over the facet across the corner varies on the scale ``a``, and the graded
nodes on that facet are spaced on the same scale, so the rectangle rule's
error at that node is a fixed fraction of the integrand that refinement
moves to a closer node instead of shrinking. Measured on a lamellar at
P = 256: up to 0.16 against the identity's 1/2. The perfectly conducting
``I/2 + K`` tolerates it, because the system is second kind and the error
stays on a node of weight ``O(h^p)``. The finite-conductivity system
applies ``K+`` to a trace inside an equation that is first kind in its
density, where it does not stay put, so ``_finite`` asks for the
subtraction. ``L-`` carries the same kind of error and needs none, because
it only ever appears inside ``V+``, whose weights suppress it. ``D_t V-``
also appears only inside ``V+`` and was not suppressed enough, so on this
mesh ``_finite`` replaces it by ``D_t`` applied to the assembled trace.
"""

from __future__ import annotations

from math import lgamma

import numpy as np
from numpy.typing import NDArray

from ...profiles import CORNER_TURN_DEGREES
from ._boundary import PhysicalBoundary
from ._greens import greens_remainder_diagonal, neumann_function
from ._kernels import (
    greens_matrix,
    kernel_geometry,
    neumann_core,
    spectral_factors,
)

__all__ = [
    "OperatorAssembler",
    "corner_subtraction",
    "single_layer_matrix",
    "double_layer_matrix",
    "adjoint_layer_matrix",
    "tangential_layer_matrix",
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
    big_x[diag] = 0.25 * (x.max() - x.min()) + float(boundary.quadrature.max())
    big_z[diag] = 0.0
    return big_x, big_z, diag


def _wrapped_neighbours(
    boundary: PhysicalBoundary, period: float
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """Separations to the two physically adjacent nodes, through the wrap.

    The wrapped neighbours cross the period boundary, so their x separation
    carries the period explicitly and does not pick up a spurious Bloch
    phase from the kernel's quasi-periodicity.
    """
    x, y = boundary.x, boundary.y
    x_next, y_next = np.roll(x, -1), np.roll(y, -1)
    x_next[-1] += period
    x_prev, y_prev = np.roll(x, 1), np.roll(y, 1)
    x_prev[0] -= period
    return x - x_next, y - y_next, x - x_prev, y - y_prev


def _principal_value(
    boundary: PhysicalBoundary,
    forward: NDArray[np.complex128],
    backward: NDArray[np.complex128],
) -> NDArray[np.complex128]:
    """The diagonal limit from the kernel at the two adjacent nodes.

    On the equal-arc mesh this is their average, which converges as
    ``O(spacing^2)`` on a smooth curve. On a graded mesh the node next to a
    corner has one neighbour just across it, where the kernel is of order
    ``1 / distance`` rather than the limit along the node's own straight
    piece, and averaging it in put an ``O(1)`` error on that diagonal entry
    (measured: up to 0.64 against the identity's 1/2). There the limit is
    taken one-sided, from the neighbour on the same piece. A neighbour is
    across a corner when its normal turns by more than
    :data:`~gratinglab.profiles.CORNER_TURN_DEGREES`, the threshold that
    placed the corner.
    """
    average = 0.5 * (forward + backward)
    if not boundary.graded:
        return average
    nx, ny = boundary.nx, boundary.ny
    limit = np.cos(np.radians(CORNER_TURN_DEGREES))
    ahead = nx * np.roll(nx, -1) + ny * np.roll(ny, -1) < limit
    behind = nx * np.roll(nx, 1) + ny * np.roll(ny, 1) < limit
    return np.where(
        ahead & ~behind, backward, np.where(behind & ~ahead, forward, average)
    )


def corner_subtraction(
    boundary: PhysicalBoundary, *, period: float
) -> NDArray[np.float64]:
    r"""The diagonal that singularity-subtracts the double layer at corners.

    Writing ``K phi(s_i) = \int K (phi - phi(s_i)) ds' + phi(s_i) \int K ds'``
    puts the rule's corner error onto an integrand that vanishes at the
    target, and leaves ``\int K ds'`` to be supplied exactly. Its singular
    part is the free-space Laplace double layer
    ``n' . (x' - x) / (2 pi r^2)``, whose integral over a straight segment
    is minus the angle the segment subtends at ``x``, over ``2 pi``. The
    returned ``c_i`` is that exact integral minus the rectangle rule's
    value, over the polygon and the images in the two neighbouring periods,
    so adding ``diag(c)`` to the discrete ``K`` is the subtraction. Farther
    images and the remainder ``K - K_Laplace``, which is continuous, are
    integrated accurately by the rule and need nothing.

    No Bloch phase enters: across the period wrap the kernel's phase and the
    quasi-periodic density's phase cancel in the product the correction
    stands in for. The target's own segment contributes zero (its principal
    value on a straight piece), and it is found as the segment nearest the
    node rather than by testing collinearity, which rounding in the
    coordinates of a node 1e-9 nm from a vertex defeats (measured: a
    misread segment adds a whole 1/2). Zero everywhere on the equal-arc
    mesh, which is left exactly as it was.
    """
    points = len(boundary.x)
    if boundary.polygon is None:
        return np.zeros(points)
    vx, vy = boundary.polygon
    big_x = np.concatenate((vx - period, vx, vx + period, [vx[0] + 2.0 * period]))
    big_y = np.concatenate((vy, vy, vy, [vy[0]]))
    ax, ay, bx, by = big_x[:-1], big_y[:-1], big_x[1:], big_y[1:]
    x, y = boundary.x[:, None], boundary.y[:, None]

    to_a_x, to_a_y, to_b_x, to_b_y = ax - x, ay - y, bx - x, by - y
    angle = np.arctan2(
        to_a_x * to_b_y - to_a_y * to_b_x, to_a_x * to_b_x + to_a_y * to_b_y
    )
    ex, ey = bx - ax, by - ay
    along = np.clip(((x - ax) * ex + (y - ay) * ey) / (ex**2 + ey**2), 0.0, 1.0)
    own = np.argmin(np.hypot(ax + along * ex - x, ay + along * ey - y), axis=1)
    angle[np.arange(points), own] = 0.0
    exact = -angle.sum(axis=1) / (2.0 * np.pi)

    rule = np.zeros(points)
    weights = boundary.quadrature
    for image in (-1, 0, 1):
        dx = boundary.x[None, :] + image * period - boundary.x[:, None]
        dy = boundary.y[None, :] - boundary.y[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            kernel = (boundary.nx[None, :] * dx + boundary.ny[None, :] * dy) / (
                2.0 * np.pi * (dx**2 + dy**2)
            )
        if image == 0:
            # A straight piece's own principal value, as in ``exact``.
            np.fill_diagonal(kernel, 0.0)
        rule += kernel @ weights
    return exact - rule


class OperatorAssembler:
    """One boundary's kernel machinery, shared across a whole scan.

    The assembler owns the wavelength-independent half of the separable
    kernel assembly (``_kernels.KernelGeometry``, built once) and memoizes
    the per-(k, alpha0) spectral factors and Neumann cores. The operators
    that share a medium at one wavelength (``V``/``K`` on the vacuum side,
    ``V``/ ``L``/``D_t V`` on the metal side) therefore pay for each
    spectral sum once. The solver clears the memo at every wavelength to
    keep memory flat, while the geometry (and its cached Kummer tails)
    persists for the scan.
    """

    def __init__(
        self, boundary: PhysicalBoundary, *, period: float, terms: int
    ) -> None:
        self.boundary = boundary
        self.period = period
        self.terms = terms
        self.geometry = kernel_geometry(boundary, period=period, terms=terms)
        self._factors: dict = {}
        self._cores: dict = {}
        self._hilbert: NDArray[np.complex128] | None = None
        self._subtraction: NDArray[np.float64] | None = None

    @property
    def sgn(self) -> NDArray[np.float64]:
        """``sign(y_j - y_l)``, the Neumann compositions' third factor."""
        return self.geometry.sgn

    def factors(self, *, k: complex, alpha0: float):
        key = (complex(k), float(alpha0))
        if key not in self._factors:
            self._factors[key] = spectral_factors(
                self.geometry, k=k, alpha0=alpha0
            )
        return self._factors[key]

    def greens(self, *, k: complex, alpha0: float) -> NDArray[np.complex128]:
        """The full Green's-function matrix, with a meaningless diagonal
        (overwritten by the caller's analytic limit)."""
        return greens_matrix(self.geometry, self.factors(k=k, alpha0=alpha0))

    def neumann(
        self, *, k: complex, alpha0: float
    ) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
        """The Neumann kernel's ``(A, B)`` components, memoized. The double,
        adjoint, and tangential layers at one (k, alpha0) are all composed
        from this single evaluation."""
        key = (complex(k), float(alpha0))
        if key not in self._cores:
            self._cores[key] = neumann_core(
                self.geometry, self.factors(k=k, alpha0=alpha0)
            )
        return self._cores[key]

    def hilbert(self) -> NDArray[np.complex128]:
        """The periodic Hilbert-multiplier matrix, which depends on geometry
        only. Fourier mode ``p`` is multiplied by ``-i sgn(p) / 2``, and the
        Nyquist mode is zeroed."""
        if self._hilbert is None:
            points = len(self.boundary.x)
            modes = np.fft.fftfreq(points, d=1.0 / points)
            multiplier = -0.5j * np.sign(modes)
            if points % 2 == 0:
                # The unpaired Nyquist mode has no signed representation.
                multiplier[points // 2] = 0.0
            self._hilbert = np.fft.ifft(
                multiplier[:, None] * np.fft.fft(np.eye(points), axis=0),
                axis=0,
            )
        return self._hilbert

    def subtraction(self) -> NDArray[np.float64]:
        """:func:`corner_subtraction` for this boundary, which depends on
        geometry only and so is built once per scan."""
        if self._subtraction is None:
            self._subtraction = corner_subtraction(
                self.boundary, period=self.period
            )
        return self._subtraction

    def clear(self) -> None:
        """Drop the per-wavelength memo; the geometry stays."""
        self._factors.clear()
        self._cores.clear()


def single_layer_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
    assembler: OperatorAssembler | None = None,
) -> NDArray[np.complex128]:
    """The discretised single-layer operator ``V``."""
    weights, speed = boundary.quadrature, boundary.speed
    points = len(boundary.x)
    diag = np.eye(points, dtype=bool)

    if assembler is None:
        assembler = OperatorAssembler(boundary, period=period, terms=terms)
    matrix = assembler.greens(k=k, alpha0=alpha0) * weights[None, :]

    remainder = greens_remainder_diagonal(
        k=k, alpha0=alpha0, period=period, terms=terms
    )
    comparison_sum = (lgamma(points) - (points - 1) * np.log(points)) / np.pi
    matrix[diag] = (
        weights
        * (
            remainder
            + np.log(2.0 * np.pi * speed / period) / (2.0 * np.pi)
            - comparison_sum
        )
        - speed / np.pi
    )
    return matrix


def double_layer_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
    assembler: OperatorAssembler | None = None,
    subtract_corners: bool = False,
) -> NDArray[np.complex128]:
    """The discretised double-layer operator ``K`` (no identity part).

    ``subtract_corners`` adds :func:`corner_subtraction` to the diagonal,
    which is exact for a constant density and changes nothing on the
    equal-arc mesh. The perfectly conducting solve leaves it off.
    """
    weights = boundary.quadrature
    nx, ny = boundary.nx, boundary.ny
    diag = np.eye(len(boundary.x), dtype=bool)

    if assembler is None:
        assembler = OperatorAssembler(boundary, period=period, terms=terms)
    a_matrix, b_matrix = assembler.neumann(k=k, alpha0=alpha0)
    kernel = a_matrix * nx[None, :] + assembler.sgn * b_matrix * ny[None, :]

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
    kernel[diag] = _principal_value(boundary, forward, backward)

    matrix = kernel * weights[None, :]
    if subtract_corners and boundary.graded:
        matrix[diag] += assembler.subtraction()
    return matrix


def adjoint_layer_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
    assembler: OperatorAssembler | None = None,
) -> NDArray[np.complex128]:
    """The discretised adjoint double-layer operator ``L``.

    The kernel is the normal derivative of ``G`` at the observation point:
    ``n(s) . grad G = -``\\ [the ``N``-form with the normal taken at ``s``],
    so the same accelerated kernel serves with the normals broadcast along
    rows instead of columns and the sign flipped. The diagonal is the same
    adjacent-node principal-value average, with the node's own normal.
    """
    weights = boundary.quadrature
    nx, ny = boundary.nx, boundary.ny
    diag = np.eye(len(boundary.x), dtype=bool)

    if assembler is None:
        assembler = OperatorAssembler(boundary, period=period, terms=terms)
    a_matrix, b_matrix = assembler.neumann(k=k, alpha0=alpha0)
    kernel = -(a_matrix * nx[:, None] + assembler.sgn * b_matrix * ny[:, None])

    dx_next, dy_next, dx_prev, dy_prev = _wrapped_neighbours(boundary, period)
    forward = -neumann_function(
        dx_next, dy_next, nx, ny, k=k, alpha0=alpha0, period=period, terms=terms
    )
    backward = -neumann_function(
        dx_prev, dy_prev, nx, ny, k=k, alpha0=alpha0, period=period, terms=terms
    )
    kernel[diag] = _principal_value(boundary, forward, backward)

    return kernel * weights[None, :]


def tangential_layer_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
    assembler: OperatorAssembler | None = None,
) -> NDArray[np.complex128]:
    r"""The discretised operator ``D_t V``, which is ``d/ds`` of a single-layer trace.

    The tangential derivative of a single-layer potential is continuous
    across the boundary (only the normal derivative jumps), so the
    arc-length derivative of the trace is the plain principal-value integral

    .. math::
        \partial_s (V \sigma)(s)
            = \int t(s) \cdot \nabla_s G(s, s')\, \sigma(s')\, ds',

    with no jump term, where ``t = (n_y, -n_x)`` is the unit tangent in the
    ``+s`` direction and ``t . grad_s = -``\ [the ``N``-form with the tangent
    at ``s``], which is the same sign flip that :func:`adjoint_layer_matrix`
    makes.

    The operator is usable because the derivative is taken analytically on
    the kernel. The alternative, assembling ``V`` and then applying a
    discrete ``d/ds``, composes an unbounded operator (discrete norm
    ``O(N)``) with a quadrature carrying ``O(N^{-p})`` error. This leaves an
    ``O(1)`` error that refinement does not remove, and on a boundary with
    any appreciable high-frequency content it destroys the solve entirely.

    Unlike ``K`` and ``L``, whose normal component cancels the singularity,
    this kernel is genuinely Cauchy-singular: ``t . grad_s G ~ 1/(2 pi u)``
    for ``u = s - s'``. The rectangular rule is only ``O(h)`` on it, so the
    singularity is removed in the way ``single_layer_matrix`` removes its
    logarithm, by subtracting a periodic comparison kernel whose action is
    known in closed form. Here that is the Bloch-phased periodic Hilbert
    kernel

    .. math::
        C(s, s') = e^{i \alpha_0 (x(s) - x(s'))}\,
                   \frac{1}{2L} \cot\frac{\pi (s - s')}{L},

    which carries the same quasi-periodicity and the same ``1/(2 pi u)``
    behaviour. Its principal-value integral is diagonal in the Fourier basis
    of the periodic factor (mode ``p`` is multiplied by ``-i sgn(p) / 2``),
    so it is applied spectrally and exactly. The remainder ``K - C`` is
    continuous, so the plain rectangular rule recovers its usual accuracy,
    with the diagonal taken as the adjacent-node average as elsewhere.
    """
    speed = boundary.speed
    points = len(boundary.x)
    diag = np.eye(points, dtype=bool)
    tx, ty = boundary.ny, -boundary.nx

    def kernel_at(dx, dy, vx, vy):
        return -neumann_function(
            dx, dy, vx, vy, k=k, alpha0=alpha0, period=period, terms=terms
        )

    if assembler is None:
        assembler = OperatorAssembler(boundary, period=period, terms=terms)
    a_matrix, b_matrix = assembler.neumann(k=k, alpha0=alpha0)
    kernel = -(a_matrix * tx[:, None] + assembler.sgn * b_matrix * ty[:, None])

    # Everything below is in the mesh parameter tau, whose nodes are equally
    # spaced on any mesh. The integrand is the kernel times f = sigma s',
    # smooth in tau where sigma is not, and the comparison is
    # cot(pi (tau_j - tau_i)) / (2 s'_j), which matches the kernel as
    # tau' -> tau_j because u ~ s'_j (tau_j - tau'). On the equal-arc mesh
    # (s' = L) this is the arc-length cot(pi u / L) / 2L.
    phase = np.exp(1j * alpha0 * boundary.x)
    offset = (np.arange(points)[:, None] - np.arange(points)[None, :]) / points
    with np.errstate(divide="ignore", invalid="ignore"):
        comparison = (
            (phase[:, None] / phase[None, :]) / 2.0 / np.tan(np.pi * offset)
        )
    comparison[diag] = 0.0

    remainder = (kernel - comparison / speed[:, None]) * speed[None, :]

    # The remainder is continuous; its diagonal is the adjacent-node average,
    # taken on the same difference so the comparison's own 1/u cancels there
    # too (its adjacent values are equal and opposite, hence absent).
    dx_next, dy_next, dx_prev, dy_prev = _wrapped_neighbours(boundary, period)
    step = 1.0 / (2.0 * np.tan(np.pi / points))
    comparison_forward = -step * np.exp(1j * alpha0 * dx_next) / speed
    comparison_backward = step * np.exp(1j * alpha0 * dx_prev) / speed
    forward = (kernel_at(dx_next, dy_next, tx, ty) - comparison_forward) * np.roll(
        speed, -1
    )
    backward = (kernel_at(dx_prev, dy_prev, tx, ty) - comparison_backward) * np.roll(
        speed, 1
    )
    remainder[diag] = _principal_value(boundary, forward, backward)

    # The comparison kernel's principal value, applied exactly to f: strip
    # the Bloch phase, multiply Fourier mode p by -i sgn(p) / 2, restore it,
    # then divide by s'_j.
    hilbert = assembler.hilbert()

    return remainder / points + (phase[:, None] / phase[None, :]) * hilbert * (
        speed[None, :] / speed[:, None]
    )


def tangential_derivative_matrix(
    boundary: PhysicalBoundary, *, alpha0: float
) -> NDArray[np.complex128]:
    """The discretised tangential derivative ``D_t`` for quasi-periodic traces.

    Applied to node values of an ``alpha0``-quasi-periodic function of arc
    length (such as the trace of a single-layer potential), the matrix
    returns ``d/ds``. The Bloch phase is stripped, the periodic factor is
    differentiated spectrally on the uniform arc-length grid, the phase is
    restored, and the product-rule term ``i alpha0 x'(s)`` is added with
    ``x'(s) = n_y``. The matrix depends on geometry only, and no wavelength
    enters.
    """
    points = len(boundary.x)
    # d/dtau spectrally on the equally spaced tau grid; d/ds = d/dtau / s'.
    wavenumber = 2j * np.pi * np.fft.fftfreq(points, d=1.0 / points)
    if points % 2 == 0:
        # The unpaired Nyquist mode has no odd (derivative) representation.
        wavenumber[points // 2] = 0.0
    spectral = np.fft.ifft(
        wavenumber[:, None] * np.fft.fft(np.eye(points), axis=0), axis=0
    ).real / boundary.speed[:, None]
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
    assembler: OperatorAssembler | None = None,
) -> NDArray[np.complex128]:
    """The TE (Dirichlet) system matrix, which is the single-layer operator itself."""
    return single_layer_matrix(
        boundary, k=k, alpha0=alpha0, period=period, terms=terms,
        assembler=assembler,
    )


def neumann_matrix(
    boundary: PhysicalBoundary,
    *,
    k: complex,
    alpha0: float,
    period: float,
    terms: int,
    assembler: OperatorAssembler | None = None,
) -> NDArray[np.complex128]:
    """``I/2 + K`` for the TM (Neumann) problem."""
    return 0.5 * np.eye(len(boundary.x)) + double_layer_matrix(
        boundary, k=k, alpha0=alpha0, period=period, terms=terms,
        assembler=assembler,
    )
