"""The operator vocabulary of the finite-conductivity system.

The perfectly conducting solve needs only ``V`` (TE) and ``I/2 + K`` (TM);
the coupled system of Goray & Schmidt (2010) adds the adjoint double layer
``L`` and the tangential derivative ``D_t``. The Stage-0 contract: the new
operators are pinned by kernel-level identities -- the single-layer jump
relation, a Bloch-reversal transpose identity, exactness on trig
polynomials -- before any physics is built on them, and the milestone-1
names remain thin compositions of the shared pieces.
"""

import numpy as np
import pytest

from gratinglab.problem import Problem
from gratinglab.profiles import Sinusoidal
from gratinglab.solvers.integral._boundary import PhysicalBoundary, physical_boundary
from gratinglab.solvers.integral._greens import greens_function
from gratinglab.solvers.integral._nystrom import (
    adjoint_layer_matrix,
    dirichlet_matrix,
    double_layer_matrix,
    neumann_matrix,
    single_layer_matrix,
    tangential_derivative_matrix,
    tangential_layer_matrix,
)

PERIOD = 600.0
WAVELENGTH = 1200.0
K_T = 2.0 * np.pi / WAVELENGTH
ALPHA0 = -K_T * 0.3

SINUSOID = physical_boundary(
    Problem(period=PERIOD, profile=Sinusoidal(depth_fraction=0.3)), 384
)


def smooth_density(boundary, alpha0):
    """A smooth quasi-periodic density: Bloch phase times low harmonics of
    the arc-length parameter (the nodes are the uniform grid in s)."""
    j = np.arange(len(boundary.x))
    turns = 2.0 * np.pi * j / len(j)
    periodic = 1.0 + 0.3 * np.cos(turns) + 0.2 * np.sin(2.0 * turns)
    return np.exp(1j * alpha0 * boundary.x) * periodic


class TestMilestoneOneNames:
    """The TE/TM system matrices are compositions of the shared operators,
    bit for bit -- the guarantee that the finite-conductivity refactor
    cannot have moved the perfectly conducting answers."""

    def test_dirichlet_matrix_is_the_single_layer(self):
        kwargs = dict(k=K_T, alpha0=ALPHA0, period=PERIOD, terms=60)
        assert np.array_equal(
            dirichlet_matrix(SINUSOID, **kwargs),
            single_layer_matrix(SINUSOID, **kwargs),
        )

    def test_neumann_matrix_is_half_identity_plus_double_layer(self):
        kwargs = dict(k=K_T, alpha0=ALPHA0, period=PERIOD, terms=60)
        expected = 0.5 * np.eye(len(SINUSOID.x)) + double_layer_matrix(
            SINUSOID, **kwargs
        )
        assert np.array_equal(neumann_matrix(SINUSOID, **kwargs), expected)


class TestAdjointLayer:
    def test_bloch_reversal_transpose_identity(self):
        """``L(alpha0) = K(-alpha0)^T`` off the diagonal: reversing the
        Bloch phase and transposing swaps observation and source point in
        the kernel exactly (G(-X,-Z; -alpha0) = G(X,Z; alpha0)). This pins
        the sign flip and the row-wise normals against the already-tested
        double layer."""
        kwargs = dict(k=K_T, period=PERIOD, terms=60)
        left = adjoint_layer_matrix(SINUSOID, alpha0=ALPHA0, **kwargs)
        right = double_layer_matrix(SINUSOID, alpha0=-ALPHA0, **kwargs).T
        off = ~np.eye(len(SINUSOID.x), dtype=bool)
        assert np.allclose(left[off], right[off], rtol=0.0, atol=1e-10)

    @pytest.mark.parametrize("side", [+1.0, -1.0], ids=["vacuum", "metal"])
    def test_single_layer_jump_relation(self, side):
        """``d(V phi)/dn`` from the side the outward normal points into is
        ``(L + I/2) phi``; from the other side ``(L - I/2) phi``. Probed by
        a one-sided third-order finite difference of the potential along
        the normal, anchored on the (continuous) trace ``V phi``. This is
        the test that pins L's principal-value diagonal and the half-jump
        the finite-conductivity system matches across the interface."""
        boundary = SINUSOID
        terms = 192
        kwargs = dict(k=K_T, alpha0=ALPHA0, period=PERIOD, terms=terms)
        phi = smooth_density(boundary, ALPHA0)

        trace = single_layer_matrix(boundary, **kwargs) @ phi
        expected = (
            adjoint_layer_matrix(boundary, **kwargs) @ phi
            + side * 0.5 * phi
        )

        probes = np.arange(0, len(boundary.x), 48)
        h = 3.0 * boundary.spacing

        def potential(step):
            px = boundary.x[probes] + side * step * boundary.nx[probes]
            py = boundary.y[probes] + side * step * boundary.ny[probes]
            big_x = px[:, None] - boundary.x[None, :]
            big_z = py[:, None] - boundary.y[None, :]
            kernel = greens_function(big_x, big_z, **kwargs)
            return boundary.spacing * (kernel @ phi)

        u0 = trace[probes]
        u1, u2, u3 = potential(h), potential(2.0 * h), potential(3.0 * h)
        derivative = side * (
            -11.0 * u0 + 18.0 * u1 - 9.0 * u2 + 2.0 * u3
        ) / (6.0 * h)

        assert np.allclose(derivative, expected[probes], rtol=0.0, atol=1e-3)


class TestTangentialDerivative:
    def test_exact_on_a_flat_boundary(self):
        """On a flat boundary arc length is x, so a quasi-periodic trig
        polynomial has a closed-form tangential derivative and the spectral
        matrix must reproduce it to machine precision."""
        points = 64
        x = np.linspace(0.0, PERIOD, points, endpoint=False)
        boundary = PhysicalBoundary(
            x=x,
            y=np.zeros(points),
            nx=np.zeros(points),
            ny=np.ones(points),
            arc_length=PERIOD,
        )
        omega = 2.0 * np.pi * 3.0 / PERIOD
        q = np.exp(1j * omega * x) + 0.5 * np.exp(-2j * omega * x)
        dq = 1j * omega * np.exp(1j * omega * x) - 1j * omega * np.exp(
            -2j * omega * x
        )
        f = np.exp(1j * ALPHA0 * x) * q
        expected = np.exp(1j * ALPHA0 * x) * (dq + 1j * ALPHA0 * q)
        got = tangential_derivative_matrix(boundary, alpha0=ALPHA0) @ f
        assert np.allclose(got, expected, rtol=0.0, atol=1e-10)

    def test_matches_finite_differences_on_a_sinusoid(self):
        """On a curved boundary the product rule's ``i alpha0 x'(s)`` term
        carries the tangent convention ``x'(s) = n_y``; a fourth-order
        finite difference on the arc-length grid (wrapped with the Bloch
        factor) sees the true ``x'(s)`` through the samples and would catch
        a wrong convention at ~1e-4 against a ~1e-6 discretisation floor."""
        boundary = SINUSOID
        f = smooth_density(boundary, ALPHA0)
        bloch = np.exp(1j * ALPHA0 * PERIOD)

        def shifted(offset):
            rolled = np.roll(f, -offset)
            if offset > 0:
                rolled[-offset:] *= bloch
            elif offset < 0:
                rolled[:-offset] /= bloch
            return rolled

        spacing = boundary.spacing
        fd = (
            8.0 * (shifted(1) - shifted(-1)) - (shifted(2) - shifted(-2))
        ) / (12.0 * spacing)

        got = tangential_derivative_matrix(boundary, alpha0=ALPHA0) @ f
        scale = np.abs(fd).max()
        assert np.allclose(got, fd, rtol=0.0, atol=1e-4 * scale)


class TestTangentialLayer:
    """``D_t V`` is one operator, not a discrete ``d/ds`` applied to ``V``.

    Composing an assembled ``V`` with a differentiation matrix multiplies
    ``V``'s quadrature error by that matrix's ``O(N)`` norm, so the product
    carries an ``O(1)`` error that refinement never removes. These pin the
    single-operator form instead: exactness against a closed form where one
    exists, and convergence where it does not.
    """

    @staticmethod
    def _flat(points, period):
        """A flat interface, where arc length is ``x`` and the tangent is
        ``+x``, so ``V`` acts diagonally on every quasi-periodic mode."""
        x = np.linspace(0.0, period, points, endpoint=False)
        return PhysicalBoundary(
            x=x,
            y=np.zeros(points),
            nx=np.zeros(points),
            ny=np.ones(points),
            arc_length=period,
        )

    @pytest.mark.parametrize("k", [0.15, 0.11163 + 0.20861j])
    @pytest.mark.parametrize("order", [0, 1, 5])
    def test_flat_interface_matches_the_closed_form(self, k, order):
        """On a flat interface ``V`` maps ``exp(i alpha_m x)`` to a multiple
        of itself, so ``d/ds (V f) = i alpha_m (V f)`` exactly -- including
        for the complex metal-side wavenumber, where the kernel is strongly
        evanescent."""
        period, alpha0 = PERIOD, -0.04
        errors = []
        for points in (128, 256):
            boundary = self._flat(points, period)
            shared = dict(
                k=k, alpha0=alpha0, period=period, terms=points // 2
            )
            single = single_layer_matrix(boundary, **shared)
            tangential = tangential_layer_matrix(boundary, **shared)

            alpha_m = alpha0 + 2.0 * np.pi * order / period
            mode = np.exp(1j * alpha_m * boundary.x)
            reference = single @ mode
            errors.append(
                np.abs(tangential @ mode - 1j * alpha_m * reference).max()
                / np.abs(reference).max()
            )

        assert errors[0] < 1e-2
        # Halving the spacing buys ~8x (the remainder quadrature is cubic).
        # The plain rectangular rule on this Cauchy kernel manages only a
        # factor of two, which is the failure this guards; 4 separates them.
        assert errors[1] < errors[0] / 4.0

    def test_agrees_with_spectral_differentiation_on_a_smooth_boundary(self):
        """Where the trace *is* smooth in arc length, the two routes must
        agree -- which is what makes the discrete one so plausible, and why
        only a high-frequency boundary exposes the difference."""
        shared = dict(k=K_T, alpha0=ALPHA0, period=PERIOD, terms=192)
        single = single_layer_matrix(SINUSOID, **shared)
        tangential = tangential_layer_matrix(SINUSOID, **shared)
        discrete = tangential_derivative_matrix(SINUSOID, alpha0=ALPHA0)

        density = smooth_density(SINUSOID, ALPHA0)
        analytic = tangential @ density
        spectral = discrete @ (single @ density)

        assert np.abs(analytic - spectral).max() < 1e-3 * np.abs(analytic).max()
