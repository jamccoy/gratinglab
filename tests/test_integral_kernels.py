"""The separable (GEMM) kernel assembly against the pointwise reference.

``_kernels`` rebuilds the sums of ``_greens`` as BLAS products with
geometry-cached Kummer tails; nothing about the physics may move. Every
case compares the assembled ``(P, P)`` matrices off-diagonal against the
pointwise kernels on the same separations, to near machine precision,
across the sign patterns (monotone, oscillating, all-flat) and wavenumbers
(real, X-ray-like, visible-metal-like) the solver actually meets -- plus
both overflow regimes of the split factors: the guard forced low so the
fallback carries almost everything, and a genuinely deep groove where it
trips on its own.
"""

import numpy as np
import pytest

from gratinglab.problem import Problem
from gratinglab.profiles import Sinusoidal
from gratinglab.solvers.integral._boundary import PhysicalBoundary, physical_boundary
from gratinglab.solvers.integral import _kernels
from gratinglab.solvers.integral._greens import greens_function, neumann_function
from gratinglab.solvers.integral._kernels import (
    greens_matrix,
    kernel_geometry,
    neumann_core,
    spectral_factors,
)
from gratinglab.solvers.integral._nystrom import _separations

PERIOD = 600.0
K_T = 2.0 * np.pi / 1200.0
ALPHA0 = -K_T * 0.3

TOLERANCE = 1e-12


def sinusoid(points=96, depth=0.3):
    return physical_boundary(
        Problem(period=PERIOD, profile=Sinusoidal(depth_fraction=depth)), points
    )


def flat(points=96):
    """Every off-diagonal ``Z`` is exactly zero -- the maximal tie case."""
    x = np.linspace(0.0, PERIOD, points, endpoint=False)
    return PhysicalBoundary(
        x=x,
        y=np.zeros(points),
        nx=np.zeros(points),
        ny=np.ones(points),
        arc_length=PERIOD,
    )


def two_level(points=96):
    """A smoothed lamellar-like curve: long equal-height plateaus give
    off-diagonal ``Z = 0`` on both sides of the diagonal, and the steps give
    an arbitrary (non-monotone) sign pattern."""
    t = np.linspace(0.0, 1.0, points, endpoint=False)
    y = 120.0 * (np.tanh(12.0 * np.sin(2.0 * np.pi * t)) + 1.0) / 2.0
    y = np.round(y, 6)  # flatten the plateaus to exact ties
    x = t * PERIOD
    angle = np.arctan2(np.gradient(y), np.gradient(x))
    return PhysicalBoundary(
        x=x,
        y=y,
        nx=-np.sin(angle),
        ny=np.cos(angle),
        arc_length=float(np.sum(np.hypot(np.gradient(x), np.gradient(y)))),
    )


BOUNDARIES = {
    "sinusoid": sinusoid(),
    "deep-sinusoid": sinusoid(depth=1.1),
    "flat": flat(),
    "two-level": two_level(),
}

WAVENUMBERS = {
    "real": K_T,
    "xray": K_T * (0.9 + 0.2j),
    "gold": K_T * (0.2 + 3.0j),
}


def reference_and_new(boundary, k, terms=48):
    """Both paths on the same separations, plus the off-diagonal mask."""
    big_x, big_z, diag = _separations(boundary)
    shared = dict(k=k, alpha0=ALPHA0, period=PERIOD, terms=terms)
    geometry = kernel_geometry(boundary, period=PERIOD, terms=terms)
    factors = spectral_factors(geometry, k=k, alpha0=ALPHA0)
    return big_x, big_z, ~diag, shared, geometry, factors


def assert_close(new, reference, off, tolerance=TOLERANCE):
    # The flat boundary with a vertical normal nulls the Neumann kernel
    # exactly (sgn = 0 everywhere off-diagonal); 1.0 keeps the comparison
    # meaningful as an absolute one there.
    scale = max(np.abs(reference[off]).max(), 1.0e-300)
    worst = np.abs(new[off] - reference[off]).max()
    assert worst < tolerance * max(scale, tolerance), (
        f"max deviation {worst / scale:.3e} relative"
    )


@pytest.mark.parametrize("name", BOUNDARIES)
@pytest.mark.parametrize("label", WAVENUMBERS)
class TestAgainstPointwiseKernels:
    def test_greens(self, name, label):
        boundary, k = BOUNDARIES[name], WAVENUMBERS[label]
        big_x, big_z, off, shared, geometry, factors = reference_and_new(
            boundary, k
        )
        reference = greens_function(big_x, big_z, **shared)
        assert_close(greens_matrix(geometry, factors), reference, off)

    def test_neumann_source_side(self, name, label):
        """The double-layer composition: the vector along columns."""
        boundary, k = BOUNDARIES[name], WAVENUMBERS[label]
        big_x, big_z, off, shared, geometry, factors = reference_and_new(
            boundary, k
        )
        nx, ny = boundary.nx, boundary.ny
        reference = neumann_function(
            big_x, big_z, nx[None, :], ny[None, :], **shared
        )
        a_matrix, b_matrix = neumann_core(geometry, factors)
        new = a_matrix * nx[None, :] + geometry.sgn * b_matrix * ny[None, :]
        assert_close(new, reference, off)

    def test_neumann_observation_side(self, name, label):
        """The adjoint/tangential composition: the vector along rows."""
        boundary, k = BOUNDARIES[name], WAVENUMBERS[label]
        big_x, big_z, off, shared, geometry, factors = reference_and_new(
            boundary, k
        )
        nx, ny = boundary.nx, boundary.ny
        reference = neumann_function(
            big_x, big_z, nx[:, None], ny[:, None], **shared
        )
        a_matrix, b_matrix = neumann_core(geometry, factors)
        new = a_matrix * nx[:, None] + geometry.sgn * b_matrix * ny[:, None]
        assert_close(new, reference, off)


class TestOverflowGuard:
    def test_forced_fallback_matches(self, monkeypatch):
        """With the guard forced to ~0 the broadcast fallback carries every
        order with any evanescence; the assembly must not care which path
        an order took."""
        boundary, k = BOUNDARIES["sinusoid"], WAVENUMBERS["gold"]
        big_x, big_z, off, shared, geometry, _ = reference_and_new(boundary, k)
        monkeypatch.setattr(_kernels, "_EXP_GUARD", 1.0)
        factors = spectral_factors(geometry, k=k, alpha0=ALPHA0)
        assert factors.alpha_hard.size > 0
        reference = greens_function(big_x, big_z, **shared)
        assert_close(greens_matrix(geometry, factors), reference, off)

    def test_deep_groove_trips_the_guard_naturally(self):
        """``Im(gamma_M) h / 2`` past the guard is where a naive split would
        produce ``0 * inf``; the hybrid must sail through and still match
        the pointwise reference, which never splits."""
        boundary = sinusoid(points=64, depth=1.4)
        terms = 400  # Im(gamma_M) h / 2 ~ pi * terms * h / period > 690
        k = WAVENUMBERS["real"]
        big_x, big_z, diag = _separations(boundary)
        geometry = kernel_geometry(boundary, period=PERIOD, terms=terms)
        factors = spectral_factors(geometry, k=k, alpha0=ALPHA0)
        assert factors.alpha_hard.size > 0, "case does not exercise the guard"
        reference = greens_function(
            big_x, big_z, k=k, alpha0=ALPHA0, period=PERIOD, terms=terms
        )
        new = greens_matrix(geometry, factors)
        assert np.isfinite(new[~diag]).all()
        assert_close(new, reference, ~diag)
