"""The quasi-periodic Green's function under the integral method.

The accelerated forms are validated against the raw spectral sum, which is
slow but unambiguous. Where the raw sum converges too slowly to be a
reference (on the boundary, ``Z = 0``), the test is self-consistency: the
acceleration must make the answer independent of the truncation.
"""

import numpy as np
import pytest

from gratinglab.solvers.integral._greens import (
    greens_function,
    greens_remainder_diagonal,
    neumann_function,
    wavenumbers,
)

PERIOD = 600.0
WAVELENGTH = 550.0
K_T = 2.0 * np.pi / WAVELENGTH
ALPHA0 = -K_T * np.sin(np.radians(30.0))


def brute_force_greens(x, z, *, k, alpha0, period, terms):
    """The raw spectral sum, no acceleration. Reference only."""
    m = np.arange(-terms, terms + 1)
    alpha_m, gamma_m = wavenumbers(k, alpha0, period, m)
    phase = np.exp(
        1j * alpha_m * np.asarray(x)[..., None]
        + 1j * gamma_m * np.abs(z)[..., None]
    )
    return (phase / (2j * period * gamma_m)).sum(axis=-1)


def brute_force_neumann(x, z, nx, ny, *, k, alpha0, period, terms):
    m = np.arange(-terms, terms + 1)
    alpha_m, gamma_m = wavenumbers(k, alpha0, period, m)
    phase = np.exp(
        1j * alpha_m * np.asarray(x)[..., None]
        + 1j * gamma_m * np.abs(z)[..., None]
    )
    bracket = (alpha_m / gamma_m) * np.asarray(nx)[..., None] + (
        np.sign(z) * np.asarray(ny)
    )[..., None]
    return -(phase * bracket).sum(axis=-1) / (2.0 * period)


class TestGreensFunction:
    def test_accelerated_matches_brute_force_off_the_boundary(self):
        """With K|Z| ~ 0.5 the raw sum converges exponentially, so it is a
        genuine reference; the accelerated form must agree to near machine
        precision with far fewer terms."""
        x = np.array([13.0, 150.0, -290.0, 411.0])
        z = np.array([55.0, 47.0, 61.0, -80.0])
        slow = brute_force_greens(
            x, z, k=K_T, alpha0=ALPHA0, period=PERIOD, terms=400
        )
        fast = greens_function(
            x, z, k=K_T, alpha0=ALPHA0, period=PERIOD, terms=40
        )
        assert np.allclose(fast, slow, rtol=0.0, atol=1e-12)

    def test_truncation_independence_on_the_boundary(self):
        """At Z = 0 the raw sum is useless (1/m, oscillatory), which is the
        whole reason the acceleration exists. The paired remainder decays like
        m^-2 here (m^-3 only on the diagonal, M&P eq. 4.87), so the measured
        contract is: truncation error at M = 500 stays below 1e-5 and keeps
        shrinking."""
        x = np.array([9.0, 30.0, 111.0, 517.0])
        z = np.zeros_like(x)
        coarse = greens_function(
            x, z, k=K_T, alpha0=ALPHA0, period=PERIOD, terms=500
        )
        fine = greens_function(
            x, z, k=K_T, alpha0=ALPHA0, period=PERIOD, terms=8000
        )
        assert np.abs(coarse - fine).max() < 1e-5

    def test_quasi_periodicity(self):
        """G(X + d, Z) = exp(i alpha0 d) G(X, Z) -- exercises the closed-form
        logs and prefactors, which must carry the Bloch phase exactly."""
        x = np.array([25.0, 260.0, -140.0])
        z = np.array([12.0, 0.0, 33.0])
        base = greens_function(x, z, k=K_T, alpha0=ALPHA0, period=PERIOD, terms=80)
        shifted = greens_function(
            x + PERIOD, z, k=K_T, alpha0=ALPHA0, period=PERIOD, terms=80
        )
        assert np.allclose(shifted, base * np.exp(1j * ALPHA0 * PERIOD), atol=1e-9)

    def test_diagonal_limit_matches_the_remainder_series(self):
        """Approaching the singular point along the boundary, G minus its
        singular asymptote converges to the closed diagonal value
        (M&P eq. 4.86)."""
        expected = greens_remainder_diagonal(
            k=K_T, alpha0=ALPHA0, period=PERIOD, terms=4000
        )
        eps = 1e-4  # nm; well below any feature scale
        x = np.array([eps])
        z = np.array([0.0])
        big_k = 2.0 * np.pi / PERIOD
        g = greens_function(x, z, k=K_T, alpha0=ALPHA0, period=PERIOD, terms=3000)
        g_inf = (
            np.log1p(-np.exp(big_k * (1j * x)))
            * np.exp(ALPHA0 * (1j * x))
            + np.log1p(-np.exp(-big_k * (1j * x)))
            * np.exp(ALPHA0 * (1j * x))
        ) / (4.0 * np.pi)
        assert complex((g - g_inf)[0]) == pytest.approx(expected, abs=1e-6)

    def test_satisfies_helmholtz_away_from_the_source(self):
        """(d2/dx2 + d2/dy2 + k^2) G = 0 off the source line, by central
        finite differences on the accelerated form."""
        h = 1e-3
        x0, z0 = 170.0, 40.0

        def g(x, z):
            return greens_function(
                np.array([x]), np.array([z]),
                k=K_T, alpha0=ALPHA0, period=PERIOD, terms=200,
            )[0]

        centre = g(x0, z0)
        lap = (
            g(x0 + h, z0) + g(x0 - h, z0) + g(x0, z0 + h) + g(x0, z0 - h)
            - 4.0 * centre
        ) / h**2
        residual = lap + K_T**2 * centre
        assert abs(residual) / abs(K_T**2 * centre) < 1e-4


class TestNeumannKernel:
    def test_accelerated_matches_brute_force_off_the_boundary(self):
        x = np.array([22.0, 180.0, -310.0])
        z = np.array([48.0, -66.0, 52.0])
        nx = np.array([0.3, -0.5, 0.1])
        ny = np.sqrt(1.0 - nx**2)
        slow = brute_force_neumann(
            x, z, nx, ny, k=K_T, alpha0=ALPHA0, period=PERIOD, terms=6000
        )
        fast = neumann_function(
            x, z, nx, ny, k=K_T, alpha0=ALPHA0, period=PERIOD, terms=60
        )
        assert np.allclose(fast, slow, rtol=0.0, atol=1e-9)

    def test_is_the_normal_derivative_of_g(self):
        """N must equal the directional derivative of G along the source
        normal -- checked by finite differences, which ties the analytic
        kernel to the function it claims to differentiate."""
        h = 1e-4
        x, z = 140.0, 35.0
        nx, ny = 0.6, 0.8

        def g(dx, dz):
            return greens_function(
                np.array([x - dx]), np.array([z - dz]),
                k=K_T, alpha0=ALPHA0, period=PERIOD, terms=300,
            )[0]

        derivative = (g(h * nx, h * ny) - g(-h * nx, -h * ny)) / (2.0 * h)
        kernel = neumann_function(
            np.array([x]), np.array([z]), np.array([nx]), np.array([ny]),
            k=K_T, alpha0=ALPHA0, period=PERIOD, terms=300,
        )[0]
        assert kernel == pytest.approx(derivative, rel=1e-6)

    def test_vanishes_on_a_flat_boundary(self):
        """Flat mirror: vertical normals, zero height difference -- the
        double-layer kernel is identically zero, which is what makes the
        TM flat-mirror solution exactly psi = 2 psi_i."""
        x = np.array([50.0, 200.0, 470.0])
        z = np.zeros_like(x)
        kernel = neumann_function(
            x, z, np.zeros_like(x), np.ones_like(x),
            k=K_T, alpha0=ALPHA0, period=PERIOD, terms=100,
        )
        assert np.allclose(kernel, 0.0, atol=1e-12)


class TestComplexWavenumber:
    """The metal-side kernels of the finite-conductivity system are this
    same code with ``k`` scaled by a complex refractive index. Off the
    boundary the raw spectral sum still converges exponentially (every term
    decays on the pinned branch), so it stays a genuine reference."""

    X = np.array([13.0, 150.0, -290.0, 411.0])
    Z = np.array([55.0, 47.0, 61.0, -80.0])

    @pytest.mark.parametrize(
        "index", [0.9 + 0.2j, 0.2 + 3.0j], ids=["xray-like", "visible-gold"]
    )
    def test_greens_matches_brute_force(self, index):
        k = K_T * index
        slow = brute_force_greens(
            self.X, self.Z, k=k, alpha0=ALPHA0, period=PERIOD, terms=400
        )
        fast = greens_function(
            self.X, self.Z, k=k, alpha0=ALPHA0, period=PERIOD, terms=80
        )
        assert np.allclose(fast, slow, rtol=0.0, atol=1e-9)

    @pytest.mark.parametrize(
        "index", [0.9 + 0.2j, 0.2 + 3.0j], ids=["xray-like", "visible-gold"]
    )
    def test_neumann_matches_brute_force(self, index):
        k = K_T * index
        nx = np.array([0.3, -0.5, 0.1, 0.4])
        ny = np.sqrt(1.0 - nx**2)
        slow = brute_force_neumann(
            self.X, self.Z, nx, ny, k=k, alpha0=ALPHA0, period=PERIOD, terms=6000
        )
        fast = neumann_function(
            self.X, self.Z, nx, ny, k=k, alpha0=ALPHA0, period=PERIOD, terms=80
        )
        assert np.allclose(fast, slow, rtol=0.0, atol=1e-8)


class TestWavenumbers:
    def test_branch_is_upper_half_plane(self):
        m = np.arange(-40, 41)
        _, gamma = wavenumbers(K_T, ALPHA0, PERIOD, m)
        assert (gamma.imag >= -1e-15).all()

    def test_propagating_orders_are_real(self):
        alpha, gamma = wavenumbers(K_T, ALPHA0, PERIOD, np.array([0, 1]))
        assert np.allclose(gamma.imag, 0.0)
        assert (gamma.real > 0.0).all()

    def test_complex_k_keeps_decaying_branch(self):
        """The finite-conductivity seam: a lossy medium's wavenumber must
        still give upward-decaying terms."""
        _, gamma = wavenumbers(
            K_T * (0.99 + 0.05j), ALPHA0, PERIOD, np.arange(-30, 31)
        )
        assert (gamma.imag >= -1e-15).all()

    def test_strongly_complex_k_keeps_decaying_branch(self):
        """A visible-metal wavenumber (n ~ 0.2 + 3i) lands numpy's sqrt on
        the wrong sheet for some orders; the pin must catch every one."""
        _, gamma = wavenumbers(
            K_T * (0.2 + 3.0j), ALPHA0, PERIOD, np.arange(-60, 61)
        )
        assert (gamma.imag >= -1e-15).all()
