r"""Scalar solver.

The heart of this file is analytic validation: the general Fourier integral is
checked against three exact results it does not implement -- the sawtooth
``sinc^2``, the binary phase grating, and the sinusoid's Bessel form. Getting
all three right from one code path is strong evidence the numerics are sound.
"""

import numpy as np
import pytest
from scipy.special import jv

from gratinglab.geometry import blaze_wavelength, cos_beta, is_propagating, sin_beta
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, FromProfileData, Lamellar, Sinusoidal
from gratinglab.solvers import UnsupportedConfiguration, get_solver, scalar
from gratinglab.solvers.scalar import interference_factor

UNPOL = "unpolarized"


def phase_amplitude(problem, illumination, wavelength, orders):
    r"""``k * depth * sin(gamma) * [cos(alpha) + cos(beta_m)]`` per order."""
    sines = sin_beta(
        orders,
        wavelength,
        problem.period,
        illumination.sin_alpha,
        illumination.sin_gamma,
    )
    live = is_propagating(sines)
    cosines = cos_beta(sines[live])
    phi = (
        (2 * np.pi / wavelength)
        * problem.depth
        * illumination.sin_gamma
        * (illumination.cos_alpha + cosines)
    )
    return live, phi


class TestSawtoothAgainstClosedForm:
    r"""Appendix D / ISSI eq. (15): :math:`E_m = \mathrm{sinc}^2(\phi/2 - m\pi)`."""

    CASES = [
        # (period, blaze, alpha, gamma) -- visible in-plane, then soft X-ray off-plane
        (1400.0, 30.0, 30.0, 90.0),
        (1400.0, 12.0, 5.0, 90.0),
        (160.0, 30.0, 25.0, 1.5),
        (315.15, 29.5, 28.0, 1.5),
    ]

    @pytest.mark.parametrize("period,blaze,alpha,gamma", CASES)
    def test_matches_sinc_squared(self, period, blaze, alpha, gamma):
        problem = Problem(period=period, profile=Blazed(blaze_angle=blaze))
        ill = Illumination(alpha_deg=alpha, gamma_deg=gamma, polarization=UNPOL)
        wavelengths = np.linspace(0.02 * period, 0.5 * period, 9)

        scan = scalar.solve(problem, ill, wavelengths, quadrature_points=8192)

        for row, wavelength in enumerate(wavelengths):
            live, phi = phase_amplitude(problem, ill, wavelength, scan.orders)
            expected = np.sinc(phi / (2 * np.pi) - scan.orders[live]) ** 2
            assert np.allclose(scan.efficiency[row][live], expected, atol=1e-8), (
                f"lambda={wavelength}, period={period}"
            )

    def test_evanescent_orders_are_zero(self):
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, np.linspace(400.0, 900.0, 11))
        assert (scan.efficiency[~scan.propagating] == 0.0).all()


class TestBinaryPhaseGratingAgainstClosedForm:
    r"""Appendix-D.tex:511.

    Derived directly: for a two-level phase grating of duty ``w`` and phase
    ``phi``, :math:`|G_m|^2 = 4\sin^2(\phi/2)\,w^2\,\mathrm{sinc}^2(mw)` for
    ``m != 0``.
    """

    @pytest.mark.parametrize("duty", [0.25, 0.5, 0.7])
    def test_matches_closed_form(self, duty):
        problem = Problem(
            period=1400.0, profile=Lamellar(depth_fraction=0.15, duty_cycle=duty)
        )
        ill = Illumination.classical(alpha=10.0, polarization=UNPOL)
        wavelengths = np.linspace(400.0, 700.0, 7)
        scan = scalar.solve(problem, ill, wavelengths, quadrature_points=16384)

        for row, wavelength in enumerate(wavelengths):
            live, phi = phase_amplitude(problem, ill, wavelength, scan.orders)
            orders = scan.orders[live]
            nonzero = orders != 0
            expected = (
                4.0
                * np.sin(phi[nonzero] / 2.0) ** 2
                * duty**2
                * np.sinc(orders[nonzero] * duty) ** 2
            )
            got = scan.efficiency[row][live][nonzero]
            assert np.allclose(got, expected, atol=2e-4), f"lambda={wavelength}"

    def test_fifty_percent_duty_suppresses_even_orders(self):
        """W/d = 50% makes sinc^2(m/2) vanish for even m != 0."""
        problem = Problem(
            period=1400.0, profile=Lamellar(depth_fraction=0.15, duty_cycle=0.5)
        )
        ill = Illumination.classical(alpha=0.0, polarization=UNPOL)
        row = scalar.solve(problem, ill, [500.0], quadrature_points=8192).at(500.0)
        for order in row.propagating_orders:
            if order != 0 and order % 2 == 0:
                assert row[int(order)] < 1e-6, f"order {order} not suppressed"


class TestSinusoidAgainstBesselForm:
    r"""A sinusoidal phase grating gives :math:`E_m = J_m^2(A/2)`.

    From Jacobi-Anger applied to :math:`\exp[-i(A/2)\cos 2\pi t]`. Independent
    of Appendix D, so it checks the machinery from a different direction.
    """

    @pytest.mark.parametrize("depth_fraction", [0.05, 0.15, 0.3])
    def test_matches_bessel_squared(self, depth_fraction):
        problem = Problem(
            period=1400.0, profile=Sinusoidal(depth_fraction=depth_fraction)
        )
        ill = Illumination.classical(alpha=8.0, polarization=UNPOL)
        wavelengths = np.linspace(450.0, 750.0, 7)
        scan = scalar.solve(problem, ill, wavelengths, quadrature_points=16384)

        for row, wavelength in enumerate(wavelengths):
            live, phi = phase_amplitude(problem, ill, wavelength, scan.orders)
            expected = jv(scan.orders[live], phi / 2.0) ** 2
            assert np.allclose(scan.efficiency[row][live], expected, atol=1e-8)


class TestBlazeBehaviour:
    def test_peak_lands_at_the_blaze_wavelength(self):
        """Ties the solver to geometry.blaze_wavelength, which is separately tested."""
        period, blaze_deg, gamma_deg = 160.0, 12.0, 1.5
        problem = Problem(period=period, profile=Blazed(blaze_angle=blaze_deg))
        # Littrow: alpha = delta maximises the facet graze angle.
        ill = Illumination.offplane(
            graze=gamma_deg, azimuth=blaze_deg, polarization=UNPOL
        )

        for order in (2, 3, 4):
            predicted = float(
                blaze_wavelength(
                    order,
                    period,
                    np.radians(blaze_deg),
                    np.radians(blaze_deg),
                    np.radians(gamma_deg),
                )
            )
            window = np.linspace(0.75 * predicted, 1.25 * predicted, 241)
            scan = scalar.solve(problem, ill, window, quadrature_points=4096)
            found = window[int(np.argmax(scan.order(order)))]
            assert found == pytest.approx(predicted, rel=0.02), f"order {order}"

    def test_efficiency_reaches_unity_at_perfect_blaze(self):
        """An ideal sawtooth at its blaze condition puts everything in one order."""
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        peak = scalar.solve(problem, ill, [700.0], quadrature_points=8192).at(700.0)
        assert peak[2] == pytest.approx(1.0, abs=1e-6)


class TestNumerics:
    def test_converges_with_quadrature_points(self):
        problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        reference = scalar.solve(problem, ill, [2.4], quadrature_points=32768).at(2.4)
        errors = [
            np.abs(
                scalar.solve(problem, ill, [2.4], quadrature_points=n).at(2.4).efficiency
                - reference.efficiency
            ).max()
            for n in (256, 1024, 4096)
        ]
        assert errors == sorted(errors, reverse=True), errors

    def test_rejects_quadrature_below_nyquist(self):
        """Too few samples cannot represent the highest propagating order."""
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        with pytest.raises(ValueError, match="Nyquist"):
            scalar.solve(problem, ill, [20.0], quadrature_points=32)

    def test_rejects_bad_wavelengths(self):
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        with pytest.raises(ValueError, match="positive"):
            scalar.solve(problem, ill, [-500.0])

    def test_accepts_a_scalar_wavelength(self):
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        assert len(scalar.solve(problem, ill, 700.0)) == 1


class TestObliquityOption:
    def test_changes_the_distribution_between_orders(self):
        """conventions.md 5: the thesis and ISSI forms are not equivalent."""
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=45.0, polarization=UNPOL)
        plain = scalar.solve(problem, ill, [600.0], obliquity=False).at(600.0)
        tilted = scalar.solve(problem, ill, [600.0], obliquity=True).at(600.0)
        assert not np.allclose(plain.efficiency, tilted.efficiency)

    def test_is_recorded_on_the_provenance(self):
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [600.0], obliquity=True)
        assert scan.provenance.notes["obliquity"] is True

    def test_default_is_the_issi_form(self):
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        assert scalar.solve(problem, ill, [600.0]).provenance.notes["obliquity"] is False


class TestProvenance:
    def test_records_method_and_knob(self):
        problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [2.4], quadrature_points=4096)
        assert scan.provenance.method == "scalar"
        assert scan.provenance.truncation == 4096
        assert scan.provenance.wall_time_s > 0
        # Never claims convergence it has not demonstrated.
        assert scan.provenance.converged is None
        assert not scan.provenance.is_defensible

    def test_warns_when_scalar_theory_is_out_of_its_regime(self):
        """lambda/period ~ 0.5 is well past where Kirchhoff theory is reliable."""
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [700.0])
        assert any("lambda/period" in w for w in scan.provenance.warnings)

    def test_does_not_warn_in_the_soft_xray_regime(self):
        """lambda/period ~ 0.015 is exactly where scalar theory is comfortable."""
        problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [2.4])
        assert not any("lambda/period" in w for w in scan.provenance.warnings)

    def test_warns_that_polarization_is_neglected(self):
        problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
        scan = scalar.solve(
            problem, Illumination.offplane(graze=1.5, azimuth=25.0), [2.4]
        )
        assert any("neglects polarization" in w for w in scan.provenance.warnings)

    def test_warns_that_results_are_relative_without_a_coating(self):
        problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [2.4])
        assert scan.provenance.notes["normalization"] == "relative"

    def test_warns_on_a_rough_surface_past_the_fraunhofer_criterion(self):
        problem = Problem(
            period=160.0, profile=Blazed(blaze_angle=30.0), roughness=5.0
        )
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [0.5])
        assert any("Fraunhofer" in w for w in scan.provenance.warnings)


class TestCapabilities:
    def test_declares_itself_not_rigorous(self):
        assert not scalar.capabilities.rigorous

    def test_declares_its_accuracy_knob(self):
        """The convergence harness reads this."""
        assert scalar.capabilities.accuracy_knob == "quadrature_points"

    def test_refuses_an_undercut_profile(self):
        undercut = FromProfileData(
            t=(0.0, 0.5, 0.3, 0.9), y=(0.0, 0.2, 0.4, 0.1)
        )
        problem = Problem(period=160.0, profile=undercut)
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        with pytest.raises(UnsupportedConfiguration, match="undercut"):
            scalar.solve(problem, ill, [2.4])

    def test_is_in_the_registry(self):
        assert get_solver("scalar") is scalar

    def test_unknown_solver_raises_with_a_helpful_list(self):
        with pytest.raises(KeyError, match="available: .*scalar"):
            get_solver("nope")


class TestInterferenceFactor:
    def test_is_unity_at_an_exact_order(self):
        """Which is why it is not applied to order efficiencies."""
        for m in (-2, 0, 1, 5):
            assert interference_factor(m * np.pi, n_grooves=100) == pytest.approx(1.0)

    def test_narrows_as_groove_count_grows(self):
        offset = 0.01
        wide = interference_factor(offset, n_grooves=10)
        narrow = interference_factor(offset, n_grooves=200)
        assert narrow < wide

    def test_is_bounded(self):
        s = np.linspace(-np.pi, np.pi, 501)
        values = interference_factor(s, n_grooves=50)
        assert (values >= 0).all() and (values <= 1.0 + 1e-12).all()


class TestMeasuredProfile:
    def test_runs_on_a_measured_boundary(self):
        """The capability the closed forms cannot provide."""
        t = np.linspace(0.0, 1.0, 200, endpoint=False)
        measured = FromProfileData(
            t=tuple(t), y=tuple(0.3 * t + 0.01 * np.sin(6 * np.pi * t))
        )
        problem = Problem(period=160.0, profile=measured)
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, np.linspace(1.0, 5.0, 21))
        assert scan.efficiency.max() > 0
        assert np.isfinite(scan.efficiency).all()

    def test_a_sampled_sawtooth_reproduces_the_analytic_sawtooth(self):
        """Sampling an ideal sawtooth densely must give the same answer."""
        ideal = Blazed(blaze_angle=30.0)
        t = np.linspace(0.0, 1.0, 4000, endpoint=False)
        sampled = FromProfileData(t=tuple(t), y=tuple(ideal.height(t)))

        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        analytic = scalar.solve(
            Problem(period=1400.0, profile=ideal), ill, [700.0], quadrature_points=4096
        )
        numeric = scalar.solve(
            Problem(period=1400.0, profile=sampled), ill, [700.0], quadrature_points=4096
        )
        assert np.allclose(analytic.efficiency, numeric.efficiency, atol=2e-3)
