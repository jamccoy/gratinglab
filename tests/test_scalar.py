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
from gratinglab.solvers.scalar import _reflecting_graze, interference_factor

UNPOL = "unpolarized"

# Every closed-form check runs in BOTH mounts.
#
# In-plane cases cannot detect an error in the sin(gamma) handling, because
# sin(gamma) = 1 there makes any power of it identical. Mutation testing showed
# this concretely: breaking sin(gamma) -> sin(gamma)**2 left 34 of 37 scalar
# tests passing, and only the off-plane cases failed. Off-plane is the primary
# application, so it must be exercised by every analytic check, not a subset.
MOUNTS = [
    # (period_nm, illumination, wavelengths_nm)
    (1400.0, Illumination.classical(alpha=10.0, polarization=UNPOL),
     np.linspace(400.0, 700.0, 7)),
    (315.15, Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL),
     np.linspace(1.0, 5.0, 7)),
]
MOUNT_IDS = ["in-plane", "off-plane"]


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


def obliquity(problem, illumination, wavelength, orders):
    r"""``4 c_a c_b / (c_a + c_b)^2`` per propagating order.

    The closed forms below are expressions for :math:`|G_m|^2` -- the bare
    Fourier coefficient -- and the solver returns that times the flux
    projection, so every prediction here has to carry the factor too.

    Written out longhand rather than importing
    :func:`gratinglab.geometry.flux_obliquity`, deliberately. These are the
    project's strongest quadrature checks, and routing them through the same
    function the solver calls would make a mutation of it invisible to all
    sixteen of them. ``test_perturbation.py`` is what says the factor is the
    *right* one; this only says the solver applies the one it claims to.
    """
    sines = sin_beta(
        orders,
        wavelength,
        problem.period,
        illumination.sin_alpha,
        illumination.sin_gamma,
    )
    cosines = cos_beta(sines[is_propagating(sines)])
    c_a = illumination.cos_alpha
    return 4.0 * c_a * cosines / (c_a + cosines) ** 2


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
            expected = expected * obliquity(problem, ill, wavelength, scan.orders)
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
    @pytest.mark.parametrize("mount", MOUNTS, ids=MOUNT_IDS)
    def test_matches_closed_form(self, duty, mount):
        period, ill, wavelengths = mount
        problem = Problem(
            period=period, profile=Lamellar(depth_fraction=0.15, duty_cycle=duty)
        )
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
            expected = (
                expected * obliquity(problem, ill, wavelength, scan.orders)[nonzero]
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
    @pytest.mark.parametrize("mount", MOUNTS, ids=MOUNT_IDS)
    def test_matches_bessel_squared(self, depth_fraction, mount):
        period, ill, wavelengths = mount
        problem = Problem(
            period=period, profile=Sinusoidal(depth_fraction=depth_fraction)
        )
        scan = scalar.solve(problem, ill, wavelengths, quadrature_points=16384)

        for row, wavelength in enumerate(wavelengths):
            live, phi = phase_amplitude(problem, ill, wavelength, scan.orders)
            expected = jv(scan.orders[live], phi / 2.0) ** 2
            expected = expected * obliquity(problem, ill, wavelength, scan.orders)
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

    def test_efficiency_reaches_unity_at_perfect_blaze_in_littrow(self):
        """An ideal sawtooth at its blaze condition puts everything in one order.

        **In Littrow.** ``alpha = blaze_angle`` puts the blaze direction back
        along the incidence direction, so ``cos(beta_2) == cos(alpha)``, the
        flux obliquity is exactly 1, and unity is reachable. The companion
        below is what happens when it is not -- which is the general case, and
        is why this test grew a qualifier in M16-C.
        """
        problem = Problem(period=1400.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=30.0, polarization=UNPOL)
        peak = scalar.solve(problem, ill, [700.0], quadrature_points=8192).at(700.0)
        assert peak[2] == pytest.approx(1.0, abs=1e-6)

    def test_but_away_from_littrow_it_reaches_the_obliquity_ceiling_instead(self):
        """Off Littrow the blaze order lands somewhere else, and cannot be 1.

        A perfectly blazed facet still directs everything into one order, but
        that order leaves along a different direction from the one the light
        arrived on, and the flux through a plane parallel to the surface is
        smaller by ``4 c_a c_b / (c_a + c_b)^2``.

        Unity-at-blaze was never a physical result -- a rigorous calculation
        does not give exactly 1 there either, and gives different answers for
        TE and TM. It was a property of the unfactored ``|G_m|^2``, and it is
        the one place the flux factor is visible without going to shallow
        grooves, so it is worth pinning to the number rather than loosening
        the tolerance.
        """
        period, blaze_deg, alpha_deg, order = 315.15, 29.5, 19.99, 2
        problem = Problem(period=period, profile=Blazed(blaze_angle=blaze_deg))
        ill = Illumination(alpha_deg=alpha_deg, gamma_deg=1.25, polarization=UNPOL)

        at_blaze = float(
            blaze_wavelength(
                order,
                period,
                np.radians(blaze_deg),
                np.radians(alpha_deg),
                np.radians(1.25),
            )
        )
        peak = scalar.solve(problem, ill, [at_blaze], quadrature_points=16384)
        ceiling = obliquity(problem, ill, at_blaze, peak.orders)
        expected = float(ceiling[list(peak.orders).index(order)])

        assert peak.at(at_blaze)[order] == pytest.approx(expected, abs=1e-6)
        # And the ceiling is genuinely below 1, or this asserts nothing.
        assert 0.98 < expected < 0.999


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


class TestEnergyBehaviour:
    r"""The formulation does not conserve energy, and that is a *choice*.

    Efficiency is :math:`|G_m|^2` and nothing else -- no obliquity factor, no
    renormalisation. Because the phase carries :math:`\cos\beta_m`, the
    :math:`G_m` are not Fourier coefficients of any single function, so
    Parseval does not apply and the sum drifts from unity.

    That is kept deliberately: the symmetric :math:`\cos\alpha + \cos\beta_m`
    is exactly what makes the result *reciprocal*, and the alternative that
    conserves energy violates reciprocity instead. These tests pin the
    behaviour so the drift can never be mistaken for an implementation bug --
    and so a future "fix" that renormalises it away fails loudly.
    """

    def test_shallow_limit_conserves_energy_exactly(self):
        """The check that would catch a genuine implementation bug.

        As depth goes to zero the transmittance tends to 1, all power lands in
        order 0, and the sum must tend to 1 regardless of formulation. If this
        ever fails, the problem is the code, not the physics.
        """
        ill = Illumination.classical(alpha=10.0, polarization=UNPOL)
        for depth, tolerance in ((1e-4, 1e-5), (1e-3, 1e-4)):
            problem = Problem(
                period=1400.0, profile=Sinusoidal(depth_fraction=depth)
            )
            row = scalar.solve(
                problem, ill, [500.0], quadrature_points=16384
            ).at(500.0)
            assert row.total == pytest.approx(1.0, abs=tolerance), f"depth={depth}"

    def test_deviation_grows_with_groove_depth(self):
        """Pins the mechanism: the drift tracks phase excursion across the
        groove, not lambda/period and not the propagating-order count."""
        ill = Illumination.classical(alpha=10.0, polarization=UNPOL)
        deviations = []
        for depth in (0.001, 0.005, 0.02, 0.05, 0.10):
            problem = Problem(
                period=1400.0, profile=Sinusoidal(depth_fraction=depth)
            )
            row = scalar.solve(
                problem, ill, [500.0], quadrature_points=16384
            ).at(500.0)
            deviations.append(abs(row.total - 1.0))

        assert deviations == sorted(deviations), deviations
        assert deviations[0] < 1e-4
        assert deviations[-1] > 0.05

    def test_the_deviation_is_reported_not_hidden(self):
        problem = Problem(period=315.15, profile=Blazed(blaze_angle=29.5))
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, np.linspace(1.0, 5.0, 15))

        assert abs(scan.total - 1.0).max() > 0.01
        assert any("summed efficiency" in w for w in scan.provenance.warnings)

    def test_efficiency_is_never_renormalised(self):
        r"""Appendix D normalises by :math:`\sum_m E_m`; that is wrong and must
        not creep back. If it had, the sum would be exactly 1 by construction."""
        problem = Problem(period=315.15, profile=Blazed(blaze_angle=29.5))
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        totals = scalar.solve(problem, ill, np.linspace(1.0, 5.0, 15)).total
        assert not np.allclose(totals, 1.0)


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
        assert any("sin gamma" in w for w in scan.provenance.warnings)

    def test_warns_when_the_cone_closes_even_at_small_lambda_over_period(self):
        """The guard tests the *reduced* ratio lambda/(period sin gamma).

        This geometry used to be the "comfortable" no-warning case: at
        lambda/period ~ 0.015 it looks safely inside the regime. But at a
        1.5 deg graze the structure is probed at the reduced wavelength
        lambda/sin(gamma), and the reduced ratio is ~0.57 -- and
        ``tests/test_cross_method.py::TestDivergenceAsTheConeCloses`` measures
        the scalar-vs-integral discrepancy growing along exactly this axis.
        The relabelling of this case from comfortable to warned is the point
        of the guard, not collateral damage.
        """
        problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [2.4])
        assert any("sin gamma" in w for w in scan.provenance.warnings)

    def test_does_not_warn_in_plane_at_the_same_lambda_over_period(self):
        """The same lambda/period with sin(gamma) = 1 genuinely is in regime,
        so the pair of these tests pins the guard to the reduced ratio."""
        problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.classical(alpha=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [2.4])
        assert not any("sin gamma" in w for w in scan.provenance.warnings)

    def test_does_not_warn_off_plane_when_the_reduced_ratio_is_small(self):
        """An off-plane mount is not warned about for being off-plane -- only
        for a reduced ratio the theory cannot carry."""
        problem = Problem(period=2000.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.offplane(graze=2.5, azimuth=20.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [2.0])
        assert not any("sin gamma" in w for w in scan.provenance.warnings)

    def test_warns_that_polarization_is_neglected(self):
        problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
        scan = scalar.solve(
            problem, Illumination.offplane(graze=1.5, azimuth=25.0), [2.4]
        )
        assert any("neglects polarization" in w for w in scan.provenance.warnings)

    def test_no_coating_is_the_default_mode_not_a_warning(self):
        """No coating is the normal, expected default -- relative efficiency
        is a correct result, not a deficiency. It must not appear in
        `warnings`, which is reserved for real validity concerns (regression
        guard: this used to say 'could not be evaluated', which read as an
        error and made a fully successful run look broken)."""
        problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0))
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [2.4])
        assert scan.provenance.notes["normalization"] == "relative"
        assert not any("coating" in w for w in scan.provenance.warnings)

    def test_an_unknown_coating_raises_rather_than_being_ignored(self):
        """The M8 bug at its source. `Problem.coating` was a free-form string
        that nothing read, so `coating="unobtanium"` relabelled the result
        "absolute" while leaving every number unchanged -- a scan that lied
        about what it was. A name that cannot be resolved is a question this
        solver cannot answer, not a default to fall back on."""
        from gratinglab.materials import UnknownMaterial

        problem = Problem(
            period=160.0, profile=Blazed(blaze_angle=30.0), coating="unobtanium"
        )
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        with pytest.raises(UnknownMaterial, match="Au"):
            scalar.solve(problem, ill, [2.4])

    def test_a_known_coating_is_recorded_with_its_source(self):
        """Non-vacuity for the raise, and what makes a result able to say where
        its optical constants came from."""
        problem = Problem(
            period=160.0, profile=Blazed(blaze_angle=30.0), coating="Au"
        )
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        scan = scalar.solve(problem, ill, [2.4])
        assert scan.provenance.notes["coating"].startswith("Au")
        assert "CXRO" in scan.provenance.notes["coating"]

    def test_the_absolute_label_and_the_numbers_move_together(self):
        """The property that killed the M8 bug. The label used to be
        `"absolute" if problem.coating is not None`, so any string flipped it
        while every number stayed put. Asserting the pair is what makes the
        label mean something: a run labelled absolute must *differ* from the
        same run without a coating.

        The other half of the guard is
        `test_an_unknown_coating_raises_rather_than_being_ignored` -- a name
        that resolves to nothing can no longer reach this point at all.
        """
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        profile = Blazed(blaze_angle=30.0)
        bare = scalar.solve(Problem(period=160.0, profile=profile), ill, [2.4])
        gold = scalar.solve(
            Problem(period=160.0, profile=profile, coating="Au"), ill, [2.4]
        )

        assert bare.provenance.notes["normalization"] == "relative"
        assert gold.provenance.notes["normalization"] == "absolute"
        assert not np.array_equal(bare.efficiency, gold.efficiency)

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

class TestFittedFacetAngle:
    """A measured profile that arrived with a facet fit stops guessing.

    Without one, `_reflecting_graze` falls back to the mean surface -- which
    for a sawtooth at grazing incidence is the wrong plane by the whole blaze
    angle. The fit is the number the metrology half already computes.
    """

    @staticmethod
    def _ramp(blaze_angle=None):
        t = np.linspace(0.0, 1.0, 256, endpoint=False)
        return FromProfileData(
            t=tuple(t), y=tuple(0.5 * t), blaze_angle=blaze_angle
        )

    def test_a_fit_reaches_the_exact_facet_branch(self):
        """Same geometry as a declared Blazed gives the same graze angle."""
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        declared, _ = _reflecting_graze(
            Problem(period=315.15, profile=Blazed(blaze_angle=29.5)), ill
        )
        fitted, _ = _reflecting_graze(
            Problem(period=315.15, profile=self._ramp(29.5)), ill
        )
        assert fitted == pytest.approx(declared)

    def test_the_description_says_fitted_not_declared(self):
        """Only one of the two has an uncertainty attached to it."""
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        _, declared = _reflecting_graze(
            Problem(period=315.15, profile=Blazed(blaze_angle=29.5)), ill
        )
        _, fitted = _reflecting_graze(
            Problem(period=315.15, profile=self._ramp(29.5)), ill
        )
        assert "declared" in declared and "fitted" not in declared
        assert "fitted" in fitted and "declared" not in fitted

    def test_without_a_fit_it_still_falls_back_to_the_mean_surface(self):
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        zeta, description = _reflecting_graze(
            Problem(period=315.15, profile=self._ramp()), ill
        )
        assert "mean surface" in description
        assert zeta != pytest.approx(
            _reflecting_graze(Problem(period=315.15, profile=self._ramp(29.5)), ill)[0]
        )

    def test_it_changes_the_facet_model_but_not_the_local_one(self):
        """The default model already resolves slope across the groove cycle.

        `local` reads the profile directly and never consults the facet angle,
        so carrying a fit must move `facet` and leave `local` alone. Asserting
        both directions is what keeps the claim about this honest.
        """
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        wavelengths = [2.0, 3.0]

        def totals(profile, model):
            return scalar.solve(
                Problem(period=315.15, profile=profile, coating="Au"),
                ill, wavelengths, quadrature_points=2048, reflectivity_model=model,
            ).total

        assert not np.allclose(
            totals(self._ramp(), "facet"), totals(self._ramp(29.5), "facet")
        )
        assert np.allclose(
            totals(self._ramp(), "local"), totals(self._ramp(29.5), "local")
        )


class TestAbsoluteEfficiency:
    """M15-D: the first thing in this project that changes an efficiency value.

    `docs/theory/scalar.md` section 1 has always listed "reflectivity applied
    separately as a scale factor" as an assumption, and ended the table with
    "efficiencies are **relative**, not absolute, until a materials layer
    supplies R_F". This is that layer arriving.

    **These pin ``reflectivity_model="facet"`` specifically.** M16-D made the
    groove-resolved model the default and kept this one so a prior result stays
    reproducible; this class is what says it really is unchanged. The new
    default has its own class below.
    """

    WAVELENGTHS = np.linspace(1.0, 5.0, 9)

    def _pair(self, profile, model="facet", **problem_kwargs):
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        bare = Problem(period=315.15, profile=profile, **problem_kwargs)
        gold = Problem(period=315.15, profile=profile, coating="Au", **problem_kwargs)
        return (
            scalar.solve(bare, ill, self.WAVELENGTHS, reflectivity_model=model),
            scalar.solve(gold, ill, self.WAVELENGTHS, reflectivity_model=model),
            ill,
        )

    def test_the_absolute_result_is_the_relative_one_times_the_reflectivity(self):
        """Not a second source of truth: the factor is recomputed here from the
        materials layer directly, so if the solver ever grew its own Fresnel
        arithmetic this would drift."""
        from gratinglab.geometry import facet_graze
        from gratinglab.materials import lookup
        from gratinglab.materials.fresnel import reflectivity

        profile = Blazed(blaze_angle=29.5, antiblaze_angle=70.5)
        bare, gold, ill = self._pair(profile)

        expected = reflectivity(
            lookup("Au").n(self.WAVELENGTHS),
            facet_graze(ill.gamma, np.radians(29.5), ill.alpha),
            polarization="unpolarized",
        )
        assert gold.efficiency == pytest.approx(bare.efficiency * expected[:, None])

    def test_reflection_can_only_ever_cost_efficiency(self):
        """A passive surface cannot amplify. True order by order, not just in
        the sum."""
        bare, gold, _ = self._pair(Blazed(blaze_angle=29.5, antiblaze_angle=70.5))
        assert np.all(gold.efficiency <= bare.efficiency + 1e-15)

    def test_and_it_actually_costs_some(self):
        """Non-vacuity: a reflectivity of 1.0 everywhere would satisfy the test
        above and mean nothing happened. Au at 1.5 degrees graze keeps about
        70% across this band."""
        bare, gold, _ = self._pair(Blazed(blaze_angle=29.5, antiblaze_angle=70.5))
        ratio = gold.total / bare.total
        assert np.all(ratio < 0.8)
        assert np.all(ratio > 0.6)

    def test_the_label_finally_says_absolute(self):
        bare, gold, _ = self._pair(Blazed(blaze_angle=29.5, antiblaze_angle=70.5))
        assert bare.provenance.notes["normalization"] == "relative"
        assert gold.provenance.notes["normalization"] == "absolute"

    def test_a_blazed_profile_reflects_off_its_active_facet(self):
        _, gold, _ = self._pair(Blazed(blaze_angle=29.5, antiblaze_angle=70.5))
        assert "active facet" in gold.provenance.notes["reflectivity_graze"]
        assert "29.5" in gold.provenance.notes["reflectivity_graze"]

    @pytest.mark.parametrize(
        "profile, wording",
        [
            (Lamellar(depth_fraction=0.3, duty_cycle=0.5), "exact for its flat"),
            (Sinusoidal(depth_fraction=0.15), "approximation"),
        ],
    )
    def test_a_profile_with_no_facet_angle_uses_the_mean_surface_and_says_so(
        self, profile, wording
    ):
        """Refusing these would not be more honest -- a sinusoid at grazing
        incidence does reflect. The approximation is made and recorded, and the
        wording distinguishes the lamellar case (exact for its flat tops and
        bottoms) from the varying-slope one."""
        _, gold, _ = self._pair(profile)
        note = gold.provenance.notes["reflectivity_graze"]
        assert "mean surface" in note
        assert wording in note

    def test_the_facet_model_breaks_reciprocity_which_is_why_it_is_not_default(self):
        """The defect M16-D was written to repair, pinned so it stays visible.

        `facet` evaluates R at a graze computed from alpha alone, so swapping
        alpha and beta_m changes the factor and E_m(alpha) != E_m(beta_m). The
        violation reaches 3e-2 -- against 1e-16 for the same run with no
        coating, which is how it went unnoticed: `check_reciprocity` had only
        ever been pointed at bare problems.

        Kept as an assertion rather than a note because "the old model is
        available for reproducibility" has to include reproducing what was
        wrong with it.
        """
        from gratinglab.checks import check_reciprocity

        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        problem = Problem(
            period=315.15, profile=Sinusoidal(depth_fraction=0.15), coating="Au"
        )
        report = check_reciprocity(
            scalar, problem, ill, [2.0, 4.0], quadrature_points=8192,
            reflectivity_model="facet",
        )
        assert not report.passed
        assert report.max_violation > 1e-3

    def test_the_mean_surface_graze_is_facet_graze_with_no_tilt(self):
        """Not a separate formula. `facet_graze(gamma, 0, alpha)` is exactly
        `arcsin(|k_i . n|)`, so the flat-facet case reuses a function that is
        already covered rather than introducing an untested one."""
        from gratinglab.geometry import facet_graze

        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        from_cosines = np.arcsin(abs(ill.direction_cosines[1]))
        assert facet_graze(ill.gamma, 0.0, ill.alpha) == pytest.approx(from_cosines)

    def test_no_coating_leaves_the_numbers_exactly_alone(self):
        """The default path must be untouched by this milestone -- bitwise, not
        approximately."""
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        problem = Problem(
            period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5)
        )
        again = scalar.solve(problem, ill, self.WAVELENGTHS)
        assert np.array_equal(
            scalar.solve(problem, ill, self.WAVELENGTHS).efficiency, again.efficiency
        )
        assert "reflectivity_graze" not in again.provenance.notes

    def test_a_scan_outside_the_table_raises_rather_than_extrapolating(self):
        """The range guard reaching a caller. Au is tabulated to 6.2 nm; a scan
        to 20 nm has no optical constants behind it, and inventing some would
        make every number in the result unfounded."""
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        problem = Problem(
            period=315.15,
            profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5),
            coating="Au",
        )
        with pytest.raises(ValueError, match="tabulated over"):
            scalar.solve(problem, ill, np.linspace(10.0, 20.0, 5))

    def test_absorption_shows_up_as_an_energy_deficit(self):
        """And is correctly *not* flagged as unphysical: `check_energy_balance`
        defaults to requiring only `sum <= 1`, because a deficit is ordinary
        and power going into the material is exactly that."""
        from gratinglab.checks import check_energy_balance

        bare, gold, _ = self._pair(Blazed(blaze_angle=29.5, antiblaze_angle=70.5))
        report = check_energy_balance(gold)
        assert report.passed
        assert report.max_deficit > check_energy_balance(bare).max_deficit


class TestGrooveCycleResolvedReflectivity:
    r"""M16-D: reflectivity stops being one number for the whole period.

    The M15 model evaluated :math:`R` once, at the active-facet angle, and
    applied it to every order alike -- justified in `docs/theory/scalar.md` on
    the grounds that "scalar theory has no mechanism by which one order could
    reflect differently from its neighbour". There is such a mechanism. A groove
    whose reflectivity varies across the cycle is an **amplitude** grating as
    well as a phase grating, and the two transform together.
    """

    ILL = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
    PERIOD = 315.15
    TASTE = Blazed(blaze_angle=29.5, antiblaze_angle=70.5)

    def _solve(self, profile, model, wavelengths, coating="Au", n=16384):
        problem = Problem(period=self.PERIOD, profile=profile, coating=coating)
        return scalar.solve(
            problem, self.ILL, wavelengths,
            quadrature_points=n, reflectivity_model=model,
        )

    def test_it_is_the_default(self):
        scan = self._solve(self.TASTE, "local", [3.0])
        default = scalar.solve(
            Problem(period=self.PERIOD, profile=self.TASTE, coating="Au"),
            self.ILL, [3.0], quadrature_points=16384,
        )
        assert default.provenance.notes["reflectivity_model"] == "local"
        assert default.efficiency == pytest.approx(scan.efficiency)

    def test_an_unknown_model_is_refused(self):
        with pytest.raises(ValueError, match="reflectivity_model must be"):
            self._solve(self.TASTE, "mean-surface", [3.0])

    def test_the_local_weight_is_the_geometric_mean_of_the_two_directions(self):
        r"""The exact prediction, on the profile where it can be written down.

        An ideal sawtooth has one facet angle, so :math:`\zeta_i` and
        :math:`\zeta_{d,m}` are each a single number per order and the whole
        groove-resolved integral collapses to

        .. math:: \mathscr{E}_m = \sqrt{R(\zeta_i)R(\zeta_{d,m})}\;
                  \mathscr{E}_m^{\text{bare}}

        Recomputed here from the materials layer, so this is a closed-form check
        on the machinery rather than a restatement of it -- and it is what shows
        the reflectivity has genuinely become order-dependent, since
        :math:`\zeta_{d,m}` carries the order and :math:`\zeta_i` does not.
        """
        from gratinglab.geometry import beta, facet_graze, sin_beta
        from gratinglab.materials import lookup
        from gratinglab.materials.fresnel import reflectivity

        wavelength, tilt = 3.0, np.radians(29.5)
        profile = Blazed(blaze_angle=29.5)
        bare = scalar.solve(
            Problem(period=self.PERIOD, profile=profile),
            self.ILL, [wavelength], quadrature_points=32768,
        )
        local = self._solve(profile, "local", [wavelength], n=32768)

        gold = lookup("Au").n(wavelength)
        graze_in = facet_graze(self.ILL.gamma, tilt, self.ILL.alpha)
        r_in = float(reflectivity(gold, graze_in))

        ratios = []
        for column, order in enumerate(local.orders):
            # Orders below 1e-3 are excluded, and the reason is specific rather
            # than tolerance-shopping. The collapse above is exact only where
            # the tilt is genuinely constant, and an ideal sawtooth's vertical
            # facet has no analytic slope, so `_facet_tilt` falls back to a
            # finite difference that smears the two samples either side of the
            # drop. That perturbs the integral by a fixed small amount, which is
            # a larger *fraction* of a weak order than a strong one: at 1e-3 the
            # agreement is 2e-3, at 3.5e-4 it degrades to 6e-3.
            if bare.efficiency[0, column] < 1e-3:
                continue
            sine = float(sin_beta(order, wavelength, self.PERIOD,
                                  self.ILL.sin_alpha, self.ILL.sin_gamma))
            graze_out = float(
                np.arcsin(np.sin(self.ILL.gamma) * np.cos(tilt - float(beta(sine))))
            )
            expected = np.sqrt(r_in * float(reflectivity(gold, graze_out)))
            got = local.efficiency[0, column] / bare.efficiency[0, column]
            assert got == pytest.approx(expected, rel=2e-3), f"order {order}"
            ratios.append(got)

        # Non-vacuity: if the factor were order-independent this would be flat,
        # and the whole change would be a rename. It spans 0.72 to 0.82.
        assert max(ratios) - min(ratios) > 0.05

    def test_the_facet_model_applies_one_flat_factor_and_the_local_one_does_not(self):
        """The same statement from the other side: `facet` really is flat."""
        wavelength = 3.0
        bare = scalar.solve(
            Problem(period=self.PERIOD, profile=self.TASTE),
            self.ILL, [wavelength], quadrature_points=16384,
        )
        live = bare.efficiency[0] > 1e-4
        flat = self._solve(self.TASTE, "facet", [wavelength]).efficiency[0][live]
        resolved = self._solve(self.TASTE, "local", [wavelength]).efficiency[0][live]

        assert np.ptp(flat / bare.efficiency[0][live]) < 1e-12
        assert np.ptp(resolved / bare.efficiency[0][live]) > 0.2

    def test_the_antiblaze_facet_is_shadowed_and_the_provenance_says_how_much(self):
        """The concrete thing the M15 model was getting wrong.

        On the reference geometry the anti-blaze facet is 16.7% of the period
        and faces away from the beam entirely, yet `facet` applies the active
        facet's reflectivity across the whole period.
        """
        scan = self._solve(self.TASTE, "local", [3.0])
        assert scan.provenance.notes["shadowed_fraction"] == pytest.approx(
            1.0 - self.TASTE.apex, abs=1e-3
        )
        assert scan.provenance.notes["shadowed_fraction"] == pytest.approx(0.167, abs=0.01)
        assert "83.3% of the period" in scan.provenance.notes["reflectivity_graze"]

    def test_an_order_no_facet_can_radiate_into_is_named_not_silently_zeroed(self):
        """Zero and passed-off are different states and must not look alike."""
        scan = self._solve(self.TASTE, "local", np.linspace(1.0, 5.0, 9))
        suppressed = [w for w in scan.provenance.warnings if "faces both" in w]
        assert len(suppressed) == 1
        assert "not passing-off" in suppressed[0]
        # The orders it names really are zero, and really do still propagate.
        zeroed = scan.efficiency == 0.0
        assert (zeroed & scan.propagating).any()

    def test_a_flat_profile_has_nothing_to_resolve_and_all_models_agree(self):
        """Reduction check. A groove with no slope anywhere presents one angle,
        so the three models cannot differ -- if they do, the machinery is
        inventing structure that is not in the profile."""
        flat = Lamellar(depth_fraction=1e-9, duty_cycle=0.5)
        wavelengths = [2.0, 4.0]
        got = [self._solve(flat, m, wavelengths).efficiency
               for m in ("facet", "average", "local")]
        assert got[1] == pytest.approx(got[0], rel=1e-6)
        assert got[2] == pytest.approx(got[0], rel=1e-6)

    def test_no_coating_means_no_model_and_the_numbers_are_untouched(self):
        """`reflectivity_model` must not perturb a relative run by any route."""
        reference = None
        for model in ("facet", "average", "local"):
            scan = self._solve(self.TASTE, model, [2.0, 4.0], coating=None)
            assert "reflectivity_model" not in scan.provenance.notes
            if reference is None:
                reference = scan.efficiency
            else:
                assert (scan.efficiency == reference).all()

    @pytest.mark.parametrize("model", ["facet", "average", "local"])
    def test_the_total_can_never_exceed_the_perfectly_reflecting_one(self, model):
        """The conservation statement that survives. A passive surface cannot
        add power to the scan."""
        wavelengths = np.linspace(1.0, 5.0, 9)
        bare = scalar.solve(
            Problem(period=self.PERIOD, profile=Sinusoidal(depth_fraction=0.15)),
            self.ILL, wavelengths, quadrature_points=16384,
        )
        coated = self._solve(Sinusoidal(depth_fraction=0.15), model, wavelengths)
        assert np.all(coated.efficiency.sum(axis=1) <= bare.efficiency.sum(axis=1) + 1e-12)

    def test_but_an_individual_order_can_gain(self):
        """And this is not a defect, which is why it is asserted rather than
        guarded against.

        `TestAbsoluteEfficiency.test_reflection_can_only_ever_cost_efficiency`
        holds for `facet`, where reflection is a scalar multiplier below 1. It
        is **false** for a resolved groove: varying reflectivity across the
        cycle is an amplitude grating, and an amplitude grating diffracts into
        directions where a pure phase grating cancels. Order by order the
        result can exceed the perfectly-reflecting one; the sum, above, cannot.
        """
        wavelengths = np.linspace(1.0, 5.0, 17)
        profile = Sinusoidal(depth_fraction=0.15)
        bare = scalar.solve(
            Problem(period=self.PERIOD, profile=profile),
            self.ILL, wavelengths, quadrature_points=16384,
        )
        local = self._solve(profile, "local", wavelengths)
        live = bare.efficiency > 1e-6
        assert (local.efficiency[live] / bare.efficiency[live]).max() > 1.0

    def test_only_the_symmetrised_model_is_reciprocal(self):
        """The property that actually distinguishes the three, and the reason
        `local` is the default.

        It is tempting to read this change as "resolving the groove fixes
        reciprocity". It does not. `average` resolves the groove just as fully
        as `local` and still violates reciprocity by 1.5e-2, because its
        zeta(t) is built from alpha alone. What repairs the invariant is
        symmetrising in the *exit* direction, which only `local` does.

        Parametrising this would hide the point; the contrast is the test.
        """
        from gratinglab.checks import check_reciprocity

        problem = Problem(
            period=self.PERIOD, profile=self.TASTE, coating="Au"
        )
        violations = {
            model: check_reciprocity(
                scalar, problem, self.ILL, [2.0, 4.0],
                quadrature_points=8192, reflectivity_model=model,
            ).max_violation
            for model in ("facet", "average", "local")
        }
        assert violations["local"] < 1e-12
        assert violations["facet"] > 1e-4
        assert violations["average"] > 1e-4

    def test_and_no_model_disturbs_reciprocity_without_a_coating(self):
        """Non-vacuity for the above: the violations are the reflectivity's,
        not something the option plumbing introduced."""
        from gratinglab.checks import check_reciprocity

        problem = Problem(period=self.PERIOD, profile=self.TASTE)
        for model in ("facet", "average", "local"):
            report = check_reciprocity(
                scalar, problem, self.ILL, [2.0, 4.0],
                quadrature_points=8192, reflectivity_model=model,
            )
            assert report.passed, model

    def test_the_average_is_taken_over_the_whole_period_not_just_the_lit_part(self):
        r"""A shadowed facet reflects nothing, and must be averaged in as zero.

        On an ideal sawtooth the lit facet is flat, so :math:`R(\zeta(t))` is a
        single value across all of it and the groove-cycle mean reduces to

        .. math:: \langle R\rangle = (1 - f_{\text{shadow}})\,R(\zeta_i)

        exactly. That makes `average` / `facet` equal to the lit fraction and
        nothing else -- a closed form, on the one profile where it can be
        written down.

        Averaging over the lit part alone would give :math:`R(\zeta_i)` and
        delete the shadowing this model exists to see. The two differ by 17% on
        the reference geometry, and nothing else in the suite tells them apart.
        """
        wavelengths = np.linspace(1.0, 5.0, 5)
        facet = self._solve(self.TASTE, "facet", wavelengths)
        average = self._solve(self.TASTE, "average", wavelengths)
        lit = 1.0 - average.provenance.notes["shadowed_fraction"]

        live = facet.efficiency > 1e-5
        assert np.allclose(
            average.efficiency[live] / facet.efficiency[live], lit, rtol=1e-6
        )
        assert lit == pytest.approx(self.TASTE.apex, abs=1e-3)

    def test_the_two_square_roots_are_taken_separately(self):
        r"""Pins the branch convention, and how little it is worth.

        :math:`\sqrt{r_i}\sqrt{r_d}` and :math:`\sqrt{r_i r_d}` differ by a
        global sign wherever the product's argument wraps uniformly across the
        groove -- and a global sign cancels out of the norm-squared integral. At
        grazing incidence :math:`\arg r` is nearly constant over the groove, so
        that is exactly the situation and the two agree to machine precision.

        Asserted rather than left implicit because the reasoning is easy to get
        backwards: an earlier draft of `_local_reflected_efficiency`'s docstring
        claimed the product wrapped destructively in this regime and that the
        separate roots were load-bearing here. They are not. Where they *do*
        matter -- non-uniform wrapping at Brewster -- the solver already warns.
        """
        from gratinglab.materials import lookup
        from gratinglab.materials.fresnel import amplitude
        from gratinglab.geometry import sin_facet_graze
        from gratinglab.solvers.scalar import _facet_tilt

        t = np.linspace(0.0, 1.0, 4096, endpoint=False)
        tilt = _facet_tilt(self.TASTE, t)
        gold = lookup("Au").n(3.0)
        sin_in = sin_facet_graze(self.ILL.gamma, tilt, self.ILL.alpha)
        lit = sin_in > 0
        r_s, _ = amplitude(gold, np.arcsin(np.clip(sin_in, -1.0, 1.0)))

        # The wrap is uniform here: arg(r) spans well under a radian.
        span = np.ptp(np.angle(r_s[lit]))
        assert span < 1.0, f"arg(r) spans {span:.3f} rad; the premise fails"

        separate = np.sqrt(r_s[lit]) * np.sqrt(r_s[lit])
        product = np.sqrt(r_s[lit] * r_s[lit])
        # Equal up to one overall sign, which the norm-squared integral drops.
        ratio = separate / product
        assert np.allclose(ratio, ratio[0], atol=1e-12)
        assert abs(abs(ratio[0]) - 1.0) < 1e-12

    @pytest.mark.parametrize("model", ["average", "local"])
    def test_the_groove_average_sees_the_shadowing(self, model):
        """Both resolved models must count a shadowed facet as reflecting
        nothing. Averaging only over the lit part would delete the effect."""
        lit = self._solve(Blazed(blaze_angle=29.5), model, [3.0])
        shadowed = self._solve(self.TASTE, model, [3.0])
        assert lit.provenance.notes["shadowed_fraction"] < 0.01
        assert shadowed.provenance.notes["shadowed_fraction"] > 0.15


class TestHorizonVisibility:
    r"""M19: the masks learn to see cast shadows, opt-in.

    The facet-normal test knows only *self*-shadowing -- a facet turned away
    from the ray. It cannot see the shadow the groove apex throws across the
    trough onto surface that faces the ray perfectly well, and at the
    reference geometry that blind spot is almost entirely on the **exit**
    side: a fraction of a percent of the period toward the incident beam,
    but 10-50% of it toward individual diffracted directions.
    ``tests/test_geometry.py::TestHorizonVisible`` pins the scan itself to a
    closed form; this class pins what the solver does with it.

    ``visibility="horizon"`` stays opt-in so every published number
    reproduces bit-for-bit; the geometry says the masks are right, but like
    the geometric-mean weight they await validation against a rigorous
    finite-conductivity method.
    """

    ILL = Illumination.offplane(graze=1.25, azimuth=19.99, polarization=UNPOL)
    PERIOD = 315.15
    TASTE = Blazed(blaze_angle=29.5, antiblaze_angle=70.5)

    def _solve(self, profile, visibility, wavelengths, coating="Au", n=16384):
        problem = Problem(period=self.PERIOD, profile=profile, coating=coating)
        return scalar.solve(
            problem, self.ILL, wavelengths,
            quadrature_points=n, visibility=visibility,
        )

    def test_the_default_is_facet_normal_and_is_bit_identical(self):
        """The reproducibility guarantee: not passing the option is the same
        run it was before the option existed."""
        explicit = self._solve(self.TASTE, "facet-normal", [3.0])
        default = scalar.solve(
            Problem(period=self.PERIOD, profile=self.TASTE, coating="Au"),
            self.ILL, [3.0], quadrature_points=16384,
        )
        assert (default.efficiency == explicit.efficiency).all()
        assert default.provenance.notes["visibility"] == "facet-normal"

    def test_an_unknown_visibility_is_refused(self):
        with pytest.raises(ValueError, match="visibility must be"):
            self._solve(self.TASTE, "cast", [3.0])

    def test_the_facet_model_refuses_the_horizon_rather_than_ignoring_it(self):
        """`facet` applies one reflectivity to the whole period and has no
        per-point masks a horizon could narrow. Accepting the combination and
        doing nothing would be a silent wrong answer."""
        problem = Problem(period=self.PERIOD, profile=self.TASTE, coating="Au")
        with pytest.raises(ValueError, match="cannot act under"):
            scalar.solve(
                problem, self.ILL, [3.0],
                reflectivity_model="facet", visibility="horizon",
            )

    def test_the_incident_sliver_shows_up_in_the_shadowed_fraction(self):
        r"""The independently derived sawtooth sliver, read off the provenance.

        At :math:`\alpha = 19.99` deg the cast shadow past the trough is
        :math:`\Delta = (1-t_a)(s_a - s_r)/(s_b + s_r) \approx 0.38\%` of the
        period, on top of the 16.69% anti-blaze facet the facet-normal test
        already sees. Small on the incident side -- the point is that it is
        *nonzero* and exactly predicted.
        """
        s_b = np.tan(np.radians(29.5))
        s_a = np.tan(np.radians(70.5))
        s_r = 1.0 / np.tan(np.radians(19.99))
        sliver = (1.0 - self.TASTE.apex) * (s_a - s_r) / (s_b + s_r)

        base = self._solve(self.TASTE, "facet-normal", [3.0])
        horizon = self._solve(self.TASTE, "horizon", [3.0])
        extra = (
            horizon.provenance.notes["shadowed_fraction"]
            - base.provenance.notes["shadowed_fraction"]
        )
        assert extra == pytest.approx(sliver, rel=0.05)
        assert extra > 0.003  # non-vacuity: ~0.38% of the period

    def test_the_exit_side_is_where_the_blind_spot_was(self):
        """The blaze order at its blaze wavelength loses ~1/3 of its power to
        shadows the facet-normal test cannot see: the m=+3 exit direction
        (beta = 39 deg) is cast-shadowed over 14.7% of the period. Pinned
        loosely as the regression anchor for the measured impact."""
        base = self._solve(self.TASTE, "facet-normal", [2.226])
        horizon = self._solve(self.TASTE, "horizon", [2.226])
        assert base.at(2.226)[3] == pytest.approx(0.514, abs=0.01)
        assert horizon.at(2.226)[3] == pytest.approx(0.348, abs=0.01)

    @pytest.mark.parametrize("coating", ["Au", None])
    def test_the_horizon_keeps_reciprocity(self, coating):
        """Occlusion along a straight ray reads the same from either end, so
        the mask pair is symmetric under alpha <-> beta_m by construction --
        with or without a coating carrying Fresnel weights on top."""
        from gratinglab.checks import check_reciprocity

        problem = Problem(
            period=self.PERIOD, profile=self.TASTE, coating=coating
        )
        report = check_reciprocity(
            scalar, problem, self.ILL, [2.0, 4.0],
            quadrature_points=8192, visibility="horizon",
        )
        assert report.max_violation < 1e-12

    def test_in_plane_too(self):
        """The mount pair for anything touching the transverse geometry."""
        from gratinglab.checks import check_reciprocity

        problem = Problem(period=self.PERIOD, profile=self.TASTE, coating="Au")
        ill = Illumination.classical(alpha=25.0, polarization=UNPOL)
        report = check_reciprocity(
            scalar, problem, ill, [4.0, 5.0],
            quadrature_points=8192, visibility="horizon",
        )
        assert report.max_violation < 1e-12

    @pytest.mark.parametrize(
        "ill",
        [
            Illumination.offplane(graze=1.25, azimuth=19.99, polarization=UNPOL),
            Illumination.classical(alpha=10.0, polarization=UNPOL),
        ],
        ids=["off-plane", "in-plane"],
    )
    def test_a_profile_too_shallow_to_occlude_collapses_to_facet_normal(self, ill):
        """No slope anywhere reaches the ray's run-to-drop ratio, so the
        horizon must find nothing and the two visibilities must agree
        exactly -- the reduction that shows the option adds shadows and
        nothing else."""
        problem = Problem(
            period=self.PERIOD,
            profile=Sinusoidal(depth_fraction=0.02),
            coating="Au",
        )
        base = scalar.solve(problem, ill, [3.0], visibility="facet-normal")
        horizon = scalar.solve(problem, ill, [3.0], visibility="horizon")
        assert (base.efficiency == horizon.efficiency).all()

    def test_uncoated_horizon_is_shadowing_on_a_perfect_reflector(self):
        """The masks at unit amplitude: a relative run made comparable to the
        integral solver on grooves deep enough to shadow themselves. Opt-in
        only -- the default uncoated run stays the untouched phase integral,
        and the provenance records both the mode and the cost."""
        problem = Problem(period=self.PERIOD, profile=self.TASTE)
        bare = scalar.solve(problem, self.ILL, [2.226], quadrature_points=16384)
        shadowed = scalar.solve(
            problem, self.ILL, [2.226],
            quadrature_points=16384, visibility="horizon",
        )

        assert "visibility" not in bare.provenance.notes
        assert shadowed.provenance.notes["visibility"] == "horizon"
        assert shadowed.provenance.notes["normalization"] == "relative"
        assert shadowed.provenance.notes["shadowed_fraction"] > 0.16
        # It genuinely acts: the deep sawtooth loses real power to shadows.
        assert np.abs(bare.efficiency - shadowed.efficiency).max() > 0.1
        # And it only ever removes contributions, never invents them.
        assert shadowed.total[0] < bare.total[0]


class TestTheValidityGuardsMaterialsUnlocked:
    """M15-E. `docs/theory/scalar.md` section 7 has listed the total-external-
    reflection row as "needs a materials layer to evaluate" since it was
    written, and `Problem.roughness` fed only the Fraunhofer *warning* -- the
    factor it was added for had never been written. Both land here."""

    WAVELENGTHS = np.linspace(1.0, 5.0, 9)
    PROFILE = Blazed(blaze_angle=29.5, antiblaze_angle=70.5)

    def _solve(self, *, graze=1.5, coating="Au", **problem_kwargs):
        ill = Illumination.offplane(graze=graze, azimuth=25.0, polarization=UNPOL)
        problem = Problem(
            period=315.15, profile=self.PROFILE, coating=coating, **problem_kwargs
        )
        return scalar.solve(problem, ill, self.WAVELENGTHS)

    # -- roughness now does something ------------------------------------

    def test_roughness_costs_reflectivity(self):
        smooth = self._solve()
        rough = self._solve(roughness=0.5)
        assert np.all(rough.efficiency <= smooth.efficiency)
        assert rough.total[0] < smooth.total[0]

    def test_and_did_nothing_before_a_coating_existed(self):
        """Non-vacuity in the honest direction: with no coating there is no
        reflectivity to damp, so roughness still changes no number -- and must
        not pretend to."""
        smooth = self._solve(coating=None)
        rough = self._solve(coating=None, roughness=0.5)
        assert np.array_equal(smooth.efficiency, rough.efficiency)

    def test_the_model_is_the_callers_choice(self):
        """Both are legitimate and they differ, so the choice is exposed
        rather than buried -- see `test_fresnel.py` for where and by how much.
        """
        nc = self._solve(roughness=0.5)
        dw = scalar.solve(
            Problem(
                period=315.15, profile=self.PROFILE, coating="Au", roughness=0.5
            ),
            Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL),
            self.WAVELENGTHS,
            roughness_model="debye-waller",
        )
        assert not np.array_equal(nc.efficiency, dw.efficiency)

    def test_and_can_be_switched_off_without_losing_the_figure(self):
        """`model="none"` says "I have a roughness measurement and do not want
        it applied", rather than making the caller zero the value and discard
        it from the record."""
        none = scalar.solve(
            Problem(
                period=315.15, profile=self.PROFILE, coating="Au", roughness=0.5
            ),
            Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL),
            self.WAVELENGTHS,
            roughness_model="none",
        )
        assert np.array_equal(none.efficiency, self._solve().efficiency)

    # -- total external reflection ---------------------------------------

    def test_a_facet_past_the_critical_angle_is_reported(self):
        """The row scalar.md section 7 could not evaluate. At 8 degrees graze
        the facet sits well above Au's critical angle across this whole band,
        and reflectivity collapses -- the design is outside the regime a
        grazing-incidence instrument operates in."""
        warnings = self._solve(graze=8.0).provenance.warnings
        assert any("critical angle" in w for w in warnings)

    def test_the_warning_names_the_material_and_where_it_happens(self):
        warning = next(
            w for w in self._solve(graze=8.0).provenance.warnings
            if "critical angle" in w
        )
        assert "Au" in warning
        assert "theta_c" in warning

    def test_a_grazing_design_is_not_warned_about(self):
        """Non-vacuity: the reference geometry sits at 1.5 degrees, comfortably
        below Au's 3.1-degree critical angle at 1 nm, and must stay quiet."""
        assert not any(
            "critical angle" in w for w in self._solve().provenance.warnings
        )

    def test_without_a_coating_the_check_does_not_apply_and_says_nothing(self):
        """The M8 distinction. There is no decrement to compare against, so
        this is a check that does not apply -- not a validity concern to warn
        about. A run with no coating is not a deficient run."""
        assert not any(
            "critical angle" in w
            for w in self._solve(graze=8.0, coating=None).provenance.warnings
        )

    # -- the other two guards --------------------------------------------

    def test_the_fraunhofer_check_now_covers_every_profile(self):
        """It used to be gated on `hasattr(profile, "blaze_angle")`, so a rough
        sinusoid was never checked. An optically rough surface is rough
        whatever shape its grooves are, and it now uses the same zeta the
        reflectivity does."""
        ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization=UNPOL)
        problem = Problem(
            period=315.15,
            profile=Sinusoidal(depth_fraction=0.15),
            roughness=5.0,
        )
        scan = scalar.solve(problem, ill, [0.5])
        assert any("Fraunhofer" in w for w in scan.provenance.warnings)

    def test_the_energy_message_distinguishes_absorption_from_the_approximation(self):
        """M9's text explains the thin-element energy defect, which is about
        the *model*. With a coating, part of the deficit is ordinary absorption
        instead -- a different thing, and reporting both as the approximation
        straying would overstate the model's error."""
        note = next(
            w for w in self._solve().provenance.warnings if "summed efficiency" in w
        )
        assert "ordinary absorption" in note

    def test_and_does_not_claim_absorption_when_there_is_none(self):
        """Non-vacuity: without a coating the deficit really is all model."""
        note = next(
            w
            for w in self._solve(coating=None).provenance.warnings
            if "summed efficiency" in w
        )
        assert "ordinary absorption" not in note
