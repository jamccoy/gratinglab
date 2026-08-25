r"""Resolving power: the closed form and the line profile checking each other.

``resolving_power`` returns :math:`R = |m|N` from a two-line derivation;
``line_profile`` evaluates the finite-:math:`N` interference function
numerically. Neither is trusted on its own: the closed form pins the profile's
peak and first zeros, and the profile re-derives the Rayleigh criterion the
closed form was solved from. The conical mount is tested at the corpus geometry
(period 315.15 nm, gamma = 1.5 deg), where a formula that silently drops
``sin gamma`` is off by a factor of ~38 and cannot pass.
"""

import numpy as np
import pytest

import gratinglab.geometry
import gratinglab.solvers.scalar
from gratinglab import ResolvingPower, resolving_power
from gratinglab.geometry import interference_factor, sin_beta
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed
from gratinglab.resolution import describe, line_profile

# The corpus off-plane geometry: soft X-ray, extreme off-plane, near-Littrow.
OFFPLANE = Illumination.offplane(graze=1.5, azimuth=25.0)
IN_PLANE = Illumination.classical(alpha=10.0)


def problem(n_grooves=None, period=315.15):
    return Problem(
        period=period,
        profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5),
        n_grooves=n_grooves,
    )


class TestTheClosedForm:
    @pytest.mark.parametrize("order", [1, -1, 2, -3])
    @pytest.mark.parametrize("n", [100, 4200])
    def test_r_is_order_times_n(self, order, n):
        # 1 nm, not 2: at 2 nm the -3rd order is already evanescent in this
        # mount, and these cases must all propagate.
        rp = resolving_power(problem(n), OFFPLANE, 1.0, order)
        assert rp.r == abs(order) * n
        assert rp.delta_lambda == pytest.approx(1.0 / (abs(order) * n))
        assert rp.criterion == "rayleigh"

    def test_the_record_carries_the_question(self):
        rp = resolving_power(problem(500), IN_PLANE, 100.0, 1)
        assert isinstance(rp, ResolvingPower)
        assert (rp.order, rp.wavelength, rp.n_grooves) == (1, 100.0, 500)

    def test_r_does_not_depend_on_the_mount(self):
        # R = |m|N holds because sin(gamma) scales the line width and the
        # dispersion alike. Same problem, wildly different mounts, same R.
        r_off = resolving_power(problem(700), OFFPLANE, 2.0, 1).r
        r_in = resolving_power(problem(700), IN_PLANE, 100.0, 1).r
        assert r_off == r_in == 700


class TestTheLineProfile:
    N = 800
    WAVELENGTH = 2.0
    ORDER = 3

    def profile(self, **kwargs):
        return line_profile(
            problem(self.N), OFFPLANE, self.WAVELENGTH, self.ORDER, **kwargs
        )

    def test_the_peak_sits_at_beta_m_with_intensity_one(self):
        beta_deg, intensity = self.profile()
        center = np.degrees(
            np.arcsin(
                sin_beta(
                    self.ORDER, self.WAVELENGTH, 315.15,
                    OFFPLANE.sin_alpha, OFFPLANE.sin_gamma,
                )
            )
        )
        peak = np.argmax(intensity)
        assert intensity[peak] == pytest.approx(1.0, abs=1e-12)
        assert beta_deg[peak] == pytest.approx(center, abs=1e-9)

    def test_the_first_zero_is_at_s_equals_m_pi_plus_pi_over_n(self):
        # The analytic zero, evaluated directly: sin(beta) one first-zero
        # spacing from the center must extinguish the interference function.
        center = float(
            sin_beta(
                self.ORDER, self.WAVELENGTH, 315.15,
                OFFPLANE.sin_alpha, OFFPLANE.sin_gamma,
            )
        )
        spacing = self.WAVELENGTH / (315.15 * OFFPLANE.sin_gamma * self.N)
        for sign in (+1, -1):
            s = (
                (OFFPLANE.sin_alpha + center + sign * spacing)
                * OFFPLANE.sin_gamma
                * (315.15 * np.pi / self.WAVELENGTH)
            )
            assert interference_factor(s, self.N) < 1e-15

    def test_rayleigh_criterion_end_to_end(self):
        # Shift the wavelength by delta_lambda: the shifted line's peak must
        # land exactly on the unshifted line's first zero.
        rp = resolving_power(problem(self.N), OFFPLANE, self.WAVELENGTH, self.ORDER)
        shifted_center = float(
            sin_beta(
                self.ORDER, self.WAVELENGTH + rp.delta_lambda, 315.15,
                OFFPLANE.sin_alpha, OFFPLANE.sin_gamma,
            )
        )
        s = (
            (OFFPLANE.sin_alpha + shifted_center)
            * OFFPLANE.sin_gamma
            * (315.15 * np.pi / self.WAVELENGTH)
        )
        assert interference_factor(s, self.N) < 1e-12

    def test_the_window_scales_with_half_widths(self):
        narrow_beta, _ = self.profile(half_widths=1.0)
        wide_beta, _ = self.profile(half_widths=3.0)
        narrow_span = narrow_beta[-1] - narrow_beta[0]
        wide_span = wide_beta[-1] - wide_beta[0]
        # In sin(beta) the spans are exactly 1:3; arcsin bends them slightly.
        assert wide_span == pytest.approx(3.0 * narrow_span, rel=1e-3)

    def test_n_points_is_honoured(self):
        beta_deg, intensity = self.profile(n_points=257)
        assert len(beta_deg) == len(intensity) == 257


class TestTheConicalMount:
    def test_dropping_sin_gamma_would_fail(self):
        # At gamma = 1.5 deg the first zero sits lambda/(p sin(gamma) N) from
        # the center in sin(beta). The in-plane spacing lambda/(pN) is ~38x
        # smaller and lands well inside the central lobe, where the intensity
        # is still high -- so a formula without sin(gamma) cannot pass both.
        n, wavelength, m = 500, 2.0, 1
        center = float(
            sin_beta(m, wavelength, 315.15, OFFPLANE.sin_alpha, OFFPLANE.sin_gamma)
        )

        def intensity_at(offset):
            s = (
                (OFFPLANE.sin_alpha + center + offset)
                * OFFPLANE.sin_gamma
                * (315.15 * np.pi / wavelength)
            )
            return float(interference_factor(s, n))

        conical_spacing = wavelength / (315.15 * OFFPLANE.sin_gamma * n)
        in_plane_spacing = wavelength / (315.15 * n)
        assert intensity_at(conical_spacing) < 1e-15
        assert intensity_at(in_plane_spacing) > 0.9


class TestTheRefusals:
    def test_no_groove_count_is_refused_not_treated_as_infinite(self):
        with pytest.raises(ValueError, match="infinite-grating limit"):
            resolving_power(problem(None), OFFPLANE, 2.0, 1)

    def test_order_zero_does_not_disperse(self):
        with pytest.raises(ValueError, match="does not disperse"):
            resolving_power(problem(100), OFFPLANE, 2.0, 0)

    def test_an_evanescent_order_has_no_line(self):
        # At 300 nm off-plane at gamma = 1.5 deg, order 1 needs
        # sin(beta) = 300/(315.15 sin(1.5deg)) - sin(25deg) >> 1.
        with pytest.raises(ValueError, match="evanescent"):
            resolving_power(problem(100), OFFPLANE, 300.0, 1)

    def test_line_profile_shares_the_refusals(self):
        with pytest.raises(ValueError, match="infinite-grating limit"):
            line_profile(problem(None), OFFPLANE, 2.0, 1)


class TestTheRelocation:
    def test_scalar_still_exports_the_interference_factor(self):
        # interference_factor moved from the scalar solver to geometry; the
        # old import path must keep resolving to the same object. The module
        # is fetched from sys.modules because the package attribute
        # `gratinglab.solvers.scalar` is the registered solver *instance*,
        # which shadows the module on attribute access.
        import sys

        module = sys.modules["gratinglab.solvers.scalar"]
        assert module.interference_factor is gratinglab.geometry.interference_factor
        assert "interference_factor" in module.__all__


class TestTheDescribeHelper:
    def test_one_line_per_defined_order_and_skips_the_rest(self):
        text = describe(problem(4200), OFFPLANE, 2.0, [-1, 0, 1])
        lines = text.splitlines()
        assert len(lines) == 2  # order 0 skipped, not refused
        assert "R = 4,200" in lines[0]
        assert "N = 4200" in lines[0]

    def test_silent_when_n_grooves_is_unset(self):
        assert describe(problem(None), OFFPLANE, 2.0, [1]) == ""


class TestTheMetrologyHandoff:
    def test_a_measured_boundary_reaches_a_resolving_power(self):
        # The path the milestone exists for: an averaged AFM groove carries
        # its measured groove count through to_problem, and R follows with no
        # file format in between.
        pytest.importorskip("matplotlib")  # the metrology extra
        from gratinglab.metrology.boundary.pipeline import BoundaryProfile

        t = np.linspace(0.0, 1.0, 65)
        y = np.where(t < 0.8, t * 0.25, (1.0 - t) * 1.0)  # a sawtooth-ish groove
        y[0] = y[-1] = 0.0
        boundary = BoundaryProfile(
            x_norm=t,
            y_norm=y,
            x_avg_um=t * 0.315,
            y_avg_nm=y * 315.15,
            y_std_nm=np.zeros_like(y),
            metrics={"groove_depth": float(y.max())},
            period_nm=315.15,
            n_grooves=14,
            n_used=12,
            n_edge_rejected=2,
        )
        rp = resolving_power(boundary.to_problem(), OFFPLANE, 2.0, 3)
        assert rp.n_grooves == 12  # n_used, the grooves actually averaged
        assert rp.r == 36
