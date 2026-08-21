"""Profiles in all four representations.

The interesting assertions are the refusals: a profile that cannot be expressed
in a representation must say so, because a silently smoothed vertical facet
produces a plausible and wrong C-method answer.
"""

import numpy as np
import pytest
from pydantic import ValidationError

from gratinglab.profiles import (
    Blazed,
    BoundaryCurve,
    FromProfileData,
    Lamellar,
    Profile,
    ProfileRepresentationError,
    Sinusoidal,
)

ALL = [
    Blazed(blaze_angle=30.0),
    Blazed(blaze_angle=30.0, antiblaze_angle=80.0),
    Lamellar(depth_fraction=0.3),
    Sinusoidal(depth_fraction=0.2),
    FromProfileData(t=(0.0, 0.4, 1.0), y=(0.0, 0.25, 0.0)),
]
SMOOTH = [Sinusoidal(depth_fraction=0.2), Blazed(blaze_angle=30.0, antiblaze_angle=80.0)]


@pytest.mark.parametrize("profile", ALL, ids=lambda p: type(p).__name__)
class TestProtocol:
    def test_satisfies_the_protocol(self, profile):
        assert isinstance(profile, Profile)

    def test_height_is_periodic(self, profile):
        t = np.linspace(0.0, 1.0, 37, endpoint=False)
        assert np.allclose(profile.height(t), profile.height(t + 1.0))
        assert np.allclose(profile.height(t), profile.height(t - 3.0))

    def test_height_spans_exactly_the_declared_depth(self, profile):
        y = profile.height(np.linspace(0.0, 1.0, 4001, endpoint=False))
        assert np.ptp(y) == pytest.approx(profile.depth, rel=2e-3)
        assert y.min() >= -1e-12

    def test_height_accepts_a_scalar(self, profile):
        assert np.isscalar(float(profile.height(0.25)))


class TestBlazedGeometry:
    def test_ideal_sawtooth_depth_is_tan_blaze(self):
        """conventions.md 9: depth = period * tan(blaze) for a vertical antiblaze."""
        for angle in (5.0, 12.0, 30.0, 45.0):
            assert Blazed(blaze_angle=angle).depth == pytest.approx(
                np.tan(np.radians(angle))
            )

    def test_two_facet_depth_matches_the_thesis_example(self):
        """Chapter-2.tex:1022 -- blaze 30 deg, antiblaze 80 deg, d=160 nm -> h < 85 nm."""
        depth_nm = Blazed(blaze_angle=30.0, antiblaze_angle=80.0).depth * 160.0
        assert depth_nm == pytest.approx(83.8, abs=0.2)
        assert depth_nm < 85.0

    def test_koh_etched_silicon_facets(self):
        """29.5 deg blaze against the 70.5 deg Si(111) plane."""
        profile = Blazed(blaze_angle=29.5, antiblaze_angle=70.5)
        assert 0.0 < profile.apex < 1.0
        assert profile.height(profile.apex) == pytest.approx(profile.depth)

    def test_apex_is_where_the_maximum_is(self):
        profile = Blazed(blaze_angle=30.0, antiblaze_angle=80.0)
        t = np.linspace(0.0, 1.0, 10001, endpoint=False)
        assert t[int(np.argmax(profile.height(t)))] == pytest.approx(
            profile.apex, abs=1e-3
        )

    def test_facet_slopes_are_the_facet_angles(self):
        profile = Blazed(blaze_angle=30.0, antiblaze_angle=80.0)
        assert profile.slope(0.5 * profile.apex) == pytest.approx(
            np.tan(np.radians(30.0))
        )
        assert profile.slope(0.5 * (1 + profile.apex)) == pytest.approx(
            -np.tan(np.radians(80.0))
        )

    def test_rejects_out_of_range_angles(self):
        for kwargs in (
            dict(blaze_angle=0.0),
            dict(blaze_angle=90.0),
            dict(blaze_angle=30.0, antiblaze_angle=91.0),
        ):
            with pytest.raises(ValidationError):
                Blazed(**kwargs)


class TestRefusals:
    """The whole point of the capability declaration."""

    def test_ideal_sawtooth_refuses_to_report_a_slope(self):
        with pytest.raises(ProfileRepresentationError, match="vertical facet"):
            Blazed(blaze_angle=30.0).slope(0.5)

    def test_lamellar_refuses_to_report_a_slope(self):
        with pytest.raises(ProfileRepresentationError, match="vertical sidewalls"):
            Lamellar(depth_fraction=0.3).slope(0.5)

    def test_undercut_profile_refuses_to_slice(self):
        undercut = FromProfileData(t=(0.0, 0.5, 0.3, 1.0), y=(0.0, 0.2, 0.4, 0.0))
        assert not undercut.is_single_valued()
        with pytest.raises(ProfileRepresentationError, match="undercut"):
            undercut.slice_layers(8)

    def test_smooth_profiles_do_report_a_slope(self):
        for profile in SMOOTH:
            assert np.isfinite(profile.slope(np.linspace(0.05, 0.95, 11))).all()


class TestUndercut:
    """An undercut boundary is representable for the integral method only.

    That asymmetry is physics, not an API limitation: scalar theory and RCWA
    genuinely cannot express a boundary that is not a height function, and an
    integral-equation solver parametrises the curve so it genuinely can.
    """

    UNDERCUT = FromProfileData(
        t=(0.0, 0.30, 0.55, 0.35, 0.60, 0.90), y=(0.0, 0.15, 0.30, 0.42, 0.50, 0.10)
    )

    def test_is_detected(self):
        assert not self.UNDERCUT.is_single_valued()

    def test_height_refuses(self):
        with pytest.raises(ProfileRepresentationError, match="undercut"):
            self.UNDERCUT.height(0.5)

    def test_slicing_refuses(self):
        with pytest.raises(ProfileRepresentationError, match="undercut"):
            self.UNDERCUT.slice_layers(8)

    def test_boundary_still_works(self):
        """The whole reason to support undercut at all."""
        curve = self.UNDERCUT.boundary(128)
        assert len(curve.t) == 128
        assert np.isfinite(curve.nx).all() and np.isfinite(curve.ny).all()
        assert np.allclose(np.hypot(curve.nx, curve.ny), 1.0)
        assert curve.arc_length > 1.0

    def test_boundary_points_are_evenly_spaced_in_arc_length(self):
        ds = self.UNDERCUT.boundary(256).ds
        assert np.std(ds) / np.mean(ds) < 0.05

    def test_single_valued_data_still_takes_the_parametrised_path(self):
        """FromProfileData overrides boundary(); it must agree with height()."""
        profile = FromProfileData(t=(0.0, 0.4, 0.7, 1.0), y=(0.0, 0.25, 0.1, 0.0))
        curve = profile.boundary(512)
        assert np.allclose(curve.y, profile.height(curve.t), atol=1e-3)


@pytest.mark.parametrize("profile", ALL, ids=lambda p: type(p).__name__)
class TestLayerSlicing:
    def test_layers_tile_the_full_depth_without_gaps(self, profile):
        layers = profile.slice_layers(16)
        assert len(layers) == 16
        assert layers[0].y_lower == pytest.approx(0.0)
        assert layers[-1].y_upper == pytest.approx(profile.depth)
        for lower, upper in zip(layers[:-1], layers[1:]):
            assert lower.y_upper == pytest.approx(upper.y_lower)

    def test_fill_factor_decreases_upward(self, profile):
        """A single-valued profile narrows with height."""
        fills = [layer.fill_factor for layer in profile.slice_layers(16)]
        assert all(a >= b - 1e-9 for a, b in zip(fills, fills[1:])), fills

    def test_sliced_volume_approaches_the_true_area(self, profile):
        """The staircase must converge on the integral of height(t)."""
        exact = np.trapezoid(
            profile.height(np.linspace(0.0, 1.0, 20001)), dx=1.0 / 20000
        )
        coarse = sum(l.fill_factor * l.thickness for l in profile.slice_layers(8))
        fine = sum(l.fill_factor * l.thickness for l in profile.slice_layers(256))
        assert abs(fine - exact) < abs(coarse - exact) + 1e-9
        assert fine == pytest.approx(exact, rel=0.02)

    def test_rejects_zero_layers(self, profile):
        with pytest.raises(ValueError, match="at least one layer"):
            profile.slice_layers(0)


class TestLamellarSlicing:
    def test_fill_factor_is_the_duty_cycle_everywhere(self):
        """A rectangular groove has the same cross-section at every height."""
        for layer in Lamellar(depth_fraction=0.3, duty_cycle=0.4).slice_layers(12):
            assert layer.fill_factor == pytest.approx(0.4, abs=1e-3)


@pytest.mark.parametrize("profile", ALL, ids=lambda p: type(p).__name__)
class TestBoundaryCurve:
    def test_returns_the_requested_number_of_points(self, profile):
        curve = profile.boundary(64)
        assert isinstance(curve, BoundaryCurve)
        assert len(curve.t) == len(curve.y) == len(curve.ds) == 64

    def test_normals_are_unit_length(self, profile):
        curve = profile.boundary(64)
        assert np.allclose(np.hypot(curve.nx, curve.ny), 1.0)

    def test_normals_point_into_the_incident_medium(self, profile):
        """Outward means +y: away from the grating material."""
        assert (profile.boundary(64).ny > 0).all()

    def test_arc_length_is_at_least_the_period(self, profile):
        """A corrugated boundary is longer than the flat period it spans."""
        assert profile.boundary(512).arc_length >= 1.0

    def test_rejects_too_few_points(self, profile):
        with pytest.raises(ValueError, match="at least 3"):
            profile.boundary(2)


class TestBoundaryAccuracy:
    def test_sinusoid_arc_length_converges_on_the_analytic_value(self):
        profile = Sinusoidal(depth_fraction=0.2)
        t = np.linspace(0.0, 1.0, 200001)
        exact = np.trapezoid(np.hypot(1.0, profile.slope(t)), t)
        assert profile.boundary(4096).arc_length == pytest.approx(exact, rel=1e-4)

    def test_normal_matches_the_analytic_normal_for_a_sinusoid(self):
        profile = Sinusoidal(depth_fraction=0.2)
        curve = profile.boundary(2048)
        expected_nx = -profile.slope(curve.t) / np.hypot(1.0, profile.slope(curve.t))
        assert np.allclose(curve.nx, expected_nx, atol=1e-4)


class TestFromProfileData:
    def test_normalises_measurements(self):
        profile = FromProfileData.from_measurements(
            x_nm=np.linspace(0.0, 160.0, 50),
            y_nm=np.linspace(0.0, 80.0, 50),
            period_nm=160.0,
        )
        assert profile.depth == pytest.approx(0.5)
        assert max(profile.t) == pytest.approx(1.0)

    def test_interpolates_linearly_between_points(self):
        profile = FromProfileData(t=(0.0, 0.5, 1.0), y=(0.0, 0.4, 0.0))
        assert profile.height(0.25) == pytest.approx(0.2)
        assert profile.height(0.75) == pytest.approx(0.2)

    def test_rejects_malformed_input(self):
        for kwargs in (
            dict(t=(0.0, 1.0), y=(0.0, 0.0)),  # too few points
            dict(t=(0.0, 0.5), y=(0.0, 0.1, 0.2)),  # length mismatch
            dict(t=(0.0, 0.5, 1.5), y=(0.0, 0.1, 0.0)),  # t out of range
        ):
            with pytest.raises((ValidationError, ValueError)):
                FromProfileData(**kwargs)

    def test_rejects_negative_period(self):
        with pytest.raises(ValueError, match="must be positive"):
            FromProfileData.from_measurements([0, 1], [0, 1], period_nm=-1.0)

    def test_reproduces_an_analytic_profile_it_was_sampled_from(self):
        """Sampling a sinusoid densely and interpolating must recover it."""
        exact = Sinusoidal(depth_fraction=0.3)
        t = np.linspace(0.0, 1.0, 401)
        sampled = FromProfileData(t=tuple(t), y=tuple(exact.height(t)))
        query = np.linspace(0.0, 1.0, 97, endpoint=False)
        assert np.allclose(sampled.height(query), exact.height(query), atol=1e-5)


class TestMeasuredBlazeAngle:
    """A measured profile may carry the facet angle somebody fitted to it."""

    def test_defaults_to_absent(self):
        t = np.linspace(0.0, 1.0, 16, endpoint=False)
        assert FromProfileData(t=tuple(t), y=tuple(0.3 * t)).blaze_angle is None

    def test_from_measurements_carries_it_through(self):
        x = np.linspace(0.0, 300.0, 32)
        profile = FromProfileData.from_measurements(
            x, 0.3 * x, period_nm=300.0, blaze_angle=16.7
        )
        assert profile.blaze_angle == pytest.approx(16.7)

    @pytest.mark.parametrize("angle", [0.0, 90.0, -5.0, 91.0])
    def test_refuses_an_angle_outside_the_open_interval(self, angle):
        """A fit landing outside (0, 90) found the wrong facet, not a shallow one."""
        t = np.linspace(0.0, 1.0, 16, endpoint=False)
        with pytest.raises(ValidationError, match="blaze_angle"):
            FromProfileData(t=tuple(t), y=tuple(0.3 * t), blaze_angle=angle)

    def test_survives_a_round_trip(self):
        t = np.linspace(0.0, 1.0, 16, endpoint=False)
        original = FromProfileData(t=tuple(t), y=tuple(0.3 * t), blaze_angle=29.5)
        restored = FromProfileData(**original.model_dump())
        assert restored == original


class TestSerialization:
    @pytest.mark.parametrize("profile", ALL, ids=lambda p: type(p).__name__)
    def test_json_round_trip(self, profile):
        recovered = type(profile).model_validate_json(profile.model_dump_json())
        assert recovered == profile

    def test_is_frozen(self):
        with pytest.raises(ValidationError):
            Blazed(blaze_angle=30.0).blaze_angle = 45.0

    def test_rejects_unknown_field(self):
        with pytest.raises(ValidationError):
            Blazed(blaze_angle=30.0, depth=0.5)
