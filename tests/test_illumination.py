"""Illumination: every mount constructor must resolve to the same internal state."""

import numpy as np
import pytest
from pydantic import ValidationError

from gratinglab.illumination import Illumination


class TestConstructorAgreement:
    @pytest.mark.parametrize("theta", [0.0, 5.0, 25.0, 60.0, 85.0])
    def test_conical_at_phi_zero_is_classical(self, theta):
        """phi = 0 must reduce exactly to the in-plane mount, or the azimuth
        reference direction is defined inconsistently."""
        a = Illumination.classical(alpha=theta)
        b = Illumination.conical(theta=theta, phi=0.0)
        assert np.isclose(a.alpha_deg, b.alpha_deg)
        assert np.isclose(a.gamma_deg, b.gamma_deg)

    def test_classical_is_in_plane(self):
        assert Illumination.classical(alpha=10.0).is_in_plane
        assert np.isclose(Illumination.classical(alpha=10.0).sin_gamma, 1.0)

    def test_offplane_maps_graze_to_gamma(self):
        ill = Illumination.offplane(graze=1.5, azimuth=0.97)
        assert np.isclose(ill.gamma_deg, 1.5)
        assert np.isclose(ill.alpha_deg, 0.97)
        assert not ill.is_in_plane

    def test_conical_at_phi_ninety_puts_plane_along_grooves(self):
        """phi = 90 degrees means the plane of incidence contains the groove axis."""
        ill = Illumination.conical(theta=30.0, phi=90.0)
        assert np.isclose(ill.alpha_deg, 0.0, atol=1e-9)
        assert np.isclose(ill.gamma_deg, 60.0)


class TestDirectionCosines:
    def test_unit_length(self):
        for ill in (
            Illumination.classical(alpha=25.0),
            Illumination.offplane(graze=1.5, azimuth=0.97),
            Illumination.conical(theta=40.0, phi=35.0),
        ):
            assert np.isclose(np.linalg.norm(ill.direction_cosines), 1.0)

    def test_travels_toward_grating(self):
        """uy < 0 always: the incident wave moves toward the surface."""
        for graze in (0.5, 1.5, 30.0, 90.0):
            assert Illumination.offplane(graze=graze, azimuth=10.0).direction_cosines[1] < 0

    @pytest.mark.parametrize(
        "alpha_deg,gamma_deg",
        [(0.0, 90.0), (25.0, 1.5), (-30.0, 45.0), (60.0, 89.0), (0.97, 1.5)],
    )
    def test_round_trip(self, alpha_deg, gamma_deg):
        original = Illumination(alpha_deg=alpha_deg, gamma_deg=gamma_deg)
        ux, uy, uz = original.direction_cosines
        recovered = Illumination.from_direction_cosines(ux, uy, uz)
        assert np.isclose(recovered.alpha_deg, alpha_deg, atol=1e-9)
        assert np.isclose(recovered.gamma_deg, gamma_deg, atol=1e-9)

    def test_grazing_incidence_is_nearly_tangential(self):
        """At gamma = 1.5 deg the wave is almost parallel to the surface."""
        u = Illumination.offplane(graze=1.5, azimuth=0.97).direction_cosines
        assert abs(u[1]) < 0.03  # tiny normal component
        assert u[2] > 0.999  # almost entirely along the grooves

    def test_rejects_outgoing_direction(self):
        with pytest.raises(ValueError, match="toward the grating"):
            Illumination.from_direction_cosines(0.0, 1.0, 0.0)

    def test_rejects_zero_vector(self):
        with pytest.raises(ValueError, match="non-zero"):
            Illumination.from_direction_cosines(0.0, 0.0, 0.0)


class TestValidation:
    @pytest.mark.parametrize("gamma", [0.0, -1.0, 90.1, 180.0])
    def test_rejects_out_of_range_gamma(self, gamma):
        with pytest.raises(ValidationError):
            Illumination(alpha_deg=0.0, gamma_deg=gamma)

    @pytest.mark.parametrize("alpha", [90.0, -90.0, 120.0])
    def test_rejects_grazing_or_past_grazing_alpha(self, alpha):
        with pytest.raises(ValidationError):
            Illumination(alpha_deg=alpha, gamma_deg=90.0)

    def test_rejects_unknown_field(self):
        """extra='forbid' stops a typo'd angle from being silently ignored."""
        with pytest.raises(ValidationError):
            Illumination(alpha_deg=0.0, gamma_deg=90.0, yaw=3.0)

    def test_rejects_unknown_polarization(self):
        with pytest.raises(ValidationError):
            Illumination(alpha_deg=0.0, gamma_deg=90.0, polarization="s")

    def test_is_frozen(self):
        ill = Illumination.classical(alpha=10.0)
        with pytest.raises(ValidationError):
            ill.alpha_deg = 20.0


class TestSerialization:
    def test_json_round_trip(self):
        """The spec must survive a trip through JSON unchanged -- it is the
        interchange format for benchmark cases."""
        original = Illumination.offplane(graze=1.5, azimuth=0.97, polarization="TM")
        recovered = Illumination.model_validate_json(original.model_dump_json())
        assert recovered == original

    def test_dump_contains_no_derived_quantities(self):
        """Only the canonical state is serialised; cosines are always recomputed."""
        assert set(Illumination.classical(alpha=10.0).model_dump()) == {
            "alpha_deg",
            "gamma_deg",
            "polarization",
        }
