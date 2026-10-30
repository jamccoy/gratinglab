"""Littrow: littrow_alpha and Illumination.littrow, held to PyXFocus's grating.py.

The PyXFocus functions below are transcribed from ``jamccoy/PyXFocus_GUI``
``grating.py`` at ``62454cd`` and are the oracle, not the implementation. Its
frame is not ours: ``inc`` is the graze onto the grating SURFACE and ``psi`` the
yaw of the beam's projection on the surface, measured from the grooves. The
mapping is pinned here from direction cosines rather than assumed, because the
obvious reading, cos(inc) = sin(gamma), is wrong by 0.97 at gamma = 1.5 deg.
"""

import numpy as np
import pytest
from numpy import cos, sin, sqrt

from gratinglab.geometry import blaze_direction, littrow_alpha, sin_beta
from gratinglab.illumination import Illumination


# -- PyXFocus grating.py, verbatim but for whitespace --------------------------


def blazeYaw(inc, wave, m, d):
    return np.arcsin(m * wave / 2 / d / cos(inc))


def blazeAngle(inc, wave, m, d):
    psi = blazeYaw(inc, wave, m, d)
    beta1 = cos(inc) * cos(psi)
    alpha1 = cos(inc) * sin(psi) - m * wave / d
    return np.arcsin(alpha1 / cos(np.arcsin(beta1)))


def blazeAngle2(inc, wave, m, d):
    a = m * wave / 2 / d
    return np.arcsin(a / sqrt(1 - cos(inc) ** 2 + a**2))


def litBetaAlpha(inc, wave, m, d):
    psi = blazeYaw(inc, wave, m, d)
    beta1 = cos(inc) * cos(psi)
    alpha1 = cos(inc) * sin(psi) - m * wave / d
    return beta1, alpha1


def eta(phi0, theta0):
    return np.arcsin(np.cos(phi0) * np.cos(theta0))


def yaw(phi0, theta0):
    return np.arctan(np.sin(theta0) / np.tan(phi0))


# -- the frame mapping ---------------------------------------------------------


def pyxfocus_angles(ill: Illumination) -> tuple[float, float]:
    """(inc, psi) for one of our illuminations, read off its direction cosines."""
    ux, uy, uz = ill.direction_cosines
    return float(np.arcsin(-uy)), float(np.arctan2(-ux, uz))


# Soft X-ray off-plane, a grazing in-between, and in-plane visible light.
MOUNTS = [
    (1, 1.55, 400.0, 1.73),
    (3, 2.4, 160.0, 1.5),
    (-2, 0.8, 600.0, 2.2),
    (1, 120.0, 1000.0, 20.0),
    (1, 632.8, 1666.7, 90.0),
    (2, 500.0, 1200.0, 60.0),
]

# blazeYaw is arcsin(sin(alpha) / cos(inc)), and in-plane cos(inc) IS sin(alpha):
# the yaw is exactly 90 deg and rounding puts the argument a hair past 1, so it
# returns NaN. A singularity in the source at the mount PyXFocus never served,
# so the forms built on blazeYaw are held off-plane only.
OFFPLANE_MOUNTS = [mount for mount in MOUNTS if mount[3] < 90.0]


class TestLittrowAlpha:
    @pytest.mark.parametrize("m, lam, p, gamma_deg", MOUNTS)
    def test_the_order_returns_along_its_own_azimuth(self, m, lam, p, gamma_deg):
        g = np.radians(gamma_deg)
        a = littrow_alpha(m, lam, p, g)
        assert np.isclose(sin_beta(m, lam, p, np.sin(a), np.sin(g)), np.sin(a))

    @pytest.mark.parametrize("m, lam, p, gamma_deg", MOUNTS)
    def test_it_is_the_facet_angle_that_blazes_there(self, m, lam, p, gamma_deg):
        a = littrow_alpha(m, lam, p, np.radians(gamma_deg))
        assert np.isclose(blaze_direction(a, a), a)

    def test_zeroth_order_is_normal_incidence(self):
        assert littrow_alpha(0, 500.0, 1000.0, np.radians(30.0)) == 0.0

    def test_order_sign_mirrors_the_angle(self):
        g = np.radians(1.5)
        assert np.isclose(littrow_alpha(-2, 2.4, 160.0, g), -littrow_alpha(2, 2.4, 160.0, g))

    @pytest.mark.parametrize("m, lam, p, gamma_deg", [
        (1, 2000.0, 1000.0, 90.0),   # exactly grazing
        (1, 3000.0, 1000.0, 90.0),
        (5, 10.0, 160.0, 1.5),       # 2 p sin(gamma) is only 8.4 nm
        (-5, 10.0, 160.0, 1.5),
    ])
    def test_refuses_an_order_that_cannot_reach_littrow(self, m, lam, p, gamma_deg):
        with pytest.raises(ValueError, match="cannot reach Littrow"):
            littrow_alpha(m, lam, p, np.radians(gamma_deg))

    @pytest.mark.parametrize("kwargs, match", [
        (dict(wavelength=1.0, period=0.0, gamma=1.0), "period"),
        (dict(wavelength=0.0, period=100.0, gamma=1.0), "wavelength"),
        (dict(wavelength=1.0, period=100.0, gamma=0.0), "gamma"),
        (dict(wavelength=1.0, period=100.0, gamma=2.0), "gamma"),
    ])
    def test_refuses_bad_inputs(self, kwargs, match):
        with pytest.raises(ValueError, match=match):
            littrow_alpha(1, **kwargs)


class TestIlluminationLittrow:
    def test_in_plane_by_default(self):
        ill = Illumination.littrow(1, 632.8, 1666.7)
        assert ill.is_in_plane
        assert np.isclose(np.sin(ill.alpha), 632.8 / (2 * 1666.7))

    def test_offplane_keeps_graze_as_gamma(self):
        ill = Illumination.littrow(1, 1.55, 400.0, graze=1.73, polarization="TM")
        assert ill.gamma_deg == 1.73
        assert ill.polarization == "TM"
        assert np.isclose(ill.alpha, littrow_alpha(1, 1.55, 400.0, np.radians(1.73)))

    def test_refusal_reaches_the_caller(self):
        with pytest.raises(ValueError, match="cannot reach Littrow"):
            Illumination.littrow(5, 10.0, 160.0, graze=1.5)


class TestAgainstPyXFocus:
    @pytest.mark.parametrize("m, lam, p, gamma_deg", OFFPLANE_MOUNTS)
    def test_blazeYaw_is_our_littrow_yaw(self, m, lam, p, gamma_deg):
        inc, psi = pyxfocus_angles(Illumination.littrow(m, lam, p, graze=gamma_deg))
        assert np.isclose(blazeYaw(inc, lam, m, p), psi, atol=1e-12)

    @pytest.mark.parametrize("m, lam, p, gamma_deg", MOUNTS)
    def test_blazeAngle2_is_our_facet_angle(self, m, lam, p, gamma_deg):
        ill = Illumination.littrow(m, lam, p, graze=gamma_deg)
        inc, _ = pyxfocus_angles(ill)
        assert np.isclose(blazeAngle2(inc, lam, m, p), ill.alpha, atol=1e-12)

    @pytest.mark.parametrize("m, lam, p, gamma_deg", OFFPLANE_MOUNTS)
    def test_blazeAngle_is_the_same_angle_of_opposite_sign(self, m, lam, p, gamma_deg):
        """It divides the diffracted ray's dispersion cosine, which PyXFocus
        signs opposite to ours, so the two forms disagree in sign only."""
        ill = Illumination.littrow(m, lam, p, graze=gamma_deg)
        inc, _ = pyxfocus_angles(ill)
        assert np.isclose(blazeAngle(inc, lam, m, p), -ill.alpha, atol=1e-10)

    @pytest.mark.parametrize("m, lam, p, gamma_deg", OFFPLANE_MOUNTS)
    def test_litBetaAlpha_is_the_diffracted_direction(self, m, lam, p, gamma_deg):
        """beta1 is the groove-axis cosine, conserved, so it is cos(gamma) and
        not the cone angle its docstring names; alpha1 is our x cosine, negated."""
        ill = Illumination.littrow(m, lam, p, graze=gamma_deg)
        inc, _ = pyxfocus_angles(ill)
        beta1, alpha1 = litBetaAlpha(inc, lam, m, p)
        assert np.isclose(beta1, np.cos(ill.gamma), atol=1e-12)
        assert np.isclose(alpha1, -np.sin(ill.alpha) * ill.sin_gamma, atol=1e-12)

    def test_cos_inc_is_not_sin_gamma(self):
        """The reading roadmap §3 once gave. Only the dispersion cosine is shared."""
        ill = Illumination.littrow(1, 1.55, 400.0, graze=1.5)
        inc, psi = pyxfocus_angles(ill)
        assert not np.isclose(np.cos(inc), ill.sin_gamma, atol=0.5)
        assert np.isclose(np.cos(inc) * np.sin(psi), ill.sin_gamma * ill.sin_alpha)

    @pytest.mark.parametrize("phi0_deg, theta0_deg", [(80.0, 10.0), (60.0, 30.0), (88.5, 1.0)])
    def test_eta_and_yaw_are_representable(self, phi0_deg, theta0_deg):
        """eta and yaw are covered, not ported. Read as a unit vector with
        components (cos phi0 sin theta0, cos phi0 cos theta0, sin phi0) in
        (dispersion, normal, groove) order, they are the graze onto the surface
        and the yaw of the beam from the grooves, the two angles pinned above."""
        phi0, theta0 = np.radians(phi0_deg), np.radians(theta0_deg)
        u = np.array([np.cos(phi0) * np.sin(theta0), -np.cos(phi0) * np.cos(theta0), np.sin(phi0)])
        ill = Illumination.from_direction_cosines(-u[0], u[1], u[2])
        inc, psi = pyxfocus_angles(ill)
        assert np.isclose(eta(phi0, theta0), inc, atol=1e-12)
        assert np.isclose(yaw(phi0, theta0), psi, atol=1e-12)
