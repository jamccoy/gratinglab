"""Geometry checks, including the thesis/ISSI cross-check as executable assertions."""

import numpy as np
import pytest

from gratinglab.geometry import (
    beta,
    blaze_direction,
    blaze_wavelength,
    cos_beta,
    facet_graze,
    is_propagating,
    order_range,
    sin_beta,
)


class TestGratingEquation:
    def test_zeroth_order_is_specular(self):
        """m = 0 gives beta_0 = -alpha, independent of wavelength and mount."""
        for alpha_deg in (0.0, 10.0, 25.0, 60.0, -30.0):
            for gamma_deg in (90.0, 30.0, 1.5):
                sa = np.sin(np.radians(alpha_deg))
                sg = np.sin(np.radians(gamma_deg))
                sb = sin_beta(0, wavelength=2.4, period=160.0, sin_alpha=sa, sin_gamma=sg)
                assert np.isclose(sb, -sa)

    def test_in_plane_reduces_correctly(self):
        """sin_gamma = 1 recovers the familiar in-plane equation."""
        sa = np.sin(np.radians(10.0))
        sb = sin_beta(1, wavelength=500.0, period=1400.0, sin_alpha=sa, sin_gamma=1.0)
        assert np.isclose(sb, 500.0 / 1400.0 - sa)

    def test_off_plane_scales_as_csc_gamma(self):
        """The off-plane case is the in-plane case at lambda -> lambda * csc(gamma).

        This is the substitution introduced at Appendix-D.tex:446, and it is the
        mechanism that lets one scalar module serve both mounts.
        """
        sa, gamma_deg, lam, period = 0.3, 1.5, 2.4, 160.0
        sg = np.sin(np.radians(gamma_deg))
        off = sin_beta(3, lam, period, sa, sin_gamma=sg)
        equiv = sin_beta(3, lam / sg, period, sa, sin_gamma=1.0)
        assert np.isclose(off, equiv)

    def test_rejects_bad_inputs(self):
        for kwargs in (
            dict(period=-1.0, wavelength=1.0, sin_alpha=0.0),
            dict(period=1.0, wavelength=0.0, sin_alpha=0.0),
        ):
            with pytest.raises(ValueError):
                sin_beta(1, **kwargs)
        with pytest.raises(ValueError):
            sin_beta(1, wavelength=1.0, period=1.0, sin_alpha=0.0, sin_gamma=0.0)


class TestOrderBookkeeping:
    def test_always_contains_zero(self):
        orders = order_range(wavelength=2.4, period=160.0, sin_alpha=0.42, sin_gamma=0.026)
        assert 0 in orders

    def test_every_returned_order_propagates(self):
        """The defining property: order_range returns exactly the propagating set."""
        rng = np.random.default_rng(0)
        for _ in range(200):
            lam = rng.uniform(0.5, 700.0)
            period = rng.uniform(100.0, 2000.0)
            sa = rng.uniform(-0.99, 0.99)
            sg = rng.uniform(0.02, 1.0)
            orders = order_range(lam, period, sa, sg)
            assert is_propagating(sin_beta(orders, lam, period, sa, sg)).all()
            # and the neighbours just outside the range do not
            for edge in (orders[0] - 1, orders[-1] + 1):
                assert not is_propagating(sin_beta(edge, lam, period, sa, sg))

    def test_is_contiguous_and_ascending(self):
        orders = order_range(2.4, 315.15, 0.3, 0.026)
        assert (np.diff(orders) == 1).all()

    def test_fewer_orders_at_longer_wavelength(self):
        common = dict(period=1400.0, sin_alpha=0.2, sin_gamma=1.0)
        assert len(order_range(400.0, **common)) >= len(order_range(700.0, **common))


class TestBranches:
    def test_reflected_has_positive_cos_beta(self):
        assert cos_beta(0.5) > 0

    def test_transmitted_has_negative_cos_beta(self):
        """docs/conventions.md 4: we do not flip sin(beta); cos(beta) carries the branch."""
        assert cos_beta(0.5, transmitted=True) < 0

    def test_evanescent_is_nan_not_an_exception(self):
        assert np.isnan(cos_beta(1.5))
        assert np.isnan(beta(1.5))

    def test_beta_inverts_sin_beta(self):
        s = np.linspace(-0.99, 0.99, 51)
        assert np.allclose(np.sin(beta(s)), s)


class TestBlaze:
    """Cross-check between the McCoy thesis and the ISSI chapter.

    The two references write the blaze wavelength differently:

        thesis (Appendix-D.tex:666)  lambda_b = 2 p sin(gamma) sin(delta) cos(delta - alpha) / m
        ISSI   (following eq. 15)    m lambda_b = 2 p sin(zeta) sin(delta),
                                     sin(zeta) = sin(gamma) cos(delta - alpha)

    These are algebraically identical. Encoding that as a test means the
    agreement cannot silently rot.
    """

    CASES = [
        # (period_nm, blaze_deg, alpha_deg, gamma_deg) -- soft X-ray and visible
        (160.0, 30.0, 25.0, 1.5),
        (315.15, 29.5, 28.0, 1.5),
        (1400.0, 30.0, 30.0, 90.0),
        (160.0, 12.0, 5.0, 2.2),
    ]

    @pytest.mark.parametrize("period,blaze_deg,alpha_deg,gamma_deg", CASES)
    def test_thesis_and_issi_blaze_wavelengths_agree(
        self, period, blaze_deg, alpha_deg, gamma_deg
    ):
        delta, alpha, gamma = np.radians([blaze_deg, alpha_deg, gamma_deg])
        m = np.arange(1, 8)

        issi = blaze_wavelength(m, period, delta, alpha, gamma)
        thesis = 2 * period * np.sin(gamma) * np.sin(delta) * np.cos(delta - alpha) / m

        assert np.allclose(issi, thesis)

    @pytest.mark.parametrize("period,blaze_deg,alpha_deg,gamma_deg", CASES)
    def test_blaze_wavelength_satisfies_grating_equation_at_beta_b(
        self, period, blaze_deg, alpha_deg, gamma_deg
    ):
        """The blaze wavelength must diffract into beta_b = 2*delta - alpha.

        This ties blaze_wavelength, blaze_direction and sin_beta together: if any
        one of them drifts, this fails.
        """
        delta, alpha, gamma = np.radians([blaze_deg, alpha_deg, gamma_deg])
        beta_b = blaze_direction(delta, alpha)

        for m in range(1, 6):
            lam_b = float(blaze_wavelength(m, period, delta, alpha, gamma))
            got = sin_beta(m, lam_b, period, np.sin(alpha), np.sin(gamma))
            assert np.isclose(got, np.sin(beta_b)), f"order {m}"

    def test_facet_graze_matches_thesis_chapter2_example(self):
        """Chapter-2.tex:1046: gamma = 1.5 deg, |delta - alpha| = 5 deg -> zeta ~ 1.49 deg."""
        zeta = facet_graze(np.radians(1.5), np.radians(30.0), np.radians(25.0))
        assert np.isclose(np.degrees(zeta), 1.49, atol=0.005)

    def test_facet_graze_never_exceeds_gamma(self):
        """sin(zeta) = sin(gamma) cos(delta - alpha) <= sin(gamma), always."""
        rng = np.random.default_rng(1)
        for _ in range(100):
            gamma = np.radians(rng.uniform(0.5, 90.0))
            delta = np.radians(rng.uniform(0.0, 45.0))
            alpha = np.radians(rng.uniform(-45.0, 45.0))
            assert facet_graze(gamma, delta, alpha) <= gamma + 1e-12

    def test_littrow_maximises_facet_graze(self):
        """At alpha = delta (Littrow), zeta = gamma -- the facet is fully illuminated."""
        gamma, delta = np.radians(1.5), np.radians(30.0)
        assert np.isclose(facet_graze(gamma, delta, delta), gamma)

    def test_zeroth_order_has_no_blaze_wavelength(self):
        assert np.isinf(blaze_wavelength(0, 160.0, 0.5, 0.4, 0.03))


def _reflect(direction, normal):
    """Specular reflection of a 2D unit vector about a unit normal."""
    return direction - 2.0 * (direction @ normal) * normal


def _projected_incident(alpha):
    r"""Incident direction projected onto the (d-hat, n-hat) plane.

    From `conventions.md` §3's :math:`\mathbf{k}_i`, dropping the common
    :math:`\sin\gamma` scale and the out-of-plane :math:`\cos\gamma\,\hat{z}`.
    """
    return np.array([-np.sin(alpha), -np.cos(alpha)])


class TestFacetHandedness:
    r"""Which physical direction the profile parameter runs.

    `Profile.height(t)` is a shape in `t`; `profiles.py` never says whether
    `+t` points along `+d-hat` or `-d-hat`, and until something drew a groove
    and a ray in the same frame, nothing had to know -- every consumer takes
    :math:`\|G_m\|^2`, where the choice cancels.

    It does not cancel in a *drawing*. These tests pin the answer,
    :math:`\hat{t} = -\hat{d}`, against the two things that constrain it:
    `blaze_direction` and `facet_graze`. See `docs/conventions.md` §3 and the
    findings entry.
    """

    #: The reference off-plane case, and the one whose numbers appear in the
    #: docs: beta_b = 2(29.5) - 25 = +34.0 degrees.
    DELTA, ALPHA, GAMMA = np.radians([29.5, 25.0, 1.5])

    #: Outward normal of the active facet under t-hat = -d-hat. `Blazed.slope`
    #: rises at +tan(delta) in `t`, so the facet descends in `x`, so its
    #: outward normal leans toward +x by delta.
    @property
    def active_normal(self):
        return np.array([np.sin(self.DELTA), np.cos(self.DELTA)])

    @property
    def mirrored_normal(self):
        """What you get by assuming `t-hat = +d-hat` instead."""
        return np.array([-np.sin(self.DELTA), np.cos(self.DELTA)])

    def test_the_active_facet_reflects_into_the_blaze_direction(self):
        reflected = _reflect(_projected_incident(self.ALPHA), self.active_normal)
        expected = blaze_direction(self.DELTA, self.ALPHA)
        # k_m = k[sin(beta) x-hat + cos(beta) y-hat], so azimuth is atan2(x, y).
        assert np.isclose(np.arctan2(reflected[0], reflected[1]), expected)

    def test_the_mirrored_facet_does_not(self):
        """Non-vacuity. Without this the test above is a tautology: it would
        pass for any normal that happened to reproduce one number."""
        reflected = _reflect(_projected_incident(self.ALPHA), self.mirrored_normal)
        expected = blaze_direction(self.DELTA, self.ALPHA)
        assert not np.isclose(np.arctan2(reflected[0], reflected[1]), expected)
        # It lands on the other branch, -(2*delta + alpha) -- -84 deg, not +34.
        assert np.isclose(reflected[0], -np.sin(2 * self.DELTA + self.ALPHA))

    def test_the_incident_ray_strikes_the_active_facet_from_the_front(self):
        """A mirrored drawing puts the beam on the anti-blaze facet, arriving
        from behind the one it is supposed to illuminate."""
        assert _projected_incident(self.ALPHA) @ self.active_normal < 0.0

    def test_the_three_dimensional_graze_onto_that_facet_is_zeta(self):
        r"""`facet_graze` must be the angle to the *same* plane the blaze
        direction implies, or the two disagree about which facet is active.

        The full 3D incident direction is
        :math:`[-\sin\alpha\sin\gamma,\ -\cos\alpha\sin\gamma,\ \cos\gamma]`;
        the facet normal has no z component, so
        :math:`\sin\zeta = -\hat{k}_i \cdot \hat{n}_f`.
        """
        incident = np.array(
            [
                -np.sin(self.ALPHA) * np.sin(self.GAMMA),
                -np.cos(self.ALPHA) * np.sin(self.GAMMA),
                np.cos(self.GAMMA),
            ]
        )
        normal = np.array([*self.active_normal, 0.0])
        assert np.isclose(
            np.arcsin(-incident @ normal), facet_graze(self.GAMMA, self.DELTA, self.ALPHA)
        )

    def test_the_mirrored_facet_gives_the_wrong_graze(self):
        """The other half of the non-vacuity check."""
        incident = np.array(
            [
                -np.sin(self.ALPHA) * np.sin(self.GAMMA),
                -np.cos(self.ALPHA) * np.sin(self.GAMMA),
                np.cos(self.GAMMA),
            ]
        )
        normal = np.array([*self.mirrored_normal, 0.0])
        assert not np.isclose(
            np.arcsin(-incident @ normal), facet_graze(self.GAMMA, self.DELTA, self.ALPHA)
        )

    @pytest.mark.parametrize("blaze_deg,alpha_deg", [(29.5, 25.0), (30.0, 5.0), (12.0, -8.0)])
    def test_blazed_rises_in_t_which_is_what_makes_it_descend_in_x(
        self, blaze_deg, alpha_deg
    ):
        """Ties the convention back to the actual profile object.

        `Blazed.slope` is positive on the active facet, i.e. it rises with
        `t`. Under `t-hat = -d-hat` that is a descent in `x`, which is the
        orientation the reflection tests above require.
        """
        from gratinglab.profiles import Blazed

        profile = Blazed(blaze_angle=blaze_deg, antiblaze_angle=70.5)
        slope_on_active_facet = float(profile.slope(profile.apex * 0.5))
        assert slope_on_active_facet > 0.0
        assert np.isclose(slope_on_active_facet, np.tan(np.radians(blaze_deg)))
