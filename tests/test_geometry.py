"""Geometry checks, including the thesis/ISSI cross-check as executable assertions."""

import numpy as np
import pytest

from gratinglab.geometry import (
    beta,
    blaze_direction,
    blaze_wavelength,
    cos_beta,
    facet_graze,
    flux_obliquity,
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

    def test_the_transmitted_branch_still_inverts_sin_beta(self):
        """`docs/conventions.md` 4 in one assertion: the branch flips the sign
        of cos(beta) and leaves sin(beta) alone. `pi - b` does that; `pi + b`
        negates sin(beta) instead, and nothing noticed -- mutation testing
        found it, because `beta(..., transmitted=True)` had no test at all
        (only `cos_beta` did)."""
        s = np.linspace(-0.99, 0.99, 51)
        b = beta(s, transmitted=True)
        assert np.allclose(np.sin(b), s)
        assert np.all(np.cos(b) < 0)

    @pytest.mark.parametrize("s", [1.0, -1.0])
    def test_an_order_exactly_at_cutoff_gets_a_real_angle(self, s):
        """The passing-off boundary, which this project cares about more than
        most: `is_propagating` includes |sin(beta)| == 1, so `beta` must too,
        or an order counted as propagating has no direction to be drawn or
        summed at."""
        assert is_propagating(s)
        assert beta(s) == pytest.approx(np.sign(s) * np.pi / 2)

    def test_and_just_past_it_does_not(self):
        """Non-vacuity for the pair above: the boundary is a boundary."""
        assert not is_propagating(1.0 + 1e-12)
        assert np.isnan(beta(1.0 + 1e-12))


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


class TestFluxObliquity:
    r"""The three properties that make ``4 c_a c_b / (c_a + c_b)^2`` admissible.

    Whether it is the *correct* factor is settled in ``test_perturbation.py``,
    against a theory derived from a different starting point. What is checked
    here is the structural half: that it cannot break the invariant the scalar
    solver is built around, and that it does something.
    """

    COSINES = [(1.0, 1.0), (0.9, 0.4), (0.5, 0.99), (0.1, 0.8), (0.98, 0.02)]

    @pytest.mark.parametrize("c_a,c_b", COSINES)
    def test_it_is_symmetric_in_the_two_directions(self, c_a, c_b):
        """Why reciprocity survives it.

        `check_reciprocity` compares E_m(alpha) against E_m(beta_m). A factor
        that is unchanged when those two swap places cannot affect the
        comparison -- which is precisely what the obliquity factor of thesis
        Appendix D, ``cos(beta_m)/cos(alpha)``, fails to be.
        """
        assert flux_obliquity(c_a, c_b) == pytest.approx(flux_obliquity(c_b, c_a))

    @pytest.mark.parametrize("c_a,c_b", COSINES)
    def test_it_never_exceeds_one(self, c_a, c_b):
        """AM-GM. So it can only ever cost efficiency, never manufacture it --
        which is what lets it improve the summed-energy bound rather than
        threaten it."""
        assert flux_obliquity(c_a, c_b) <= 1.0

    def test_it_is_exactly_one_when_the_two_directions_agree(self):
        """Littrow, specular, and every zeroth order (beta_0 = -alpha, so the
        cosines are equal whatever alpha is). Exactly, not nearly: this is why
        the shallow-groove energy identity is untouched, and why the perfect
        blaze test at Littrow still reaches unity."""
        for cosine in (1.0, 0.9, 0.5, 0.01):
            assert flux_obliquity(cosine, cosine) == 1.0

    def test_it_vanishes_as_an_order_passes_off(self):
        """cos(beta_m) -> 0 is an order leaving along the surface, carrying no
        flux through a plane parallel to it. The unfactored |G_m|^2 claimed
        such an order was as bright as any other."""
        assert flux_obliquity(0.7, 1e-9) < 1e-8
        assert flux_obliquity(0.7, 0.0) == 0.0

    def test_it_is_not_vacuous_away_from_littrow(self):
        """Non-vacuity: the factor has to actually depart from 1 somewhere, or
        applying it would be a no-op dressed up as physics. At the geometry
        `test_perturbation.py` uses in-plane it is 0.61."""
        assert flux_obliquity(np.cos(np.radians(10.0)), np.cos(np.radians(76.82))) == (
            pytest.approx(0.611, rel=1e-2)
        )

    def test_it_broadcasts_over_an_array_of_orders(self):
        """The solver hands it one cos(alpha) and a vector of cos(beta_m)."""
        cosines = np.array([0.2, 0.5, 0.9])
        got = flux_obliquity(0.5, cosines)
        assert got.shape == (3,)
        assert got[1] == 1.0


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


class TestHorizonVisible:
    r"""Cast shadows on the groove, pinned to a closed form the scan never saw.

    On an ideal sawtooth (blaze slope :math:`s_b`, anti-blaze slope
    :math:`s_a`, apex at :math:`t_a`) a transverse ray at ``angle``
    :math:`\theta` clears the apex and lands on the next blaze facet a shadow
    of width

    .. math:: \Delta = (1 - t_a)\,\frac{s_a - s_r}{s_b + s_r}, \qquad
              s_r = \cot\theta

    past the trough -- zero until :math:`s_r < s_a`, i.e. until the ray dips
    below the anti-blaze slope (:math:`\theta > 19.5` deg for the 29.5/70.5
    reference groove). The trig here is derived in the test, not imported,
    so it also pins the two silent-failure traps: the travel direction
    (``angle > 0`` moves toward +t) and the nm-vs-normalised height units --
    getting either wrong reports "no cast shadow" and nothing else.
    """

    PERIOD = 315.15
    N = 16384

    def _sawtooth(self, blaze_deg=29.5, antiblaze_deg=70.5):
        from gratinglab.problem import Problem
        from gratinglab.profiles import Blazed

        profile = Blazed(blaze_angle=blaze_deg, antiblaze_angle=antiblaze_deg)
        t = np.linspace(0.0, 1.0, self.N, endpoint=False)
        problem = Problem(period=self.PERIOD, profile=profile)
        return profile, problem.height_nm(t)

    def _expected_shadow(self, blaze_deg, antiblaze_deg, angle_deg):
        s_b = np.tan(np.radians(blaze_deg))
        s_a = np.tan(np.radians(antiblaze_deg))
        s_r = 1.0 / np.tan(np.radians(angle_deg))
        apex = (1.0 / s_b) / (1.0 / s_b + 1.0 / s_a)
        sliver = (1.0 - apex) * max(s_a - s_r, 0.0) / (s_b + s_r)
        return sliver

    @pytest.mark.parametrize("angle_deg", [19.99, 25.0, 39.0, 55.3])
    def test_the_sawtooth_shadow_matches_the_closed_form(self, angle_deg):
        from gratinglab.geometry import horizon_visible

        profile, height = self._sawtooth()
        lit = horizon_visible(height, self.PERIOD, np.radians(angle_deg))
        expected = self._expected_shadow(29.5, 70.5, angle_deg)
        # The scan shadows the anti-blaze back-slope portion the ray cannot
        # see either; subtract the facet geometry to isolate the cast sliver.
        # Behind the apex the ray at angle < 90 - antiblaze grazes down the
        # anti-blaze facet; steeper rays shadow all of it (1 - apex).
        cast_onto_blaze = (1.0 - lit.mean()) - (1.0 - profile.apex)
        assert cast_onto_blaze == pytest.approx(expected, abs=2.0 / self.N * 4)

    def test_no_shadow_at_all_until_the_ray_dips_below_the_antiblaze_slope(self):
        from gratinglab.geometry import horizon_visible

        _, height = self._sawtooth()
        lit = horizon_visible(height, self.PERIOD, np.radians(15.0))
        assert lit.all()

    def test_the_mirrored_direction_sees_the_whole_groove(self):
        """At angle = -19.99 deg the ray descends steeper than either facet
        drops away from it, so it keeps intersecting terrain and nothing is
        occluded. A travel-direction sign error makes this fail, because the
        +19.99 deg case above genuinely does shadow."""
        from gratinglab.geometry import horizon_visible

        _, height = self._sawtooth()
        assert horizon_visible(height, self.PERIOD, np.radians(-19.99)).all()

    def test_matches_a_brute_force_multi_period_scan(self):
        """The O(n) closed-form treatment of whole upstream periods against
        explicit tiling, deep into the multi-period regime (near-grazing
        exit, where one apex shadows dozens of downstream periods)."""
        from gratinglab.geometry import horizon_visible

        _, height = self._sawtooth()
        t = np.arange(self.N) / self.N
        for angle_deg in (19.99, 39.0, -39.0, 80.0, 89.0):
            angle = np.radians(angle_deg)
            cot = 1.0 / np.tan(abs(angle))
            sign = 1.0 if angle > 0 else -1.0
            u = height + sign * cot * t * self.PERIOD
            reps = 4000
            if angle > 0:
                tiles = [u - k * cot * self.PERIOD for k in range(reps, -1, -1)]
                bound = np.maximum.accumulate(np.concatenate(tiles))[-self.N:]
            else:
                tiles = [u - k * cot * self.PERIOD for k in range(0, reps + 1)]
                stacked = np.concatenate(tiles)
                bound = np.maximum.accumulate(stacked[::-1])[::-1][: self.N]
            assert (
                horizon_visible(height, self.PERIOD, angle) == (u >= bound)
            ).all(), angle_deg

    def test_a_flat_mirror_is_lit_everywhere_at_any_angle(self):
        from gratinglab.geometry import horizon_visible

        flat = np.zeros(self.N)
        for angle_deg in (0.0, 30.0, -30.0, 89.0):
            assert horizon_visible(flat, self.PERIOD, np.radians(angle_deg)).all()

    def test_straight_down_shadows_nothing(self):
        from gratinglab.geometry import horizon_visible

        _, height = self._sawtooth()
        assert horizon_visible(height, self.PERIOD, 0.0).all()
