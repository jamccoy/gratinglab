r"""Fresnel reflectivity, checked against limits rather than against itself.

There is no reference table here on purpose. Every assertion is a case where
the answer is known in closed form or forced by physics -- grazing, normal
incidence, the critical angle of a lossless material, zero roughness. A test
that compared this implementation against a stored array of its own output
would pin the arithmetic and verify nothing about the formula, which is the
distinction `docs/findings.md` records from the first mutation sweep.
"""

import numpy as np
import pytest

from gratinglab import materials
from gratinglab.materials.fresnel import (
    _normal_components,
    critical_angle,
    debye_waller,
    nevot_croce,
    reflectivity,
)

#: A lossless material with a convenient critical angle: decrement 2e-3 gives
#: exactly theta_c = sqrt(2 * 2e-3) = 0.0632 rad = 3.62 deg. No absorption, so
#: reflection below theta_c must be *exactly* total -- the sharpest available
#: check, and one a real material cannot provide.
LOSSLESS = 1.0 - 2e-3 + 0j
LOSSLESS_THETA_C = np.sqrt(2 * 2e-3)


@pytest.fixture(scope="module")
def gold():
    return materials.lookup("Au")


class TestTheGrazingLimit:
    def test_reflection_is_total_at_zero_graze(self):
        """Forced, and for any material: at zero graze the normal component of
        the incident wavevector vanishes, so the amplitude ratio is -1."""
        for n in (LOSSLESS, 1.0 - 1e-3 + 1e-4j, 0.5 + 0.5j):
            assert reflectivity(n, 0.0) == pytest.approx(1.0, abs=1e-12)

    def test_and_it_falls_away_from_there(self):
        """Non-vacuity: without this, a function returning 1.0 everywhere
        would pass the test above for every material it was handed."""
        angles = np.radians([0.0, 1.0, 2.0, 5.0, 20.0])
        R = reflectivity(1.0 - 1e-3 + 1e-4j, angles)
        assert np.all(np.diff(R) < 0)
        assert R[-1] < 1e-3


class TestTheCriticalAngle:
    def test_a_lossless_material_reflects_everything_below_it(self):
        """*Exactly* everything -- this is total external reflection, not an
        approximation to it, and absorption is the only thing that spoils it."""
        below = LOSSLESS_THETA_C * np.array([0.1, 0.5, 0.9, 0.99])
        assert reflectivity(LOSSLESS, below) == pytest.approx(1.0, abs=1e-12)

    def test_and_loses_it_immediately_above(self):
        above = LOSSLESS_THETA_C * np.array([1.01, 1.1, 1.5])
        assert np.all(reflectivity(LOSSLESS, above) < 1.0)

    def test_the_falloff_is_orders_of_magnitude_not_a_tolerance(self):
        """What makes "steep" a fact. Two critical angles out, a lossless
        material has lost more than 99% of its reflectivity."""
        assert reflectivity(LOSSLESS, 2 * LOSSLESS_THETA_C) < 0.01

    def test_critical_angle_matches_the_decrement_it_came_from(self):
        assert critical_angle(2e-3) == pytest.approx(LOSSLESS_THETA_C)

    def test_a_negative_decrement_is_refused(self):
        """It would put the real part of n above 1 and delete total external
        reflection, which is the only reason grazing optics work."""
        with pytest.raises(ValueError, match="non-negative"):
            critical_angle(-1e-3)

    def test_gold_matches_its_own_table(self, gold):
        """The two paths to theta_c -- the bare relation here and the method on
        `OpticalConstants` -- must agree, or a validity guard and the number it
        guards would be computed differently."""
        decrement, _ = gold.at(1.0)
        assert critical_angle(decrement) == pytest.approx(gold.critical_angle(1.0))


class TestNormalIncidence:
    def test_it_reproduces_the_closed_form(self):
        r""":math:`R = |(n-1)/(n+1)|^2`, the textbook result, and the one place
        the general formula collapses to something checkable by hand."""
        for n in (LOSSLESS, 1.0 - 1e-3 + 1e-4j, 2.0 + 0.1j):
            expected = float(np.abs((1 - n) / (1 + n)) ** 2)
            assert reflectivity(n, np.pi / 2, polarization="s") == pytest.approx(
                expected, rel=1e-12
            )

    def test_s_and_p_are_degenerate_there(self):
        """They must be: at normal incidence there is no plane of incidence to
        distinguish them. A formula that separated them here would be wrong."""
        for n in (LOSSLESS, 1.0 - 1e-3 + 1e-4j, 2.0 + 0.1j):
            s = reflectivity(n, np.pi / 2, polarization="s")
            p = reflectivity(n, np.pi / 2, polarization="p")
            assert s == pytest.approx(p, rel=1e-9)


class TestPolarizationInTheGrazingRegime:
    """`conventions.md` §7 claims reflectivity is "nearly polarization-
    independent" at grazing incidence in the soft X-ray, citing the thesis.
    Measured here rather than left as folklore, because M15-D will lean on it
    when it maps groove-referenced TE/TM onto these facet-local s and p."""

    def test_s_and_p_agree_to_under_a_percent_at_the_reference_geometry(self, gold):
        n = gold.n(1.0)
        graze = np.radians(1.5)
        s = float(reflectivity(n, graze, polarization="s"))
        p = float(reflectivity(n, graze, polarization="p"))
        assert abs(s - p) / s < 0.01

    def test_but_they_separate_usefully_away_from_it(self, gold):
        """Non-vacuity, and the reason §7 calls it a regime-specific
        approximation rather than a general one: the same function knows the
        difference perfectly well at 45 degrees."""
        n = gold.n(1.0)
        s = float(reflectivity(n, np.radians(45.0), polarization="s"))
        p = float(reflectivity(n, np.radians(45.0), polarization="p"))
        assert abs(s - p) / s > 0.5

    def test_unpolarized_is_the_mean_and_is_computed(self, gold):
        """`conventions.md` §7: "(TE + TM) / 2, computed, never assumed"."""
        n = gold.n(2.0)
        graze = np.radians([0.5, 3.0, 30.0])
        s = reflectivity(n, graze, polarization="s")
        p = reflectivity(n, graze, polarization="p")
        both = reflectivity(n, graze, polarization="unpolarized")
        assert both == pytest.approx(0.5 * (s + p))

    def test_an_unknown_polarization_is_refused_and_names_the_convention(self):
        """Passing "TE" here is the mistake worth catching: it is a real
        polarization name in this project and it means something else."""
        with pytest.raises(ValueError, match="TE/TM"):
            reflectivity(LOSSLESS, 0.1, polarization="TE")


class TestTheBranchOfTheSquareRoot:
    def test_the_transmitted_wave_decays(self):
        r"""The one place a sign error hides silently. Under ``exp(-i\omega t)``
        the transmitted normal component needs a non-negative imaginary part;
        the other branch describes a wave *growing* with depth and returns a
        reflectivity above one."""
        for graze_deg in (0.1, 1.0, 3.0, 10.0, 89.0):
            _, k2 = _normal_components(1.0 - 1e-3 + 1e-4j, np.radians(graze_deg))
            assert k2.imag >= 0.0

    def test_and_reflectivity_never_exceeds_one(self, gold):
        """The observable consequence, across the whole table."""
        w = np.linspace(*gold.range_nm, 40)
        for graze_deg in (0.01, 0.5, 1.5, 3.0, 10.0, 80.0):
            R = reflectivity(gold.n(w), np.radians(graze_deg))
            assert np.all(R <= 1.0 + 1e-12)
            assert np.all(R >= 0.0)


class TestRoughness:
    WAVELENGTH = 2.0
    GRAZE = np.radians(1.0)

    def _R(self, sigma, model="nevot-croce"):
        return float(
            reflectivity(
                1.0 - 1e-3 + 1e-4j,
                self.GRAZE,
                roughness_nm=sigma,
                wavelength_nm=self.WAVELENGTH,
                model=model,
            )
        )

    def test_a_perfectly_smooth_surface_is_unchanged(self):
        """Exactly unchanged, not nearly: both factors are ``exp(0)``."""
        smooth = float(reflectivity(1.0 - 1e-3 + 1e-4j, self.GRAZE))
        assert self._R(0.0) == smooth
        assert self._R(0.0, "debye-waller") == smooth

    def test_both_factors_are_exactly_one_at_zero_roughness(self):
        assert nevot_croce(1.0, 1.0, 0.0) == 1.0
        assert debye_waller(1.0, 0.0) == 1.0

    def test_roughness_only_ever_costs_reflectivity(self):
        for model in ("nevot-croce", "debye-waller"):
            values = [self._R(s, model) for s in (0.0, 0.1, 0.3, 1.0)]
            assert values == sorted(values, reverse=True), model

    def test_nevot_croce_reduces_to_debye_waller_when_there_is_no_material(self):
        r"""The limit that fixes the normalisation, and nothing else does.

        As :math:`n \to 1` the interface disappears: the transmitted wave *is*
        the incident wave, :math:`k_{tz} \to k_{iz}`, and the two roughness
        models become the same expression. Every other property either function
        has -- bounded by 1, monotone in sigma, exactly 1 at sigma = 0 -- is
        satisfied by both the amplitude form and the intensity form, so this is
        the only check in the suite that can tell them apart.

        It failed before M16-B by a factor of the exponent: the code returned
        :math:`|e^{-2k_{iz}k_{tz}\sigma^2}|` where the intensity needs its
        square. The residual 1e-5 here is the 1e-6 left in ``n``, not slack.
        """
        vacuum = 1.0 - 1e-6 + 0j
        wavelength, sigma = 2.0, 0.5
        for degrees in (5.0, 15.0, 30.0, 60.0, 85.0):
            graze = np.radians(degrees)
            common = dict(roughness_nm=sigma, wavelength_nm=wavelength)
            nc = float(reflectivity(vacuum, graze, model="nevot-croce", **common))
            dw = float(reflectivity(vacuum, graze, model="debye-waller", **common))
            assert nc == pytest.approx(dw, rel=1e-4), f"graze={degrees} deg"

    def test_debye_waller_over_damps_near_and_above_the_critical_angle(self):
        """Which is the reason both are offered, and the scope of the claim.

        Debye-Waller knows only about the incident wave, so it damps the same
        whether the field penetrates or not; Névot-Croce carries the
        transmitted component and damps less where that wave is evanescent.
        The gap is ~1e-2 on Au, which is where the choice actually matters.

        Scoped on purpose. An earlier version of this test asserted Debye-
        Waller was *always* the more pessimistic, and passed -- because it
        sampled one wavelength and one synthetic material, both near
        theta_c. See the companion below for what happens far from it.
        """
        gold = materials.lookup("Au")
        for lam in (1.0, 2.0, 4.0, 6.0):
            theta_c = float(gold.critical_angle(lam))
            for fraction in (0.8, 1.0, 1.5):
                common = dict(roughness_nm=0.5, wavelength_nm=lam)
                nc = float(reflectivity(gold.n(lam), fraction * theta_c,
                                        model="nevot-croce", **common))
                dw = float(reflectivity(gold.n(lam), fraction * theta_c,
                                        model="debye-waller", **common))
                assert dw < nc, f"lambda={lam}, zeta/theta_c={fraction}"

    def test_and_far_above_it_the_two_models_converge(self):
        """The other regime -- and M16-B moved it.

        The models agree where :math:`k_{tz} \\to k_{iz}`, which is *well above*
        theta_c, not below: that is the same limit as
        ``test_nevot_croce_reduces_to_debye_waller_when_there_is_no_material``,
        reached by opening the angle instead of by removing the material. At
        3x theta_c on Au they differ by ~1e-5.

        This test used to assert convergence *below* theta_c, at 1e-3 with "no
        reliable ordering". That was reading the un-squared Névot-Croce of
        M15-B. With the intensity form the ordering below theta_c is perfectly
        reliable, and it runs the other way -- see the companion below.
        """
        gold = materials.lookup("Au")
        for lam in (1.0, 2.0, 4.0, 6.0):
            graze = 3.0 * float(gold.critical_angle(lam))
            common = dict(roughness_nm=0.5, wavelength_nm=lam)
            nc = float(reflectivity(gold.n(lam), graze,
                                    model="nevot-croce", **common))
            dw = float(reflectivity(gold.n(lam), graze,
                                    model="debye-waller", **common))
            assert abs(dw - nc) < 1e-4, f"lambda={lam}"

    def test_while_deep_below_it_debye_waller_damps_slightly_less(self):
        """The third regime, and the direction reverses.

        Deep under theta_c the transmitted wave is strongly evanescent,
        ``k_tz`` is nearly pure imaginary, and the Névot-Croce product
        ``k_iz k_tz`` picks up a real part that Debye-Waller has no term for.
        The result is a *small* gap, ~2e-3, in the opposite direction to the
        near-theta_c case -- Debye-Waller now retains marginally more
        reflectivity.

        Worth pinning precisely because it is counterintuitive next to the
        "Debye-Waller over-damps" test above. Both are true; they describe
        different angles, and the sign of the gap is what tells them apart.
        """
        gold = materials.lookup("Au")
        for lam in (1.0, 2.0, 4.0, 6.0):
            graze = 0.15 * float(gold.critical_angle(lam))
            common = dict(roughness_nm=0.5, wavelength_nm=lam)
            nc = float(reflectivity(gold.n(lam), graze,
                                    model="nevot-croce", **common))
            dw = float(reflectivity(gold.n(lam), graze,
                                    model="debye-waller", **common))
            assert 0.0 < dw - nc < 3e-3, f"lambda={lam}: dw-nc={dw - nc:.2e}"

    def test_roughness_without_a_wavelength_is_refused(self):
        """Both factors compare sigma against a wavelength. Defaulting one
        would invent the very scale the factor is about."""
        with pytest.raises(ValueError, match="wavelength_nm is required"):
            reflectivity(LOSSLESS, self.GRAZE, roughness_nm=0.5)

    def test_but_a_smooth_surface_needs_none(self):
        """Non-vacuity for the refusal above: sigma = 0 must not demand a
        wavelength it will not use."""
        assert reflectivity(LOSSLESS, self.GRAZE, roughness_nm=0.0) > 0

    def test_negative_roughness_is_refused(self):
        with pytest.raises(ValueError, match="non-negative"):
            reflectivity(
                LOSSLESS, self.GRAZE, roughness_nm=-0.1, wavelength_nm=self.WAVELENGTH
            )

    def test_an_unknown_model_is_refused(self):
        with pytest.raises(ValueError, match="nevot-croce"):
            reflectivity(
                LOSSLESS,
                self.GRAZE,
                roughness_nm=0.5,
                wavelength_nm=self.WAVELENGTH,
                model="handwaving",
            )

    def test_model_none_skips_it_entirely(self):
        """An explicit way to say "I have a roughness figure and I do not want
        it applied", rather than making the caller zero the value and lose it."""
        smooth = float(reflectivity(LOSSLESS, self.GRAZE))
        assert (
            reflectivity(
                LOSSLESS,
                self.GRAZE,
                roughness_nm=1.0,
                wavelength_nm=self.WAVELENGTH,
                model="none",
            )
            == smooth
        )


class TestShape:
    def test_it_broadcasts_over_a_wavelength_scan(self, gold):
        """How M15-D will call it: one graze angle, a whole scan of indices."""
        w = np.linspace(1.0, 5.0, 17)
        R = reflectivity(gold.n(w), np.radians(1.5))
        assert R.shape == (17,)

    def test_and_over_a_scan_of_angles(self):
        R = reflectivity(LOSSLESS, np.radians(np.linspace(0.1, 10.0, 9)))
        assert R.shape == (9,)

    def test_a_scalar_in_gives_something_float_like_out(self):
        assert float(reflectivity(LOSSLESS, 0.01)) <= 1.0


class TestTheSolverUsesThisAndNotACopyOfIt:
    """Until M15-D this class asserted the opposite -- that the solver did not
    import this module -- to hold the boundary between the pure layer and the
    commit that changes efficiency values. M15-D crossed it deliberately, so
    the guard becomes the thing it was protecting: the solver reflects using
    *this* code, not arithmetic of its own.
    """

    def test_the_solver_imports_fresnel(self):
        import ast
        import inspect
        import sys

        import gratinglab.solvers.scalar  # noqa: F401

        # The module, not the registered singleton that shadows it.
        source = inspect.getsource(sys.modules["gratinglab.solvers.scalar"])
        imported = {
            node.module
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.ImportFrom) and node.module
        }
        assert any("fresnel" in name for name in imported)

    def test_and_computes_no_fresnel_arithmetic_of_its_own(self):
        """The property that matters. A solver that grew its own amplitude
        coefficients could drift from this module silently -- and
        `test_scalar.py::TestAbsoluteEfficiency` recomputes the factor from
        here, so the two would then disagree.
        """
        import inspect
        import sys

        import gratinglab.solvers.scalar  # noqa: F401

        source = inspect.getsource(sys.modules["gratinglab.solvers.scalar"])
        for fingerprint in ("np.sqrt(n", "cos(graze)", "** 2 * k1"):
            assert fingerprint not in source
