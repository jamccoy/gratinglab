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

    def test_debye_waller_is_always_the_more_pessimistic(self):
        """Which is the reason both are offered. Debye-Waller knows only about
        the incident wave, so it damps the same whether the field penetrates or
        not; Névot-Croce carries the transmitted component and damps less where
        the transmitted wave is evanescent.

        The direction is asserted, not a magnitude: at sigma = 0.5 nm and
        lambda = 2 nm the gap runs 0.1% at half a degree to 3% at 4 degrees --
        real, one-signed, and small. Pinning a ratio would be pinning a guess.
        """
        for graze_deg in (0.5, 1.0, 2.0, 4.0):
            graze = np.radians(graze_deg)
            common = dict(roughness_nm=0.5, wavelength_nm=self.WAVELENGTH)
            nc = float(reflectivity(1.0 - 1e-3 + 1e-4j, graze,
                                    model="nevot-croce", **common))
            dw = float(reflectivity(1.0 - 1e-3 + 1e-4j, graze,
                                    model="debye-waller", **common))
            assert dw < nc, f"{graze_deg} deg"

    def test_and_the_gap_widens_away_from_grazing(self):
        """The mechanism, not just the sign: the two agree where the wave
        barely penetrates and separate as it starts to."""
        common = dict(roughness_nm=0.5, wavelength_nm=self.WAVELENGTH)

        def gap(graze_deg):
            graze = np.radians(graze_deg)
            nc = float(reflectivity(1.0 - 1e-3 + 1e-4j, graze,
                                    model="nevot-croce", **common))
            dw = float(reflectivity(1.0 - 1e-3 + 1e-4j, graze,
                                    model="debye-waller", **common))
            return nc / dw

        gaps = [gap(d) for d in (0.5, 1.0, 2.0, 4.0)]
        assert gaps == sorted(gaps)
        assert gaps[0] < 1.005 < gaps[-1]

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


class TestNothingImportsThisYet:
    def test_the_solver_does_not_reference_reflectivity(self):
        """M15-B is the pure layer. Wiring is M15-D, and it changes efficiency
        values -- the first commit in this project to do so. Keeping the steps
        apart is what makes that diff reviewable."""
        import inspect
        import sys

        import gratinglab.solvers.scalar  # noqa: F401

        # The *module*, not the registered singleton `solvers.scalar`, which
        # shadows it on the package.
        module = sys.modules["gratinglab.solvers.scalar"]
        assert "reflectivity" not in inspect.getsource(module)
