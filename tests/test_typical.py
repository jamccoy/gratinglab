"""The typical-material model against its own physics, not against a table.

No real optical constants are committed, so nothing here compares with CXRO.
What is checked instead is what the model has to satisfy whatever the facts:
the Kramers-Kronig relation its closed form claims to solve, the sum rule its
normalisation claims to meet, the limits at high and low energy, linearity in
density, additivity over atoms, passivity, and that its range and its caveats
reach a result. The measured accuracy against real tables is in
``docs/typical-materials.md``.
"""

import numpy as np
import pytest
from scipy.integrate import quad

from gratinglab import materials
from gratinglab.illumination import Illumination
from gratinglab.materials.optical import HC_EV_NM
from gratinglab.materials.typical import (
    DENSITIES,
    ELEMENTS,
    HIGH_EV,
    LOW_EV,
    scattering_factors,
)
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed
from gratinglab.solvers import get_solver, scalar


def z_star(el):
    return el.z - (el.z / 82.5) ** 2.37


def f2_at(el, energy):
    return float(scattering_factors(el, np.array([energy]))[1][0])


def edges_of(el):
    return sorted(binding for _, binding, _ in el.absorbers)


class TestTheClosedFormIsTheKramersKronigTransform:
    """f1 is written in closed form from f2. Recompute it by numerical
    principal-value integration of f2 instead, which shares nothing with the
    closed form but the definition of f2, so a sign or a factor in the algebra
    shows up here."""

    @staticmethod
    def kk_f1(el, energy):
        # f1(E) = Z* + (2/pi) P int eps f2(eps) / (E^2 - eps^2) d eps
        #       = Z* - (2/pi) P int [eps f2 / (eps + E)] / (eps - E) d eps
        # Only the piece around E needs the Cauchy weight; the rest is smooth
        # and is integrated in ln(eps), which the E^-2 tail makes short.
        points = sorted({*edges_of(el), 0.5 * energy, 2.0 * energy, 1e12})
        points = [p for p in points if p >= edges_of(el)[0]]
        total = 0.0
        for a, b in zip(points, points[1:]):
            g = lambda eps: -eps * f2_at(el, eps) / (eps + energy)  # noqa: E731
            if a < energy < b:
                total += quad(g, a, b, weight="cauchy", wvar=energy, limit=400)[0]
            else:
                total += quad(
                    lambda u: g(np.exp(u)) / (np.exp(u) - energy) * np.exp(u),
                    np.log(a), np.log(b), limit=400,
                )[0]
        return z_star(el) + 2.0 / np.pi * total

    @pytest.mark.parametrize("symbol, energy", [
        ("C", 150.0), ("C", 1000.0), ("Si", 500.0), ("Si", 5000.0), ("Au", 1500.0),
    ])
    def test_at(self, symbol, energy):
        el = ELEMENTS[symbol]
        closed = float(scattering_factors(el, np.array([energy]))[0][0])
        assert self.kk_f1(el, energy) == pytest.approx(closed, rel=1e-6, abs=1e-6)


class TestTheSumRules:
    @pytest.mark.parametrize("symbol", list(ELEMENTS))
    def test_every_electron_absorbs_once(self, symbol):
        """(2/pi) int f2 dE/E = Z: each shell's strength is normalised to its
        own occupancy, and the valence electrons are a shell too."""
        el = ELEMENTS[symbol]
        points = [*edges_of(el), np.inf]
        total = sum(
            quad(lambda e: f2_at(el, e) / e, a, b, limit=200)[0]
            for a, b in zip(points, points[1:])
        )
        assert 2.0 / np.pi * total == pytest.approx(el.z, rel=1e-6)

    @pytest.mark.parametrize("symbol", list(ELEMENTS))
    def test_far_above_every_edge_f1_is_z_star(self, symbol):
        el = ELEMENTS[symbol]
        f1, f2 = scattering_factors(el, np.array([1e9]))
        assert f1[0] == pytest.approx(z_star(el), abs=1e-6)
        assert f2[0] < 1e-6

    @pytest.mark.parametrize("symbol", list(ELEMENTS))
    def test_far_below_them_the_bound_electrons_drop_out(self, symbol):
        el = ELEMENTS[symbol]
        f1, f2 = scattering_factors(el, np.array([1e-3]))
        assert f1[0] == pytest.approx(z_star(el) - el.z, abs=1e-4)
        assert f2[0] == 0.0


class TestScaling:
    def test_the_index_is_linear_in_density(self):
        thin, thick = materials.typical("C", 2.0), materials.typical("C", 4.0)
        for a, b in zip(thin.at([0.1, 1.0, 10.0]), thick.at([0.1, 1.0, 10.0])):
            np.testing.assert_allclose(b, 2.0 * a, rtol=1e-12)

    def test_well_above_its_edges_the_decrement_goes_as_lambda_squared(self):
        """The free-electron limit, and so a critical angle proportional to
        lambda sqrt(rho): carbon between 20 and 30 keV is ~70 times its K edge."""
        c = materials.typical("C")
        lam = HC_EV_NM / np.array([20000.0, 30000.0])
        d, _ = c.at(lam)
        assert d[0] / lam[0] ** 2 == pytest.approx(d[1] / lam[1] ** 2, rel=2e-3)

    def test_a_compound_is_the_sum_of_its_atoms(self):
        rho = DENSITIES["SiO2"]
        mass = ELEMENTS["Si"].weight + 2 * ELEMENTS["O"].weight
        si = materials.typical("Si", rho * ELEMENTS["Si"].weight / mass)
        o = materials.typical("O", rho * 2 * ELEMENTS["O"].weight / mass)
        lam = np.geomspace(0.05, 40.0, 50)
        total = materials.typical("SiO2").at(lam)
        for whole, a, b in zip(total, si.at(lam), o.at(lam)):
            np.testing.assert_allclose(whole, a + b, rtol=1e-12)

    @pytest.mark.parametrize("formula", list(DENSITIES))
    def test_it_absorbs_everywhere_in_its_range(self, formula):
        """Never lossless: the rigorous solvers treat zero absorption as a
        transmitting dielectric, which no real material is in the X-ray."""
        m = materials.typical(formula)
        _, absorption = m.at(np.geomspace(*m.range_nm, 4000))
        assert (absorption > 0).all()


class TestTheRange:
    def test_it_runs_from_30_ev_to_30_kev_inclusive(self):
        m = materials.typical("Au")
        assert m.range_nm == (HC_EV_NM / HIGH_EV, HC_EV_NM / LOW_EV)
        m.at(list(m.range_nm))

    @pytest.mark.parametrize("wavelength", [0.03, 50.0, 500.0])
    def test_outside_it_is_refused_as_a_table_is(self, wavelength):
        with pytest.raises(ValueError, match="Extrapolating"):
            materials.typical("Au").at(wavelength)

    def test_the_sampled_arrays_are_the_model(self):
        m = materials.typical("Ni")
        d, b = m.at(m.wavelength_nm)
        np.testing.assert_array_equal(d, m.decrement)
        np.testing.assert_array_equal(b, m.absorption)


class TestNames:
    def test_lookup_by_prefix(self):
        au = materials.lookup("typical:Au")
        assert au.name == "typical:Au" and au.density == DENSITIES["Au"]
        assert "model" in au.source and "approximate" in au.source

    def test_a_density_after_an_at(self):
        c = materials.lookup("typical:C@2.05")
        assert c.density == 2.05 and c.name == "typical:C@2.05"

    def test_without_the_prefix_a_missing_table_still_refuses(self):
        """The model never stands in for a table silently."""
        with pytest.raises(materials.UnknownMaterial):
            materials.lookup("Pt")

    @pytest.mark.parametrize("name, message", [
        ("typical:Fe", "no facts for 'Fe'"),
        ("typical:O", "no default density"),
        ("typical:Au@heavy", "must be a number"),
        ("typical:Au@-1", "must be positive"),
        ("typical:sio2", "not a chemical formula"),
    ])
    def test_what_it_cannot_answer_is_refused(self, name, message):
        with pytest.raises(ValueError, match=message):
            materials.lookup(name)


class TestCaveats:
    def test_near_an_edge_it_says_which(self):
        (warning,) = materials.typical("C").caveats([4.4])
        assert "C K at 284.2 eV" in warning and "qualitative" in warning

    def test_between_edges_it_says_nothing(self):
        assert materials.typical("C").caveats([1.0, 2.0]) == ()

    def test_below_100_ev_it_says_so(self):
        (warning,) = materials.typical("C").caveats([13.5])
        assert "below 100 eV" in warning

    def test_gold_names_its_m5_edge(self):
        assert any("Au M5" in w for w in materials.typical("Au").caveats([0.56]))

    def test_above_gold_4f_it_names_the_delayed_onset(self):
        (warning,) = materials.typical("Au").caveats([HC_EV_NM / 200.0])
        assert "Au N7 at 84 eV" in warning and "delayed onset" in warning

    def test_but_not_above_gold_3d(self):
        """Measured: the M4,5 band shows no delayed-onset error, so it is not
        reported; 8 keV is under 4x the M5 energy."""
        assert materials.typical("Au").caveats([HC_EV_NM / 8000.0]) == ()

    def test_a_table_has_none(self):
        assert materials.lookup("Au").caveats([1.0, 4.4]) == ()


class TestTheyReachAResult:
    ill = Illumination.offplane(graze=1.5, azimuth=25.0, polarization="unpolarized")

    def problem(self):
        return Problem(period=160.0, profile=Blazed(blaze_angle=10.0), coating="typical:C")

    def test_scalar(self):
        near = scalar.solve(self.problem(), self.ill, [4.4])
        clear = scalar.solve(self.problem(), self.ill, [1.0])
        assert any("C K at 284.2 eV" in w for w in near.provenance.warnings)
        assert not any("typical-material" in w for w in clear.provenance.warnings)
        assert near.provenance.notes["coating"].startswith("typical:C (typical-material")

    def test_rcwa(self):
        scan = get_solver("rcwa").solve(self.problem(), self.ill, [4.4], truncation=5)
        assert any("C K at 284.2 eV" in w for w in scan.provenance.warnings)
