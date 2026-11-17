"""Tabulated scattering factors, and CXRO tables rescaled to a density.

No real factors are committed (``materials/data/SOURCES.md``), so every table
here is invented in CXRO's ``.nff`` layout and written into a throwaway
library. What is checked is the arithmetic a table cannot hide: the index of
one element against the relation written out by hand, additivity over atoms,
linearity in density, the range a compound inherits, and the refusals.
"""

import importlib

import numpy as np
import pytest

from gratinglab import materials
from gratinglab.gui.state import FormState, FormErrors, build
from gratinglab.materials.elements import ATOMIC_NUMBER, WEIGHTS
from gratinglab.materials.factors import (
    N_A,
    R_E,
    atoms_per_m3,
    parse_formula,
    read_nff,
)
from gratinglab.materials.optical import HC_EV_NM

typical_module = importlib.import_module("gratinglab.materials.typical")


def nff(path, energy, f1, f2, header="E(eV)\tf1\tf2"):
    rows = [f"{e:.6g}\t{a:.6g}\t{b:.6g}" for e, a, b in zip(energy, f1, f2)]
    path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def library(tmp_path, monkeypatch):
    """A library holding invented Si and O factors over 50 eV to 5 keV."""
    monkeypatch.setenv("GRATINGLAB_MATERIALS", str(tmp_path))
    e = np.geomspace(50.0, 5000.0, 200)
    nff(tmp_path / "si.nff", e, 14.0 - 2.0 / np.sqrt(e / 50.0), 3.0 * (e / 50.0) ** -1.5)
    nff(tmp_path / "O.nff", np.geomspace(30.0, 4000.0, 150), np.full(150, 8.0), np.full(150, 0.5))
    return tmp_path


class TestTheIndexIsTheDissertationRelation:
    """delta + i beta = r_e lambda^2 n (f1 + i f2) / 2 pi (McCoy, Appendix C)."""

    # Checked at tabulated nodes: between them `at` interpolates the index,
    # which carries a lambda^2, so it differs from interpolating f1 by ~1e-4.

    def test_one_element_by_hand(self, library):
        got = materials.factors("O", density=1.5)
        wl = got.wavelength_nm[[0, 70, -1]]
        n_atoms = 1.5e6 / WEIGHTS["O"] * N_A
        scale = R_E * (wl * 1e-9) ** 2 * n_atoms / (2 * np.pi)
        d, b = got.at(wl)
        assert np.allclose(d, scale * 8.0, rtol=1e-12)
        assert np.allclose(b, scale * 0.5, rtol=1e-12)

    def test_a_compound_is_the_sum_of_its_atoms(self, tmp_path, monkeypatch):
        monkeypatch.setenv("GRATINGLAB_MATERIALS", str(tmp_path))
        e = np.geomspace(50.0, 5000.0, 120)
        nff(tmp_path / "Si.nff", e, 14.0 - 2.0 / np.sqrt(e / 50.0), 3.0 * (e / 50.0) ** -1.5)
        nff(tmp_path / "O.nff", e, 8.0 - np.log(e / 50.0), 0.5 * (e / 50.0) ** -2)
        sio2 = materials.factors("SiO2", density=2.2)
        si = materials.factors("Si", density=2.2)
        o = materials.factors("O", density=2.2)
        wl = sio2.wavelength_nm  # every node, shared by both elements
        m_si, m_o = WEIGHTS["Si"], WEIGHTS["O"]
        # At one mass density, rescale each element to SiO2's formula units.
        for part in (0, 1):
            expect = (si.at(wl)[part] * m_si + 2 * o.at(wl)[part] * m_o) / (m_si + 2 * m_o)
            assert np.allclose(sio2.at(wl)[part], expect, rtol=1e-12)

    def test_linear_in_density(self, library):
        wl = [1.0, 3.0]
        a = materials.factors("SiO2", density=2.2).at(wl)
        b = materials.factors("SiO2", density=1.1).at(wl)
        assert np.allclose(a[0], 2 * b[0]) and np.allclose(a[1], 2 * b[1])

    def test_the_typical_model_uses_the_same_relation(self):
        el = typical_module.ELEMENTS["Au"]
        wl = np.array([1.0, 4.0])
        f1, f2 = typical_module.scattering_factors(el, HC_EV_NM / wl)
        n = 19.32e6 / el.weight * N_A
        scale = R_E * (wl * 1e-9) ** 2 * n / (2 * np.pi)
        d, b = materials.typical("Au").at(wl)
        assert np.allclose(d, scale * f1, rtol=1e-12) and np.allclose(b, scale * f2, rtol=1e-12)


class TestNumberDensity:
    @pytest.mark.parametrize(
        "formula, density, table",
        # McCoy dissertation, Appendix C, the table of mass and number
        # densities. It rounds; these agree with it to within 1%.
        [("C5H8O2", 1.19, 1.08e29), ("SiO2", 2.20, 6.62e28), ("Ni", 8.90, 9.14e28), ("Au", 19.32, 5.90e28)],
    )
    def test_matches_the_dissertation_table(self, formula, density, table):
        assert atoms_per_m3(formula, density) == pytest.approx(table, rel=0.01)

    def test_weights_agree_with_the_typical_model(self):
        for symbol, el in typical_module.ELEMENTS.items():
            assert WEIGHTS[symbol] == pytest.approx(el.weight, rel=1e-4)
            assert ATOMIC_NUMBER[symbol] == el.z

    def test_hydrogen_to_uranium(self):
        assert ATOMIC_NUMBER["H"] == 1 and ATOMIC_NUMBER["U"] == 92 and len(WEIGHTS) == 92


class TestReadingTables:
    def test_rows_without_f1_are_skipped(self, tmp_path):
        path = nff(tmp_path / "au.nff", [10.0, 20.0, 30.0, 40.0], [-9999.0, -9999.0, 5.0, 6.0], [1, 2, 3, 4])
        table = read_nff(path)
        assert table.symbol == "Au"
        assert list(table.energy_ev) == [30.0, 40.0]

    def test_a_stem_that_is_not_an_element_is_refused(self, tmp_path):
        path = nff(tmp_path / "gold.nff", [10.0, 20.0], [1, 2], [1, 2])
        with pytest.raises(ValueError, match="not an element"):
            read_nff(path)

    def test_formulas(self):
        assert parse_formula("C5H8O2") == (("C", 5.0), ("H", 8.0), ("O", 2.0))
        assert parse_formula("Ni0.8Cr0.2") == (("Ni", 0.8), ("Cr", 0.2))
        for bad in ("sio2", "SiO2x", "Xx2", ""):
            with pytest.raises(ValueError):
                parse_formula(bad)


class TestLookupByName:
    def test_a_name_resolves_and_says_where_it_came_from(self, library):
        got = materials.lookup("factors:SiO2@2.2")
        assert got.name == "factors:SiO2"  # 2.2 is the default density
        assert "si.nff" in got.source and "O.nff" in got.source

    def test_a_compound_covers_only_what_every_element_covers(self, library):
        low, high = materials.lookup("factors:SiO2").range_nm
        assert low == pytest.approx(HC_EV_NM / 4000.0)
        assert high == pytest.approx(HC_EV_NM / 50.0)
        with pytest.raises(ValueError):
            materials.lookup("factors:SiO2").at([HC_EV_NM / 4500.0])

    def test_a_missing_element_names_itself_and_where_to_get_it(self, library):
        with pytest.raises(materials.UnknownMaterial, match="Ti.*henke.lbl.gov"):
            materials.lookup("factors:TiO2@4.2")

    def test_no_default_density_is_refused_by_name(self, library):
        with pytest.raises(ValueError, match="factors:SiO@"):
            materials.lookup("factors:SiO")

    def test_the_installed_elements_are_listed_by_atomic_number(self, library):
        assert materials.factor_elements() == ("O", "Si")  # from O.nff and si.nff

    def test_install_files_it_under_the_symbol(self, library, tmp_path_factory):
        src = nff(tmp_path_factory.mktemp("dl") / "ti.nff", [30.0, 5000.0], [20, 22], [3, 1])
        assert materials.install(src).name == "Ti.nff"
        assert materials.lookup("factors:TiO2@4.2").range_nm


class TestRescalingACxroTable:
    @pytest.fixture
    def cxro(self, tmp_path, monkeypatch):
        monkeypatch.setenv("GRATINGLAB_MATERIALS", str(tmp_path))
        rows = "\n".join(f" {e:g} {0.01 * 200 / e:.6g} {0.005 * 200 / e:.6g}" for e in np.geomspace(200, 2000, 20))
        (tmp_path / "Pd.txt").write_text(" Pd Density=12.\n Energy(eV), Delta, Beta\n" + rows + "\n", encoding="utf-8")
        (tmp_path / "Nodensity.txt").write_text(" Synthetic\n Energy(eV), Delta, Beta\n" + rows + "\n", encoding="utf-8")
        return tmp_path

    def test_scales_both_parts(self, cxro):
        base = materials.lookup("Pd").at([2.0])
        film = materials.lookup("Pd@9")
        assert film.name == "Pd@9"
        assert film.at([2.0])[0] == pytest.approx(base[0] * 0.75)
        assert film.at([2.0])[1] == pytest.approx(base[1] * 0.75)
        assert "rescaled from 12 to 9" in film.source

    def test_a_table_without_a_density_is_refused(self, cxro):
        with pytest.raises(ValueError, match="does not state the density"):
            materials.lookup("Nodensity@3")

    def test_an_unknown_base_is_unknown(self, cxro):
        with pytest.raises(materials.UnknownMaterial):
            materials.lookup("Rh@11")


class TestTheForm:
    def test_accepts_model_and_factor_names(self, library):
        assert build(FormState(coating="typical:Au")).problem.coating == "typical:Au"
        assert build(FormState(coating="factors:SiO2@2.2")).problem.coating == "factors:SiO2@2.2"

    def test_reports_why_a_name_does_not_resolve(self, library):
        with pytest.raises(FormErrors) as caught:
            build(FormState(coating="factors:TiO2@4.2"))
        assert any("Ti" in e.message for e in caught.value.errors)
