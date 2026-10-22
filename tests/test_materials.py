r"""Optical constants: the readers, the range guard, and the lookup.

The reference values here come from **outside this repo** -- published Henke
constants and the analytic definition of the critical angle -- because a reader
checked only against its own output verifies the parser and nothing else. That
is the same reasoning that put `check_reciprocity` in `checks.py` rather than
more closed-form tests in `test_scalar.py`.
"""

import os
from pathlib import Path

import numpy as np
import pytest

from gratinglab import materials
from gratinglab.materials.optical import (
    HC_EV_NM,
    OpticalConstants,
    read_ari,
    read_cxro,
    write_ari,
)

#: A CXRO export, as downloaded: two header lines, ascending energy.
CXRO_TEXT = """\
 Xx Density=1.0
 Energy(eV), Delta, Beta
  1239.84193  0.002  0.0004
  2479.68386  0.001  0.0001
"""


@pytest.fixture
def cxro_file(tmp_path) -> Path:
    path = tmp_path / "Xx_CXRO.txt"
    path.write_text(CXRO_TEXT, encoding="utf-8")
    return path


@pytest.fixture(scope="module")
def gold() -> OpticalConstants:
    return materials.lookup("Au")


class TestTheCxroReader:
    def test_energy_becomes_wavelength(self, cxro_file):
        """1239.84193 eV is exactly 1 nm, by the definition of `HC_EV_NM`."""
        table = read_cxro(cxro_file)
        assert table.wavelength_nm == pytest.approx([0.5, 1.0])

    def test_and_the_table_comes_back_ascending_in_wavelength(self, cxro_file):
        """CXRO exports ascending *energy*, which is descending wavelength.
        `np.interp` requires ascending, and would return silent nonsense rather
        than raise if handed a reversed table."""
        table = read_cxro(cxro_file)
        assert np.all(np.diff(table.wavelength_nm) > 0)

    def test_the_rows_travel_with_their_wavelengths(self, cxro_file):
        """Non-vacuity for the reversal: sorting the wavelengths while leaving
        the columns alone would still ascend, and be wrong."""
        table = read_cxro(cxro_file)
        # 1 nm was the *second* row of the file (lower energy).
        assert table.at(1.0)[0] == pytest.approx(0.002)
        assert table.at(0.5)[0] == pytest.approx(0.001)

    def test_the_name_is_taken_from_the_filename(self, cxro_file):
        assert read_cxro(cxro_file).name == "Xx"

    def test_but_can_be_overridden(self, cxro_file):
        assert read_cxro(cxro_file, name="something else").name == "something else"

    def test_a_file_with_too_few_columns_is_refused(self, tmp_path):
        path = tmp_path / "bad.txt"
        path.write_text("h1\nh2\n1.0 2.0\n3.0 4.0\n", encoding="utf-8")
        with pytest.raises(ValueError, match="three numbers"):
            read_cxro(path)

    def test_the_header_is_detected_rather_than_counted(self, tmp_path):
        """The prototype hardcodes `skip_header=3` against a two-header-line
        file and so drops the first data row -- the lowest energy, which is the
        *longest* wavelength. Both variants must read identically, and neither
        may lose a row."""
        two = tmp_path / "two.txt"
        two.write_text(CXRO_TEXT, encoding="utf-8")
        three = tmp_path / "three.txt"
        three.write_text(" extra provenance line\n" + CXRO_TEXT, encoding="utf-8")

        a, b = read_cxro(two), read_cxro(three)
        assert len(a.wavelength_nm) == 2
        assert a.wavelength_nm == pytest.approx(b.wavelength_nm)

    def test_the_longest_wavelength_row_survives(self):
        """The row the prototype dropped, in the real vendored table. 200 eV is
        the first data line of the CXRO export and the end of the range a
        grazing-incidence scan is most likely to want."""
        gold = materials.lookup("Au")
        assert gold.range_nm[1] == pytest.approx(HC_EV_NM / 200.0, rel=1e-9)

    def test_which_the_corpus_ari_files_do_not_have(self, tmp_path):
        """Non-vacuity, and the size of the defect: the .ari tables in the
        corpus were written by the prototype and carry the truncation, so the
        two ranges differ by a real amount rather than by rounding."""
        gold = materials.lookup("Au")
        truncated = HC_EV_NM / 200.925018   # what skip_header=3 leaves
        assert gold.range_nm[1] - truncated == pytest.approx(0.0286, abs=1e-3)


class TestTheAriReader:
    def test_it_round_trips(self, tmp_path, cxro_file):
        original = read_cxro(cxro_file)
        path = write_ari(original, tmp_path / "Xx.ari")
        again = read_ari(path)
        # `%f` is six decimal places; that is the format's precision, not this
        # writer's, and it is why the tolerance is absolute rather than relative.
        assert again.wavelength_nm == pytest.approx(original.wavelength_nm, abs=5e-7)
        assert again.decrement == pytest.approx(original.decrement, abs=5e-7)

    def test_ari_columns_are_not_n_and_k(self, tmp_path, cxro_file):
        """The corpus files are named `*_optical_constants_n_k.ari` and hold
        `(wavelength, decrement, absorption)`. Reading column 1 as `n` would put
        the index at ~2e-3 instead of ~1 -- three orders of magnitude, so this
        cannot pass by accident."""
        path = write_ari(read_cxro(cxro_file), tmp_path / "Xx.ari")
        table = read_ari(path)
        assert table.decrement.max() < 0.01
        assert np.all(np.abs(table.n(1.0)) > 0.9)


class TestTheRangeIsEnforced:
    def test_asking_outside_it_raises(self, gold):
        with pytest.raises(ValueError, match="tabulated over"):
            gold.at(50.0)

    def test_the_message_names_the_range_and_what_was_asked(self, gold):
        with pytest.raises(ValueError) as caught:
            gold.at(50.0)
        assert "50" in str(caught.value)
        assert f"{gold.range_nm[1]:.4g}" in str(caught.value)

    def test_one_bad_wavelength_in_an_array_is_enough(self, gold):
        """The alternative -- returning good values for the good entries -- puts
        a silent hole in a scan."""
        with pytest.raises(ValueError):
            gold.at([1.0, 2.0, 50.0])

    def test_inside_the_range_is_fine(self, gold):
        """Non-vacuity: the guard has to let real work through."""
        low, high = gold.range_nm
        assert gold.at([low, 0.5 * (low + high), high])[0].shape == (3,)

    def test_the_endpoints_are_included(self, gold):
        """A table covers its endpoints. `<` rather than `<=` here would make
        the shortest tabulated wavelength unusable -- the same off-by-one
        mutation testing found in `geometry.beta`."""
        low, high = gold.range_nm
        assert gold.covers([low, high]).all()
        assert not gold.covers([low - 1e-9, high + 1e-9]).any()


class TestAgainstPublishedValues:
    """Checked outside this repo, so the reader is not merely self-consistent."""

    def test_gold_at_1_nm(self, gold):
        """Au at 1239.84 eV: Henke gives decrement ~1.46e-3, absorption
        ~5.1e-4. Loose tolerances on purpose -- this is checking that the
        columns are the right way round and the units are right, not
        re-deriving the table."""
        decrement, absorption = gold.at(1.0)
        assert decrement == pytest.approx(1.46e-3, rel=0.05)
        assert absorption == pytest.approx(5.1e-4, rel=0.10)

    def test_the_index_is_just_below_one(self, gold):
        r"""The defining property of X-ray optics: :math:`n = 1 - \delta + i\beta`
        with a *positive* decrement, so the real part is below 1 and total
        external reflection exists at all. A sign slip here would put it above
        1 and quietly delete the critical angle."""
        n = gold.n([1.0, 2.0, 4.0])
        assert np.all(n.real < 1.0)
        assert np.all(n.real > 0.99)
        assert np.all(n.imag > 0.0)

    def test_the_critical_angle_is_a_few_degrees_and_grows_with_wavelength(self, gold):
        r""":math:`\theta_c \approx \sqrt{2\delta}`, and the decrement rises towards
        longer wavelengths, so grazing optics get more forgiving there. About
        3.1 degrees at 1 nm for Au."""
        assert np.degrees(gold.critical_angle(1.0)) == pytest.approx(3.1, abs=0.3)
        angles = np.degrees(gold.critical_angle([1.0, 2.0, 4.0]))
        assert np.all(np.diff(angles) > 0)

    def test_it_matches_the_exact_lossless_definition_exactly(self, gold):
        r"""Not a re-implementation -- the point is that `critical_angle` reads
        the same table `at` does, so the two can never disagree.

        This pinned :math:`\sqrt{2\delta}` until the materials layer met
        refractiveindex.info, and was changed on purpose rather than loosened:
        the arccos is what Snell's law gives and what `fresnel.amplitude`
        reflects totally below, so it is the definition, and the square root
        is now the expansion the next test holds it to."""
        decrement, _ = gold.at(2.5)
        assert gold.critical_angle(2.5) == pytest.approx(
            np.arccos(1.0 - decrement), rel=1e-12
        )

    def test_and_differs_from_the_square_root_only_by_its_next_term(self, gold):
        r"""What the change from :math:`\sqrt{2\delta}` costs X-ray results:
        :math:`\arccos(1-\delta) = \sqrt{2\delta}\,(1 + \delta/12 + \dots)`, so
        across the whole Au table the two agree to that term, about 1e-4.
        Asserting the term itself, not a tolerance wide enough to hide a
        wrong formula."""
        lam = np.geomspace(*gold.range_nm, 400)
        decrement, _ = gold.at(lam)
        relative = gold.critical_angle(lam) / np.sqrt(2.0 * decrement) - 1.0
        assert relative == pytest.approx(decrement / 12.0, rel=0.01)
        assert relative.max() < 2e-3


class TestLookup:
    def test_the_vendored_table_is_there(self):
        assert "Au" in materials.available()

    def test_and_is_a_real_table(self, gold):
        assert len(gold.wavelength_nm) > 100
        assert "CXRO" in gold.source

    def test_an_unknown_name_raises_and_says_what_is_available(self):
        """`Problem.coating` used to be a free-form string nothing read, so
        setting it to anything relabelled a result as absolute while changing
        no number. Resolution is what turns the string into a claim."""
        with pytest.raises(materials.UnknownMaterial) as caught:
            materials.lookup("unobtanium")
        assert "Au" in str(caught.value)
        assert "from_cxro_file" in str(caught.value)

    def test_it_is_cached_by_name(self):
        """Same object, not merely an equal one -- a solve may resolve the same
        coating repeatedly and reparsing 500 rows each time is waste."""
        assert materials.lookup("Au") is materials.lookup("Au")

    def test_an_empty_data_directory_is_a_supported_state(self, monkeypatch):
        """An install ships no tables (see `SOURCES.md`), so an empty library
        is what every user starts with, not a broken install."""
        monkeypatch.setattr(materials, "_DATA", Path("/nonexistent"))
        assert materials.available() == ()


class TestUserLibraries:
    """A library is a directory of tables, one material per file. An install
    ships none, so this is how anybody gets a coating at all."""

    @pytest.fixture
    def library(self, tmp_path, monkeypatch):
        monkeypatch.setenv(materials.LIBRARY_ENV, str(tmp_path))
        return tmp_path

    def test_a_table_dropped_in_is_a_coating(self, library):
        (library / "Edge.txt").write_text(EDGE_TEXT, encoding="utf-8")
        assert "Edge" in materials.available()
        assert materials.lookup("Edge").name == "Edge"

    def test_the_source_names_the_file_it_came_from(self, library):
        """Provenance has to say which table answered, now that there can be
        more than one candidate for a name."""
        (library / "Edge.txt").write_text(EDGE_TEXT, encoding="utf-8")
        source = materials.lookup("Edge").source
        assert "Edge.txt" in source and str(library) in source

    def test_yours_replaces_a_bundled_table_of_the_same_name(self, library):
        (library / "Au.txt").write_text(EDGE_TEXT, encoding="utf-8")
        assert len(materials.lookup("Au").wavelength_nm) == 6
        assert materials.available().count("Au") == 1

    def test_several_directories_are_searched_in_order(self, tmp_path, monkeypatch):
        first, second = tmp_path / "first", tmp_path / "second"
        first.mkdir()
        second.mkdir()
        (second / "Edge.txt").write_text(EDGE_TEXT, encoding="utf-8")
        (second / "Only.txt").write_text(EDGE_TEXT, encoding="utf-8")
        (first / "Edge.txt").write_text(CXRO_TEXT, encoding="utf-8")
        monkeypatch.setenv(materials.LIBRARY_ENV, f"{first}{os.pathsep}{second}")
        assert len(materials.lookup("Edge").wavelength_nm) != 6   # first wins
        assert materials.lookup("Only").name == "Only"           # second searched

    def test_a_pcgrate_ari_table_is_read_too(self, library, gold):
        materials.write_ari(gold, library / "Gold2.ari")
        table = materials.lookup("Gold2")
        assert table.decrement == pytest.approx(gold.decrement, abs=1e-6)

    def test_two_tables_for_one_name_in_one_directory_are_refused(self, library, gold):
        """Neither suffix outranks the other, so picking one would be a guess."""
        (library / "Edge.txt").write_text(EDGE_TEXT, encoding="utf-8")
        materials.write_ari(gold, library / "Edge.ari")
        with pytest.raises(ValueError, match="two tables for 'Edge'"):
            materials.lookup("Edge")

    def test_an_edited_table_is_reread(self, library):
        """The cache is keyed on the file and its mtime, not the name, or an
        edit would keep serving the old numbers for the life of the process."""
        path = library / "Edge.txt"
        path.write_text(EDGE_TEXT, encoding="utf-8")
        before = materials.lookup("Edge")
        path.write_text(CXRO_TEXT, encoding="utf-8")
        os.utime(path, ns=(path.stat().st_atime_ns, path.stat().st_mtime_ns + 10**9))
        assert materials.lookup("Edge") is not before

    def test_without_the_variable_the_home_library_is_used(self, tmp_path, monkeypatch):
        monkeypatch.delenv(materials.LIBRARY_ENV)
        monkeypatch.setattr(materials, "user_dir", lambda: tmp_path)
        (tmp_path / "Edge.txt").write_text(EDGE_TEXT, encoding="utf-8")
        assert materials.lookup("Edge").name == "Edge"

    def test_a_refusal_says_where_it_looked_and_how_to_add_one(self, library):
        with pytest.raises(materials.UnknownMaterial) as caught:
            materials.lookup("Pt")
        message = str(caught.value)
        assert str(library) in message
        assert "Pt.txt" in message and materials.LIBRARY_ENV in message


class TestEveryVendoredTable:
    """Each file in `data/` is a claim made in `SOURCES.md`, so each is held to
    it rather than Au alone standing in for any that arrive later."""

    @pytest.mark.parametrize("name", materials.available())
    def test_it_loads_as_a_cxro_table(self, name):
        table = materials.lookup(name)
        assert "CXRO" in table.source
        assert len(table.wavelength_nm) == 500

    @pytest.mark.parametrize("name", materials.available())
    def test_it_spans_200_ev_to_2_kev(self, name):
        low, high = materials.lookup(name).range_nm
        assert low == pytest.approx(0.6199, abs=1e-3)
        assert high == pytest.approx(6.199, abs=1e-3)

    @pytest.mark.parametrize("name", materials.available())
    def test_its_index_is_close_to_one_and_lossy(self, name):
        """Close to one, not below it everywhere: just under an absorption edge
        the decrement goes negative, so the bound is two-sided."""
        table = materials.lookup(name)
        assert (np.abs(table.decrement) < 0.05).all()
        assert (table.absorption > 0).all()

    @pytest.mark.parametrize("name", materials.available())
    def test_it_has_a_row_in_sources(self, name):
        sources = (materials.data_dir() / "SOURCES.md").read_text(encoding="utf-8")
        assert f"| `{name}.txt` |" in sources


#: A synthetic CXRO export with a negative decrement at 572.7 eV, the shape
#: real data takes just below an L3 edge (chromium's is 574 eV). Synthetic so
#: the test needs no table the package may not ship.
EDGE_TEXT = """\
 Edge Density=1.0
 Energy(eV), Delta, Beta
  200.  0.0100  0.0010
  400.  0.0050  0.0010
  572.7  -0.0050  0.0010
  600.  0.0040  0.0010
  1000.  0.0020  0.0010
  2000.  0.0005  0.0010
"""
EDGE_NM = 1239.84 / 572.7


class TestBelowAnAbsorptionEdge:
    @pytest.fixture
    def edge(self, tmp_path):
        path = tmp_path / "Edge.txt"
        path.write_text(EDGE_TEXT, encoding="utf-8")
        return materials.from_cxro_file(path)

    def test_the_decrement_really_is_negative_there(self, edge):
        """Non-vacuity for the test below: without this the clip is untested."""
        decrement, _ = edge.at(EDGE_NM)
        assert decrement < 0

    def test_the_critical_angle_is_zero_not_nan(self, edge):
        """n above 1 deletes total external reflection, so no graze is below
        critical. NaN said nothing, and every comparison against it is False."""
        assert edge.critical_angle(EDGE_NM) == 0.0

    def test_and_it_is_still_in_the_regime_where_that_matters(self, edge):
        """n is above 1 but only just, so reflectivity has nothing to stand on
        and the scalar solver must still warn here. Being outside the regime
        is what silences the warning for a glass, and this is not a glass."""
        assert edge.near_unity(EDGE_NM)


RII = Path(__file__).parent / "data" / "refractiveindex"


@pytest.fixture(scope="module")
def au():
    return materials.from_refractiveindex_file(RII / "tabulated-nk_Au-Johnson.yml")


@pytest.fixture(scope="module")
def bk7():
    return materials.from_refractiveindex_file(RII / "formula-2-tabulated-k_N-BK7.yml")


class TestOutsideTheXRay:
    """From the UV to the IR there is usually no critical angle at all.

    The committed refractiveindex.info fixtures are the two shapes that went
    wrong when the small-angle form met them: a visible metal, whose real part
    is below 1 but whose absorption is of order 1, and a glass, whose real part
    is above 1. Both now read as having no critical angle, and both as outside
    the regime where that is worth a warning.
    """

    def test_visible_gold_has_its_real_part_below_one(self, au):
        """Non-vacuity: this is the case the square root answered with 72
        degrees and the arccos would answer with 76. Neither means anything."""
        n = complex(au.n(600.0))
        assert 0.0 < n.real < 1.0
        assert n.imag > 1.0

    def test_and_no_critical_angle(self, au):
        assert au.critical_angle(600.0) == 0.0
        assert not au.near_unity(600.0)

    def test_because_fresnel_finds_no_edge_there(self, au):
        """The physics behind the zero, checked on the exact amplitudes rather
        than asserted: from 10 degrees to normal incidence the reflectivity of
        Au at 600 nm stays between 0.90 and 0.93, with nothing at 76 degrees or
        anywhere else. A critical angle is the place where it would collapse,
        and an X-ray material loses more than 99% within two of them."""
        from gratinglab.materials.fresnel import reflectivity

        graze = np.radians(np.linspace(10.0, 90.0, 81))
        R = reflectivity(complex(au.n(600.0)), graze)
        assert R.min() > 0.90
        assert R.max() < 0.93

    def test_glass_has_no_critical_angle_either(self, bk7):
        lam = np.array([400.0, 550.0, 1000.0])
        assert np.all(bk7.n(lam).real > 1.0)
        assert np.all(bk7.critical_angle(lam) == 0.0)
        assert not bk7.near_unity(lam).any()

    def test_and_the_x_ray_table_is_inside_the_regime_everywhere(self):
        """The other side: nothing about Au's own X-ray constants moved out."""
        gold = materials.lookup("Au")
        assert gold.near_unity(np.geomspace(*gold.range_nm, 400)).all()


class TestConstruction:
    def test_a_descending_table_is_refused(self):
        with pytest.raises(ValueError, match="ascending"):
            OpticalConstants(
                name="Xx",
                wavelength_nm=np.array([2.0, 1.0]),
                decrement=np.array([1e-3, 2e-3]),
                absorption=np.array([1e-4, 2e-4]),
            )

    def test_mismatched_column_lengths_are_refused(self):
        with pytest.raises(ValueError, match="match"):
            OpticalConstants(
                name="Xx",
                wavelength_nm=np.array([1.0, 2.0]),
                decrement=np.array([1e-3]),
                absorption=np.array([1e-4, 2e-4]),
            )

    def test_a_single_point_cannot_be_interpolated_between(self):
        with pytest.raises(ValueError, match="at least two"):
            OpticalConstants(
                name="Xx",
                wavelength_nm=np.array([1.0]),
                decrement=np.array([1e-3]),
                absorption=np.array([1e-4]),
            )

    def test_the_arrays_are_frozen(self, gold):
        """Same rule as `EfficiencyScan`: a table that can be edited after
        construction is not a record of anything, and `lookup` hands the same
        cached object to every caller."""
        with pytest.raises(ValueError):
            gold.decrement[0] = 999.0


def test_hc_matches_the_prototype():
    """Kept identical to `panter1.py`'s `const`, so re-converting a corpus
    table reproduces it rather than shifting every wavelength in the fifth
    decimal place."""
    assert HC_EV_NM == 1239.84193
