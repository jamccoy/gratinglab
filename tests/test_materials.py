r"""Optical constants: the readers, the range guard, and the lookup.

The reference values here come from **outside this repo** -- published Henke
constants and the analytic definition of the critical angle -- because a reader
checked only against its own output verifies the parser and nothing else. That
is the same reasoning that put `check_reciprocity` in `checks.py` rather than
more closed-form tests in `test_scalar.py`.
"""

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
    path.write_text(CXRO_TEXT)
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
        path.write_text("h1\nh2\n1.0 2.0\n3.0 4.0\n")
        with pytest.raises(ValueError, match="three numbers"):
            read_cxro(path)

    def test_the_header_is_detected_rather_than_counted(self, tmp_path):
        """The prototype hardcodes `skip_header=3` against a two-header-line
        file and so drops the first data row -- the lowest energy, which is the
        *longest* wavelength. Both variants must read identically, and neither
        may lose a row."""
        two = tmp_path / "two.txt"
        two.write_text(CXRO_TEXT)
        three = tmp_path / "three.txt"
        three.write_text(" extra provenance line\n" + CXRO_TEXT)

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
        r""":math:`\theta_c = \sqrt{2\delta}`, and the decrement rises towards
        longer wavelengths, so grazing optics get more forgiving there. About
        3.1 degrees at 1 nm for Au."""
        assert np.degrees(gold.critical_angle(1.0)) == pytest.approx(3.1, abs=0.3)
        angles = np.degrees(gold.critical_angle([1.0, 2.0, 4.0]))
        assert np.all(np.diff(angles) > 0)

    def test_it_matches_the_square_root_definition_exactly(self, gold):
        """Not a re-implementation -- the point is that `critical_angle` reads
        the same table `at` does, so the two can never disagree."""
        decrement, _ = gold.at(2.5)
        assert gold.critical_angle(2.5) == pytest.approx(np.sqrt(2.0 * decrement))


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
        """`SOURCES.md` says the vendored data can be deleted outright if the
        redistribution terms do not hold up, and that this is a deletion rather
        than a rewrite. That claim is worth a test."""
        monkeypatch.setattr(materials, "_DATA", Path("/nonexistent"))
        assert materials.available() == ()


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
