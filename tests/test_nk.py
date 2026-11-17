"""The plain ``wavelength, n[, k]`` reader, on tables written here.

What matters most is the unit: a table read in the wrong one is a plausible
material at the wrong wavelengths, so every way of naming it is checked, and
a header that names none is checked to be refused rather than guessed.
"""

import numpy as np
import pytest

from gratinglab import materials
from gratinglab.materials.nk import read_nk
from gratinglab.materials.optical import HC_EV_NM


def table(tmp_path, text, name="film.csv"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


class TestTheUnit:
    @pytest.mark.parametrize(
        "header, scale",
        [
            ("wavelength (nm), n, k", 1.0),
            ("Wavelength [nm]\tn\tk", 1.0),
            ("lambda_um n k", 1000.0),
            ("Wavelength (µm), n, k", 1000.0),
            ("wl/micron, n, k", 1000.0),
            ("Wavelength (A), n, k", 0.1),
            ("lambda_angstrom; n; k", 0.1),
        ],
    )
    def test_wavelength_units(self, tmp_path, header, scale):
        sep = ";" if ";" in header else ","
        rows = f"\n0.5{sep}1.5{sep}0.1\n0.7{sep}1.4{sep}0.2\n"
        got = read_nk(table(tmp_path, header + rows))
        assert list(got.wavelength_nm) == pytest.approx([0.5 * scale, 0.7 * scale])

    def test_energy_is_converted_and_sorted(self, tmp_path):
        got = read_nk(table(tmp_path, "Energy (eV), n, k\n100, 0.99, 0.01\n200, 0.995, 0.002\n"))
        assert list(got.wavelength_nm) == pytest.approx([HC_EV_NM / 200, HC_EV_NM / 100])
        assert got.absorption[0] == pytest.approx(0.002)

    @pytest.mark.parametrize("header", ["wavelength, n, k", "x y z", "nm or um, n, k"])
    def test_a_header_without_exactly_one_unit_is_refused(self, tmp_path, header):
        with pytest.raises(ValueError, match="unit"):
            read_nk(table(tmp_path, header + "\n500, 1.5, 0\n600, 1.4, 0\n"))

    def test_no_header_at_all_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="no header"):
            read_nk(table(tmp_path, "500, 1.5, 0\n600, 1.4, 0\n"))


class TestTheValues:
    def test_stored_as_decrement_and_absorption(self, tmp_path):
        got = read_nk(table(tmp_path, "# a comment\nwavelength (nm), n, k\n400, 1.5, 0.25\n500, 1.4, 0.5\n"))
        assert got.name == "film"
        assert np.allclose(got.decrement, [-0.5, -0.4]) and np.allclose(got.absorption, [0.25, 0.5])
        assert complex(got.n([450.0])[0]) == pytest.approx(complex(1.45, 0.375))

    def test_two_columns_means_k_not_given(self, tmp_path):
        got = read_nk(table(tmp_path, "wavelength (nm), n\n400, 1.5\n500, 1.4\n"))
        assert not got.absorption.any()
        assert "not measured" in got.source

    def test_negative_k_is_clipped_and_recorded(self, tmp_path):
        got = read_nk(table(tmp_path, "wavelength (nm), n, k\n400, 1.5, -0.001\n500, 1.4, 0.2\n"))
        assert got.absorption[0] == 0.0
        assert "1 negative k clipped" in got.source

    def test_a_repeated_wavelength_names_both_lines(self, tmp_path):
        with pytest.raises(ValueError, match="lines 2 and 4"):
            read_nk(table(tmp_path, "wavelength (nm), n, k\n400, 1.5, 0\n500, 1.4, 0\n400, 1.6, 0\n"))

    def test_ragged_rows_and_trailing_text_are_refused(self, tmp_path):
        with pytest.raises(ValueError, match="columns"):
            read_nk(table(tmp_path, "wavelength (nm), n, k\n400, 1.5, 0\n500, 1.4\n"))
        with pytest.raises(ValueError, match="line 4"):
            read_nk(table(tmp_path, "wavelength (nm), n, k\n400, 1.5, 0\n500, 1.4, 0\nend\n"))

    def test_a_spreadsheet_bom_is_ignored(self, tmp_path):
        path = tmp_path / "film.csv"
        path.write_bytes("﻿wavelength (nm),n,k\r\n400,1.5,0\r\n500,1.4,0\r\n".encode("utf-8"))
        assert read_nk(path).range_nm == (400.0, 500.0)


class TestInALibrary:
    def test_nk_and_csv_are_coatings(self, tmp_path, monkeypatch):
        monkeypatch.setenv(materials.LIBRARY_ENV, str(tmp_path))
        table(tmp_path, "wavelength (nm), n, k\n400, 1.5, 0\n500, 1.4, 0\n", "Glass.nk")
        table(tmp_path, "wavelength (um), n, k\n0.4, 1.6, 0.1\n0.5, 1.5, 0.1\n", "Film.csv")
        assert {"Glass", "Film"} <= set(materials.available())
        assert materials.lookup("Film").range_nm == (400.0, 500.0)
        assert "Glass.nk" in materials.lookup("Glass").source

    def test_install_checks_it_reads_first(self, tmp_path, monkeypatch):
        monkeypatch.setenv(materials.LIBRARY_ENV, str(tmp_path / "lib"))
        bad = table(tmp_path, "wavelength, n, k\n400, 1.5, 0\n500, 1.4, 0\n", "Bad.csv")
        with pytest.raises(ValueError, match="unit"):
            materials.install(bad)
        assert not (tmp_path / "lib" / "Bad.csv").exists()

    def test_from_nk_file(self, tmp_path):
        path = table(tmp_path, "wavelength (nm), n, k\n400, 1.5, 0\n500, 1.4, 0\n")
        assert materials.from_nk_file(path, name="Mine").name == "Mine"
