"""refractiveindex.info YAML: the reader against the database's own numbers.

The fixtures are verbatim database files and `reference.json` was computed by
the database's own tooling, not by this package (see the README beside them).
Two glasses are also held to their published d-line indices, which neither
implementation produced.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from gratinglab import materials
from gratinglab.materials import refractiveindex as rii
from gratinglab.materials.fresnel import reflectivity

# Skipped rather than failed without the extra, like the Qt tests; CI installs
# it by name and checks it resolved, so a skip there cannot pass for a run.
yaml = pytest.importorskip("yaml")

FIXTURES = Path(__file__).parent / "data" / "refractiveindex"
REFERENCE = json.loads((FIXTURES / "reference.json").read_text(encoding="utf-8"))

#: The helium d line, where catalogue glasses quote n_d.
D_LINE_UM = 0.5875618


def read(name: str):
    return materials.from_refractiveindex_file(FIXTURES / name)


def write(tmp_path, text: str, name: str = "x.yml") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


class TestAgainstTheDatabase:
    @pytest.mark.parametrize("name", sorted(REFERENCE))
    def test_n_matches_the_database_tool(self, name):
        """Formula entries are sampled, so they are held to the sampling
        tolerance; a tabulated row is held to exact recovery."""
        ref = REFERENCE[name]
        table = read(name)
        nm = np.array(ref["wavelength_um"]) * 1000
        decrement, absorption = table.at(nm)
        tolerance = 2 * rii.FORMULA_TOLERANCE if ref["type"].startswith("formula") else 1e-12
        np.testing.assert_allclose(1 - decrement, ref["n"], rtol=0, atol=tolerance)
        if "k" in ref:
            np.testing.assert_allclose(absorption, ref["k"], rtol=0, atol=1e-12)

    def test_every_formula_number_is_covered(self):
        covered = {REFERENCE[n]["type"] for n in REFERENCE if REFERENCE[n]["type"].startswith("formula")}
        assert covered == {f"formula {i}" for i in range(1, 10)}


class TestPublishedIndices:
    """Neither the database's code nor this package's wrote these numbers."""

    def test_fused_silica_at_the_d_line(self):
        """Malitson (1965), the reference fused silica: n_d = 1.4585."""
        decrement, _ = read("formula-1_SiO2-Malitson.yml").at(D_LINE_UM * 1000)
        assert 1 - decrement == pytest.approx(1.4585, abs=5e-5)

    def test_n_bk7_at_the_d_line(self):
        """Schott's catalogue value for N-BK7: n_d = 1.5168."""
        decrement, _ = read("formula-2-tabulated-k_N-BK7.yml").at(D_LINE_UM * 1000)
        assert 1 - decrement == pytest.approx(1.5168, abs=5e-5)

    def test_and_through_fresnel_it_reflects_four_percent(self):
        """At normal incidence R = ((n - 1)/(n + 1))^2 exactly, so this ties the
        reader's sign and units to the reflectivity a solver will compute."""
        n = read("formula-2-tabulated-k_N-BK7.yml").n(D_LINE_UM * 1000)
        r = reflectivity(n, np.pi / 2)
        assert r == pytest.approx(((1.5168 - 1) / (1.5168 + 1)) ** 2, abs=1e-4)


class TestUnitsAndSigns:
    def test_wavelengths_arrive_in_nanometres(self):
        low, high = read("tabulated-nk_Au-Johnson.yml").range_nm
        assert 180 < low < 200 and 1900 < high < 2000

    def test_a_glass_has_a_negative_decrement(self):
        """n above 1 is decrement below 0: the X-ray tables' own quantity, far
        from zero. Nothing is flipped to make it look familiar."""
        assert (read("formula-1_SiO2-Malitson.yml").decrement < 0).all()

    def test_gold_is_lossy_with_positive_k(self):
        assert (read("tabulated-nk_Au-Johnson.yml").absorption > 0).all()


class TestWhatTheSourceRecords:
    def test_n_without_k_says_k_was_taken_as_zero(self):
        table = read("formula-1_SiO2-Malitson.yml")
        assert (table.absorption == 0).all()
        assert "k not given, taken as 0" in table.source

    def test_a_formula_with_tabulated_k_combines_both(self):
        table = read("formula-2-tabulated-k_N-BK7.yml")
        assert (table.absorption > 0).any()
        assert "formula 2" in table.source and "tabulated k" in table.source

    def test_the_reference_is_carried_as_plain_text(self):
        source = read("tabulated-nk_Au-Johnson.yml").source
        assert "refractiveindex.info" in source and "<" not in source


class TestFormulaSampling:
    def test_interpolation_stays_within_tolerance(self):
        """Measured on a grid the sampler never saw. Ar's gas formula is the
        kind that has a resonance near the short end of its range."""
        name = "formula-6_Ar-Peck-15C.yml"
        table = read(name)
        low, high = (v / 1000 for v in table.range_nm)
        check = np.geomspace(low, high, 20001)
        entry = yaml.safe_load((FIXTURES / name).read_text(encoding="utf-8"))["DATA"][0]
        c = [float(v) for v in entry["coefficients"].split()]
        exact = rii.dispersion("6", c, check)
        assert np.abs(exact - (1 - table.at(check * 1000)[0])).max() <= 1.05 * rii.FORMULA_TOLERANCE

    def test_a_range_containing_its_own_pole_is_refused(self, tmp_path):
        """Xe/Bideau-Mehu in the database states 0.1404 to 0.6234 um with a
        pole at 0.147. n diverges there, and no table can follow it."""
        path = write(tmp_path, """
REFERENCES: test
DATA:
  - type: formula 6
    wavelength_range: 0.1404 0.6234
    coefficients: 0 0.00322869 46.301 0.00355393 59.578 0.0606764 112.74
""")
        with pytest.raises(ValueError, match="resonance pole"):
            materials.from_refractiveindex_file(path)

    def test_a_negative_n_squared_is_refused(self, tmp_path):
        path = write(tmp_path, """
REFERENCES: test
DATA:
  - type: formula 2
    wavelength_range: 0.5 2.0
    coefficients: -3 0 0
""")
        with pytest.raises(ValueError, match="n\\^2 < 0"):
            materials.from_refractiveindex_file(path)


class TestTabulatedQuirks:
    def test_rows_are_sorted_and_a_repeated_wavelength_averaged(self, tmp_path):
        """Seventeen blocks in the database are unsorted and 33 repeat a
        wavelength, where two measurement ranges were joined."""
        path = write(tmp_path, """
REFERENCES: test
DATA:
  - type: tabulated nk
    data: |
        0.6 1.50 0.10
        0.4 1.60 0.30
        0.5 1.54 0.20
        0.5 1.56 0.20
""")
        table = materials.from_refractiveindex_file(path)
        np.testing.assert_allclose(table.wavelength_nm, [400, 500, 600])
        assert 1 - table.at(500.0)[0] == pytest.approx(1.55)
        assert "averaged 1 repeated wavelength" in table.source

    def test_a_negative_k_is_clipped_and_recorded(self, tmp_path):
        path = write(tmp_path, """
REFERENCES: test
DATA:
  - type: tabulated nk
    data: |
        0.4 1.5 0.01
        0.5 1.5 -0.02
        0.6 1.5 0.01
""")
        table = materials.from_refractiveindex_file(path)
        assert (table.absorption >= 0).all()
        assert "clipped 1 negative k to 0, the lowest -0.02" in table.source

    def test_n_and_k_on_different_ranges_keep_their_overlap(self, tmp_path):
        path = write(tmp_path, """
REFERENCES: test
DATA:
  - type: tabulated n
    data: |
        0.3 1.5
        0.9 1.4
  - type: tabulated k
    data: |
        0.5 0.1
        1.5 0.2
""")
        low, high = materials.from_refractiveindex_file(path).range_nm
        assert (low, high) == (pytest.approx(500), pytest.approx(900))


class TestRefusals:
    @pytest.mark.parametrize(
        "body, match",
        [
            ("  - type: tabulated n2\n    data: |\n        1.0 2.0e-20\n        2.0 2.1e-20\n", "Kerr"),
            ("  - type: tabulated k\n    data: |\n        0.5 0.1\n        0.6 0.1\n", "no n"),
            ("  - type: tabulated n\n    data: |\n        0.5 1.5\n        0.6 1.5\n"
             "  - type: tabulated n\n    data: |\n        0.5 1.5\n        0.6 1.5\n", "n twice"),
            ("  - type: tabulated n\n    data: |\n        0.3 1.5\n        0.4 1.5\n"
             "  - type: tabulated k\n    data: |\n        0.5 0.1\n        0.6 0.1\n", "do not overlap"),
            ("  - type: spectral fit\n    data: x\n", "unrecognised"),
        ],
        ids=["kerr-n2", "k-only", "n-twice", "no-overlap", "unknown-type"],
    )
    def test_refused_with_a_reason(self, tmp_path, body, match):
        path = write(tmp_path, "REFERENCES: test\nDATA:\n" + body)
        with pytest.raises(ValueError, match=match):
            materials.from_refractiveindex_file(path)

    def test_a_file_that_is_not_this_format_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="no DATA"):
            materials.from_refractiveindex_file(write(tmp_path, "just: text\n"))


class TestInALibrary:
    def test_a_yml_file_is_a_coating(self, tmp_path, monkeypatch):
        """A database file is named for its author, so it is saved under the
        material's name to become that coating."""
        (tmp_path / "Au.yml").write_bytes((FIXTURES / "tabulated-nk_Au-Johnson.yml").read_bytes())
        monkeypatch.setenv(materials.LIBRARY_ENV, str(tmp_path))
        table = materials.lookup("Au")
        assert "refractiveindex.info" in table.source
        assert str(tmp_path) in table.source

    def test_so_is_a_yaml_file(self, tmp_path, monkeypatch):
        (tmp_path / "Glass.yaml").write_bytes((FIXTURES / "formula-1_SiO2-Malitson.yml").read_bytes())
        monkeypatch.setenv(materials.LIBRARY_ENV, str(tmp_path))
        assert "Glass" in materials.available()
        assert materials.lookup("Glass").name == "Glass"


class TestTheExtra:
    def test_the_setup_tab_lists_an_unreadable_yaml_rather_than_failing(
        self, tmp_path, monkeypatch
    ):
        """A .yml in a library on an install without PyYAML is a row saying
        so, not a Materials page that raises."""
        from gratinglab.gui.qt.setup_tab import library_page

        (tmp_path / "Au.yml").write_bytes((FIXTURES / "tabulated-nk_Au-Johnson.yml").read_bytes())
        monkeypatch.setenv(materials.LIBRARY_ENV, str(tmp_path))
        monkeypatch.setitem(sys.modules, "yaml", None)
        page = library_page()
        assert "| `Au` | unreadable |" in page and "refractiveindex" in page

    def test_without_pyyaml_the_message_names_the_extra(self, monkeypatch):
        monkeypatch.setitem(sys.modules, "yaml", None)
        with pytest.raises(Exception, match="refractiveindex"):
            read("formula-1_SiO2-Malitson.yml")

    def test_importing_the_package_does_not_import_pyyaml(self):
        """The core install is numpy, scipy and pydantic. PyYAML is imported
        when a YAML file is read and not a moment before."""
        probe = subprocess.run(
            [sys.executable, "-c",
             "import sys, gratinglab, gratinglab.materials; "
             "assert 'yaml' not in sys.modules, 'yaml imported'; print('ok')"],
            capture_output=True, text=True, encoding="utf-8",
        )
        assert probe.returncode == 0, probe.stderr
