"""gratinglab's own boundary files: the period survives, and nothing else is lost.

A stand-in profile object is used so these need no metrology extra; the
round trip against a real ``BoundaryProfile`` is in
``tests/metrology/test_boundary_pipeline.py``.
"""

import json
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import pytest

from gratinglab import __version__
from gratinglab.io.boundary import (
    FORMAT,
    VERSION,
    read_boundary_json,
    read_profile_file,
    write_boundary_csv,
    write_boundary_json,
)
from gratinglab.io.ggp import write_ggp


def _profile(n=41, period=314.97):
    t = np.linspace(0.0, 1.0, n)
    y = np.where(t < 0.7, t / 0.7, (1 - t) / 0.3) * 0.31   # a sawtooth, in period units
    y[0] = y[-1] = 0.0
    return SimpleNamespace(
        x_norm=t, y_norm=y, period_nm=np.float64(period),
        n_grooves=6, n_used=5, n_edge_rejected=1,
        metrics={"groove_depth": np.float64(0.31), "max_angle_deg": 66.0, "n_grooves": 5},
    )


@dataclass(frozen=True)
class _Settings:
    edge_exclusion_periods: float = 0.6
    ggp_n_points: int = 41


class TestJson:
    def test_it_round_trips_exactly(self, tmp_path):
        """Full precision, unlike a .ggp's six decimals."""
        p = _profile()
        back = read_boundary_json(write_boundary_json(tmp_path / "g.json", p))
        assert back.t == tuple(p.x_norm)
        assert back.y == tuple(p.y_norm)
        assert back.period_nm == p.period_nm

    def test_the_period_and_groove_count_reach_the_problem(self, tmp_path):
        back = read_boundary_json(write_boundary_json(tmp_path / "g.json", _profile()))
        problem = back.to_problem(coating="Au")
        assert problem.period == pytest.approx(314.97)
        assert problem.n_grooves == 5
        assert problem.coating == "Au"
        assert problem.profile.depth == pytest.approx(0.31)

    def test_provenance_records_how_it_was_made(self, tmp_path):
        path = write_boundary_json(tmp_path / "g.json", _profile(),
                                   source="scan_007.txt", settings=_Settings())
        prov = read_boundary_json(path).provenance
        assert prov["gratinglab_version"] == __version__
        assert prov["source"] == "scan_007.txt"
        assert (prov["grooves_detected"], prov["grooves_averaged"],
                prov["grooves_edge_rejected"]) == (6, 5, 1)
        assert prov["settings"]["edge_exclusion_periods"] == 0.6

    def test_metrics_survive_without_numpy_types(self, tmp_path):
        back = read_boundary_json(write_boundary_json(tmp_path / "g.json", _profile()))
        assert back.metrics["groove_depth"] == pytest.approx(0.31)

    def test_a_fitted_blaze_angle_is_carried(self, tmp_path):
        back = read_boundary_json(write_boundary_json(tmp_path / "g.json", _profile(), blaze_angle=29.5))
        assert back.profile().blaze_angle == 29.5

    def test_it_is_utf8_json_naming_its_format(self, tmp_path):
        path = write_boundary_json(tmp_path / "g.json", _profile(), source="scène.txt")
        document = json.loads(path.read_bytes().decode("utf-8"))
        assert document["format"] == FORMAT and document["version"] == VERSION
        assert b"\r\n" not in path.read_bytes()


class TestJsonRefusals:
    def _write(self, tmp_path, document):
        path = tmp_path / "x.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        return path

    def test_an_unrelated_json_is_refused_by_name(self, tmp_path):
        with pytest.raises(ValueError, match="not a gratinglab boundary profile"):
            read_boundary_json(self._write(tmp_path, {"hello": 1}))

    def test_a_newer_version_says_to_upgrade(self, tmp_path):
        doc = {"format": FORMAT, "version": VERSION + 1, "period_nm": 300.0,
               "profile": {"t": [0, 1], "y": [0, 0]}}
        with pytest.raises(ValueError, match="upgrade"):
            read_boundary_json(self._write(tmp_path, doc))

    @pytest.mark.parametrize("period", [0, -5, None, "315"])
    def test_a_bad_period_is_refused(self, tmp_path, period):
        doc = {"format": FORMAT, "version": VERSION, "period_nm": period,
               "profile": {"t": [0, 1], "y": [0, 0]}}
        with pytest.raises(ValueError, match="period_nm"):
            read_boundary_json(self._write(tmp_path, doc))

    def test_unequal_columns_are_refused(self, tmp_path):
        doc = {"format": FORMAT, "version": VERSION, "period_nm": 300.0,
               "profile": {"t": [0, 0.5, 1], "y": [0, 0]}}
        with pytest.raises(ValueError, match="equal length"):
            read_boundary_json(self._write(tmp_path, doc))

    def test_not_json_at_all(self, tmp_path):
        path = tmp_path / "x.json"
        path.write_text("3 0 - Polygonal type\n", encoding="utf-8")
        with pytest.raises(ValueError, match="not a JSON file"):
            read_boundary_json(path)


class TestCsv:
    def test_columns_are_fractions_and_nanometres(self, tmp_path):
        p = _profile()
        path = write_boundary_csv(tmp_path / "g.csv", p, source="scan_007.txt")
        text = path.read_bytes().decode("utf-8")
        header = [line for line in text.splitlines() if line.startswith("#")]
        assert any("period_nm: 314.97" in line for line in header)
        assert any("5 averaged of 6 detected" in line for line in header)
        assert any("scan_007.txt" in line for line in header)
        table = np.genfromtxt(path, delimiter=",", names=True, comments="#", encoding="utf-8")
        np.testing.assert_allclose(table["t"], p.x_norm)
        np.testing.assert_allclose(table["y"], p.y_norm)
        np.testing.assert_allclose(table["x_nm"], p.x_norm * 314.97)
        np.testing.assert_allclose(table["y_nm"], p.y_norm * 314.97)

    def test_pandas_reads_it_too(self, tmp_path):
        """The two readers a student is likeliest to reach for disagree about
        comment lines; header-first is the layout both accept."""
        pd = pytest.importorskip("pandas")
        p = _profile()
        frame = pd.read_csv(write_boundary_csv(tmp_path / "g.csv", p), comment="#")
        assert list(frame.columns) == ["t", "y", "x_nm", "y_nm"]
        np.testing.assert_allclose(frame["y_nm"], p.y_norm * 314.97)


class TestEitherFormat:
    def test_a_json_reads_as_a_profile(self, tmp_path):
        p = _profile()
        profile = read_profile_file(write_boundary_json(tmp_path / "g.json", p))
        assert profile.t == tuple(p.x_norm)

    def test_a_ggp_still_does(self, tmp_path):
        p = _profile()
        path = write_ggp(tmp_path / "g.ggp", t=p.x_norm, y=p.y_norm)
        np.testing.assert_allclose(read_profile_file(path).t, p.x_norm, atol=1e-6)
