"""PCGrate .ggp boundary profiles.

Read liberally (three header variants exist on disk), write strictly (only one
of them loads in PCGrate).
"""

import numpy as np
import pytest

from gratinglab.io.ggp import GGP_HEADER, read_ggp, read_profile, write_ggp
from gratinglab.profiles import FromProfileData, Sinusoidal

from .conftest import reference_dir

CANONICAL = "3 0 - Polygonal type\nPeriod: 1 PSC: 1\n"
POINTS = "0.0 0.0\n0.25 0.2\n0.5 0.3\n0.75 0.1\n1.0 0.0\n"


def write_variant(tmp_path, header, name="p.ggp"):
    path = tmp_path / name
    path.write_text(header + POINTS)
    return path


class TestHeaderVariants:
    def test_canonical(self, tmp_path):
        ggp = read_ggp(write_variant(tmp_path, CANONICAL))
        assert ggp.header_variant == "canonical"
        assert ggp.format_valid
        assert len(ggp.t) == 5

    def test_commented_canonical_parses_but_is_flagged(self, tmp_path):
        """`np.savetxt(header=...)` produces this; PCGrate rejects it."""
        ggp = read_ggp(
            write_variant(tmp_path, "# 3 0 - Polygonal type\n# Period: 1 PSC: 1\n")
        )
        assert ggp.header_variant == "commented-canonical"
        assert not ggp.format_valid
        assert len(ggp.t) == 5

    def test_commented_columns_parses_but_is_flagged(self, tmp_path):
        ggp = read_ggp(
            write_variant(tmp_path, "# X(normalized_0-1) Y(normalized)\n")
        )
        assert ggp.header_variant == "commented-columns"
        assert not ggp.format_valid

    def test_no_header_at_all(self, tmp_path):
        ggp = read_ggp(write_variant(tmp_path, ""))
        assert ggp.header_variant == "none"
        assert not ggp.format_valid

    def test_all_variants_give_identical_data(self, tmp_path):
        datasets = [
            read_ggp(write_variant(tmp_path, header, f"v{i}.ggp")).y
            for i, header in enumerate(
                [CANONICAL, "# 3 0 - Polygonal type\n# Period: 1 PSC: 1\n", ""]
            )
        ]
        for other in datasets[1:]:
            assert np.allclose(datasets[0], other)


class TestRejection:
    def test_rejects_a_truncated_data_line(self, tmp_path):
        """A silently dropped point would shift the profile and change every
        efficiency computed from it."""
        path = tmp_path / "bad.ggp"
        path.write_text(CANONICAL + "0.0 0.0\n0.25 0.2\n0.5\n0.75 0.1\n")
        with pytest.raises(ValueError, match="expected an .x y. pair"):
            read_ggp(path)

    def test_rejects_non_numeric_data(self, tmp_path):
        path = tmp_path / "bad.ggp"
        path.write_text(CANONICAL + "0.0 0.0\n0.25 abc\n0.5 0.3\n")
        with pytest.raises(ValueError, match="expected an .x y. pair"):
            read_ggp(path)

    def test_rejects_empty_file(self, tmp_path):
        path = tmp_path / "empty.ggp"
        path.write_text(CANONICAL)
        with pytest.raises(ValueError, match="no data points"):
            read_ggp(path)

    def test_rejects_too_few_points(self, tmp_path):
        path = tmp_path / "short.ggp"
        path.write_text(CANONICAL + "0.0 0.0\n1.0 0.0\n")
        with pytest.raises(ValueError, match="at least 3"):
            read_ggp(path)

    def test_rejects_header_appearing_after_data(self, tmp_path):
        path = tmp_path / "mixed.ggp"
        path.write_text(CANONICAL + "0.0 0.0\n0.5 0.3\nPeriod: 1\n1.0 0.0\n")
        with pytest.raises(ValueError, match="must precede all data"):
            read_ggp(path)


class TestWriting:
    def test_writes_the_canonical_header(self, tmp_path):
        path = write_ggp(tmp_path / "out.ggp", t=[0.0, 0.5, 1.0], y=[0.0, 0.3, 0.0])
        assert path.read_text().startswith(GGP_HEADER + "\n")
        assert not path.read_text().startswith("#")

    def test_round_trip_preserves_points(self, tmp_path):
        profile = FromProfileData(t=(0.0, 0.3, 0.6, 1.0), y=(0.0, 0.2, 0.35, 0.0))
        recovered = read_profile(write_ggp(tmp_path / "rt.ggp", profile))
        assert np.allclose(recovered.t, profile.t, atol=1e-6)
        assert np.allclose(recovered.y, profile.y, atol=1e-6)

    def test_refuses_to_write_an_undercut_profile(self, tmp_path):
        """PCGrate's polygonal border expects a monotonic boundary."""
        undercut = FromProfileData(t=(0.0, 0.5, 0.3, 0.9), y=(0.0, 0.2, 0.4, 0.1))
        with pytest.raises(ValueError, match="undercut"):
            write_ggp(tmp_path / "u.ggp", undercut)

    def test_rejects_ambiguous_arguments(self, tmp_path):
        profile = FromProfileData(t=(0.0, 0.5, 1.0), y=(0.0, 0.3, 0.0))
        with pytest.raises(ValueError, match="not both or neither"):
            write_ggp(tmp_path / "x.ggp", profile, t=[0, 1], y=[0, 1])
        with pytest.raises(ValueError, match="not both or neither"):
            write_ggp(tmp_path / "x.ggp")

    def test_rejects_mismatched_lengths(self, tmp_path):
        with pytest.raises(ValueError, match="differ in length"):
            write_ggp(tmp_path / "x.ggp", t=[0, 1], y=[0, 1, 2])

    def test_matches_afm_blaze_meas_line_format(self, tmp_path):
        """Byte compatibility with the upstream producer of these files."""
        path = write_ggp(tmp_path / "fmt.ggp", t=[0.0, 0.5], y=[0.0, 0.25])
        assert path.read_text() == GGP_HEADER + "\n0.000000 0.000000\n0.500000 0.250000\n"


class TestProfileIntegration:
    def test_becomes_a_usable_profile(self, tmp_path):
        """A .ggp must support every representation an analytic profile does."""
        exact = Sinusoidal(depth_fraction=0.3)
        t = np.linspace(0.0, 1.0, 200, endpoint=False)
        path = write_ggp(tmp_path / "sin.ggp", t=t, y=exact.height(t))

        profile = read_profile(path)
        query = np.linspace(0.0, 1.0, 41, endpoint=False)
        assert np.allclose(profile.height(query), exact.height(query), atol=2e-4)
        assert profile.depth == pytest.approx(exact.depth, rel=1e-3)
        assert len(profile.slice_layers(8)) == 8
        assert len(profile.boundary(64).t) == 64

    def test_heights_are_referenced_to_the_valley(self, tmp_path):
        path = write_ggp(tmp_path / "off.ggp", t=[0.0, 0.5, 1.0], y=[0.7, 1.0, 0.7])
        profile = read_profile(path)
        assert min(profile.y) == pytest.approx(0.0)
        assert profile.depth == pytest.approx(0.3)


REF = reference_dir()
requires_corpus = pytest.mark.skipif(REF is None, reason="reference corpus unavailable")


@requires_corpus
class TestRealFiles:
    def test_panter1_profile(self, ref_dir):
        """The exact AFM profile the panter1 PCGrate run used."""
        path = ref_dir / "PCGrateProjects" / "panter1.ggp"
        if not path.exists():
            pytest.skip("panter1.ggp not in corpus")

        ggp = read_ggp(path)
        assert ggp.format_valid
        assert ggp.header_variant == "canonical"
        assert len(ggp.t) == 82
        assert ggp.depth == pytest.approx(0.3119, abs=1e-4)

        profile = read_profile(path)
        assert profile.is_single_valued()
        assert profile.apex == pytest.approx(0.383, abs=0.01)
        # A rounded AFM apex, not an ideal sawtooth: the peak is interior.
        assert 0.0 < profile.apex < 1.0

    def test_known_defect_afm_test_ggp_is_truncated(self, ref_dir):
        """`OGRE/AFM_test.ggp` line 164 has a missing x value.

        Recorded so it is not rediscovered. If the file is ever regenerated,
        this test fails and should be deleted.
        """
        path = ref_dir / "OGRE" / "AFM_test.ggp"
        if not path.exists():
            pytest.skip("AFM_test.ggp not in corpus")
        with pytest.raises(ValueError, match="164.*expected an"):
            read_ggp(path)
