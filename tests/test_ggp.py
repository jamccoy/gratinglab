"""PCGrate .ggp boundary profiles.

Read liberally (five control-string variants exist on disk), write one form:
the points alone by default, the two control strings the User's Guide of
PCGrate-S(X) 6.1 documents in 10.2.2 (p. 63) on request.

The vendor ships 35 polygonal files of its own, and `TestPcgrateOwnFiles`
checks this module against them rather than against what this module used to
write. They are not committed; see `corpus.pcgrate_profiles_dir`.
"""

import numpy as np
import pytest

from gratinglab.io.ggp import GGP_HEADER, read_ggp, read_profile, write_ggp
from gratinglab.profiles import FromProfileData, Sinusoidal

from corpus import pcgrate_profiles_dir, reference_dir

CANONICAL = "3 0 - 20 harmonics - Polygonal type\nPeriod: 1 PSC: 1\n"
SHORT = "3 0 - Polygonal type\nPeriod: 1 PSC: 1\n"
POINTS = "0.0 0.0\n0.25 0.2\n0.5 0.3\n0.75 0.1\n1.0 0.0\n"


def write_variant(tmp_path, header, name="p.ggp"):
    path = tmp_path / name
    path.write_text(header + POINTS, encoding="utf-8")
    return path


class TestHeaderVariants:
    def test_canonical(self, tmp_path):
        """The documented pair: 10.2.2 step 2, and every file PCGrate ships."""
        ggp = read_ggp(write_variant(tmp_path, CANONICAL))
        assert ggp.header_variant == "canonical"
        assert ggp.format_valid
        assert len(ggp.t) == 5

    def test_the_short_type_line_parses_but_is_flagged(self, tmp_path):
        """What this module wrote until 2026-10-30.

        It lacks the `- 20 harmonics -` segment, so it is a control string "in
        another format": legal through the Border Profile Editor, rejected by a
        direct load. This file was labelled `canonical` and valid until the
        guide was read.
        """
        ggp = read_ggp(write_variant(tmp_path, SHORT))
        assert ggp.header_variant == "short-type-line"
        assert not ggp.format_valid
        assert len(ggp.t) == 5

    def test_commented_type_line_parses_but_is_flagged(self, tmp_path):
        """`np.savetxt(header=...)` produces this; a direct load rejects it."""
        ggp = read_ggp(
            write_variant(
                tmp_path,
                "# 3 0 - 20 harmonics - Polygonal type\n# Period: 1 PSC: 1\n",
            )
        )
        assert ggp.header_variant == "commented-type-line"
        assert not ggp.format_valid
        assert len(ggp.t) == 5

    def test_a_canonical_type_line_without_a_period_is_flagged(self, tmp_path):
        """10.2.3's second error: the second control string is empty."""
        ggp = read_ggp(write_variant(tmp_path, "3 0 - 20 harmonics - Polygonal type\n"))
        assert ggp.header_variant == "canonical"
        assert not ggp.format_valid

    def test_a_tab_is_an_acceptable_separator(self, tmp_path):
        """10.2.2 says space, and 24 of the 35 files PCGrate ships use a tab."""
        path = tmp_path / "tabbed.ggp"
        path.write_text(CANONICAL + POINTS.replace(" ", "\t"), encoding="utf-8")
        ggp = read_ggp(path)
        assert ggp.format_valid
        assert list(ggp.t) == [0.0, 0.25, 0.5, 0.75, 1.0]

    def test_commented_columns_parses_but_is_flagged(self, tmp_path):
        ggp = read_ggp(
            write_variant(tmp_path, "# X(normalized_0-1) Y(normalized)\n")
        )
        assert ggp.header_variant == "commented-columns"
        assert not ggp.format_valid

    def test_no_header_at_all(self, tmp_path):
        """The form the writer produces by default."""
        ggp = read_ggp(write_variant(tmp_path, ""))
        assert ggp.header_variant == "none"
        assert ggp.format_valid

    def test_all_variants_give_identical_data(self, tmp_path):
        datasets = [
            read_ggp(write_variant(tmp_path, header, f"v{i}.ggp")).y
            for i, header in enumerate(
                [CANONICAL, SHORT, "# 3 0 - Polygonal type\n# Period: 1 PSC: 1\n", ""]
            )
        ]
        for other in datasets[1:]:
            assert np.allclose(datasets[0], other)


class TestEncoding:
    """PCGrate writes cp1252 and `write_ggp` writes UTF-8, and which one a
    reader assumed used to be the platform's locale. A degree sign in a header
    is where the two differ."""

    HEADER = "3 0 - 20 harmonics - Polygonal type\nPeriod: 1 PSC: 1 blaze 29.5°\n"

    @pytest.mark.parametrize("encoding", ["utf-8", "cp1252"])
    def test_a_non_ascii_header_reads_the_same_either_way(self, tmp_path, encoding):
        path = tmp_path / "p.ggp"
        path.write_bytes((self.HEADER + POINTS).encode(encoding))
        reference = tmp_path / "ascii.ggp"
        reference.write_text(CANONICAL + POINTS, encoding="utf-8")
        np.testing.assert_array_equal(read_ggp(path).t, read_ggp(reference).t)
        np.testing.assert_array_equal(read_ggp(path).y, read_ggp(reference).y)


class TestRejection:
    def test_rejects_a_truncated_data_line(self, tmp_path):
        """A silently dropped point would shift the profile and change every
        efficiency computed from it."""
        path = tmp_path / "bad.ggp"
        path.write_text(CANONICAL + "0.0 0.0\n0.25 0.2\n0.5\n0.75 0.1\n", encoding="utf-8")
        with pytest.raises(ValueError, match="expected an .x y. pair"):
            read_ggp(path)

    def test_rejects_non_numeric_data(self, tmp_path):
        path = tmp_path / "bad.ggp"
        path.write_text(CANONICAL + "0.0 0.0\n0.25 abc\n0.5 0.3\n", encoding="utf-8")
        with pytest.raises(ValueError, match="expected an .x y. pair"):
            read_ggp(path)

    def test_rejects_empty_file(self, tmp_path):
        path = tmp_path / "empty.ggp"
        path.write_text(CANONICAL, encoding="utf-8")
        with pytest.raises(ValueError, match="no data points"):
            read_ggp(path)

    def test_rejects_too_few_points(self, tmp_path):
        path = tmp_path / "short.ggp"
        path.write_text(CANONICAL + "0.0 0.0\n1.0 0.0\n", encoding="utf-8")
        with pytest.raises(ValueError, match="at least 3"):
            read_ggp(path)

    def test_rejects_header_appearing_after_data(self, tmp_path):
        path = tmp_path / "mixed.ggp"
        path.write_text(CANONICAL + "0.0 0.0\n0.5 0.3\nPeriod: 1\n1.0 0.0\n", encoding="utf-8")
        with pytest.raises(ValueError, match="must precede all data"):
            read_ggp(path)


class TestWriting:
    def test_writes_no_header_by_default(self, tmp_path):
        """Aditi: PCGrate's border-profile import does not want one (2026-10-29)."""
        path = write_ggp(tmp_path / "out.ggp", t=[0.0, 0.5, 1.0], y=[0.0, 0.3, 0.0])
        assert path.read_text(encoding="utf-8").startswith("0.000000 0.000000\n")
        assert read_ggp(path).header_variant == "none"

    def test_writes_the_documented_header_on_request(self, tmp_path):
        path = write_ggp(tmp_path / "out.ggp", t=[0.0, 0.5, 1.0], y=[0.0, 0.3, 0.0],
                         header=True)
        assert path.read_text(encoding="utf-8").startswith(GGP_HEADER + "\n")
        assert not path.read_text(encoding="utf-8").startswith("#")

    def test_round_trip_preserves_points(self, tmp_path):
        profile = FromProfileData(t=(0.0, 0.3, 0.6, 1.0), y=(0.0, 0.2, 0.35, 0.0))
        recovered = read_profile(write_ggp(tmp_path / "rt.ggp", profile))
        assert np.allclose(recovered.t, profile.t, atol=1e-6)
        assert np.allclose(recovered.y, profile.y, atol=1e-6)

    def test_the_file_runs_along_the_dispersion_direction(self, tmp_path):
        """A ``.ggp``'s x runs along +d-hat and t against it, so a facet
        that rises with x in the file descends with t in the profile. The
        PCGrate TASTE match is what fixes the file's direction; see
        docs/findings.md. A symmetric groove could not tell the two apart."""
        path = write_ggp(tmp_path / "saw.ggp", t=[0.0, 0.8, 1.0], y=[0.0, 0.4, 0.0])
        profile = read_profile(path)
        np.testing.assert_allclose(profile.t, [0.0, 0.2, 1.0], atol=1e-6)
        np.testing.assert_allclose(profile.y, [0.0, 0.4, 0.0], atol=1e-6)
        written = write_ggp(tmp_path / "back.ggp", profile).read_text(encoding="utf-8")
        assert "0.800000 0.400000" in written

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

    def test_writes_the_exact_pcgrate_line_format(self, tmp_path):
        """Six decimals, one space, no break after the last point.

        PCGrate rejects a file ending in a line break (found by Aditi, 2026-10-27),
        which is why the last point is not followed by one.

        This began as a byte-compatibility check against `afm_blaze_meas`,
        which shipped a second `write_ggp`. That copy is gone and this module
        is now the only implementation, so there is no longer another writer
        to agree with. The assertion outlives the reason for it: the format
        is still a fact about PCGrate, and still worth pinning.
        """
        path = write_ggp(tmp_path / "fmt.ggp", t=[0.0, 0.5], y=[0.0, 0.25])
        assert path.read_text(encoding="utf-8") == "0.000000 0.000000\n0.500000 0.250000"
        assert not path.read_bytes().endswith((b"\n", b"\r"))

    def test_header_can_be_left_out_explicitly(self, tmp_path):
        path = write_ggp(tmp_path / "bare.ggp", t=[0.0, 0.5, 1.0], y=[0.0, 0.25, 0.0],
                         header=False)
        assert path.read_text(encoding="utf-8") == \
            "0.000000 0.000000\n0.500000 0.250000\n1.000000 0.000000"
        ggp = read_ggp(path)
        assert ggp.header_variant == "none"
        assert list(ggp.t) == [0.0, 0.5, 1.0]


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
        # `short-type-line`, and no longer valid: this file carries
        # `3 0 - Polygonal type`, which this project wrote and which lacks the
        # `- 20 harmonics -` segment PCGrate's first control string needs. It
        # was labelled `canonical` until the guide was read. The points are
        # unaffected, and so is every number below.
        assert ggp.header_variant == "short-type-line"
        assert not ggp.format_valid
        assert len(ggp.t) == 82
        assert ggp.depth == pytest.approx(0.3119, abs=1e-4)

        profile = read_profile(path)
        assert profile.is_single_valued()
        # 0.617 and not the 0.383 this asserted from M1 until 2026-10-29: the
        # file's x runs along +d-hat and `to_profile` reflects it into t, so the
        # apex sits at 1 - 0.383. The expectation was written before that
        # reflection existed and was never moved with it. The two numbers
        # summing to one is the whole of the difference.
        assert profile.apex == pytest.approx(0.617, abs=0.01)
        # A rounded AFM apex, not an ideal sawtooth: the peak is interior.
        assert 0.0 < profile.apex < 1.0
        # And the reflection is invertible, so the file's own apex is still
        # where the file puts it.
        x, y = profile.dispersion_axis()
        assert x[y.argmax()] == pytest.approx(0.383, abs=0.01)

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


PCGRATE = pcgrate_profiles_dir()
requires_pcgrate = pytest.mark.skipif(
    PCGRATE is None,
    reason="PCGrate's own BorderProfiles unavailable; "
           "set GRATINGLAB_PCGRATE_PROFILES_DIR",
)


def vendor_polygonal_files():
    """The border profiles PCGrate ships that are of the polygonal type.

    Selected on the type code at the start of the first control string -- `3 0`
    is polygonal, while the trapezoidal, trigonometric and sine families start
    `0 n`, `1 n` and `2 0`. Deliberately not selected on the whole caption:
    the caption is what `TestPcgrateOwnFiles` is checking, and a selector that
    read it would make that check circular.
    """
    if PCGRATE is None:
        return []
    files = []
    for path in sorted(PCGRATE.glob("*.ggp")):
        first = path.read_bytes().decode("latin-1").splitlines()[0]
        if first.split()[:2] == ["3", "0"]:
            files.append(path)
    return files


@requires_pcgrate
class TestPcgrateOwnFiles:
    """The vendor's own files, which are the authority on what it writes.

    Until 2026-10-30 every one of them read back as `commented-canonical` and
    `format_valid is False`, despite containing no `#` at all: `_classify`
    compared the first control string against the short form this project
    wrote, and the real one carries a `- 20 harmonics -` segment the short form
    omits. Both the label and the verdict were backwards on the only files that
    settle the question.
    """

    def test_every_polygonal_file_is_canonical_and_valid(self):
        files = vendor_polygonal_files()
        assert len(files) >= 30, f"only {len(files)} polygonal files; corpus thinned?"
        wrong = {}
        for path in files:
            ggp = read_ggp(path)
            if (ggp.header_variant, ggp.format_valid) != ("canonical", True):
                wrong[path.name] = (ggp.header_variant, ggp.format_valid)
        assert not wrong

    def test_the_first_control_string_is_the_one_we_write(self, tmp_path):
        """Byte for byte, and identical across all 35."""
        firsts = {
            path.read_bytes().split(b"\n")[0].rstrip(b"\r")
            for path in vendor_polygonal_files()
        }
        assert firsts == {b"3 0 - 20 harmonics - Polygonal type"}

        path = write_ggp(tmp_path / "ours.ggp", t=[0.0, 0.5, 1.0], y=[0.0, 0.3, 0.0],
                         header=True)
        assert path.read_bytes().split(b"\n")[0] == firsts.pop()

    def test_the_second_control_string_always_carries_a_period(self):
        """With or without the PSC, which 4 of the 35 leave out."""
        seconds = [
            path.read_bytes().decode("latin-1").splitlines()[1]
            for path in vendor_polygonal_files()
        ]
        assert all(line.startswith("Period:") for line in seconds)

    def test_every_polygonal_file_obeys_the_documented_geometry(self):
        """The rules of 10.2.3 are a closed list; the vendor's files keep it."""
        offenders = {
            path.name: violations
            for path in vendor_polygonal_files()
            if (violations := documented_violations(
                path.read_bytes().decode("latin-1"), control_strings=2))
        }
        assert not offenders


def documented_violations(text, *, control_strings):
    """The errors of 10.2.3 (p. 63) that `text` would be rejected for.

    The thirteen the guide lists, less the three about the control strings
    themselves, which `read_ggp` classifies instead. A tab separator is
    accepted: the guide asks for a space, and 24 of the files PCGrate ships use
    a tab, so the program reads more than it documents.
    """
    lines = text.splitlines()[control_strings:]
    data = [line for line in lines if line.strip()]
    problems = []

    if not 2 <= len(data) <= 5000:
        problems.append(f"{len(data)} data strings, outside 2 .. 5000")

    points = []
    for number, line in enumerate(data, start=1):
        if "," in line:
            problems.append(f"line {number}: a comma, not a decimal point")
            continue
        fields = line.split()
        if len(fields) != 2:
            problems.append(f"line {number}: {len(fields)} numbers, not 2")
            continue
        points.append((float(fields[0]), float(fields[1])))

    if not points:
        return problems or ["no data points"]

    x = [point[0] for point in points]
    y = [point[1] for point in points]
    if min(x) < 0.0 or max(x) > 1.0:
        problems.append(f"x outside 0 .. 1: [{min(x):g}, {max(x):g}]")
    if min(y) < 0.0 or max(y) > 100.0:
        problems.append(f"y outside 0 .. 100: [{min(y):g}, {max(y):g}]")
    if any(later < earlier for earlier, later in zip(x, x[1:])):
        problems.append("x decreases")
    if (x[0], y[0]) != (0.0, 0.0):
        problems.append(f"first point {(x[0], y[0])}, not (0, 0)")
    if (x[-1], y[-1]) != (1.0, 0.0):
        problems.append(f"last point {(x[-1], y[-1])}, not (1, 0)")
    return problems


class TestWrittenFilesObeyTheDocumentedGeometry:
    """What `write_ggp` produces, against the closed list of 10.2.3 (p. 63).

    Nothing here was failing: the writer already satisfied every rule. The
    point is that the rules are now known rather than guessed, so a future
    change to the format, the point count or the resampling cannot quietly
    break one.
    """

    @staticmethod
    def groove(n_points):
        """A groove that starts and ends on the documented endpoints."""
        x = np.linspace(0.0, 1.0, n_points)
        return x, 0.3 * np.sin(np.pi * x) ** 2

    @pytest.mark.parametrize("n_points", [2, 3, 100, 2000, 5000])
    @pytest.mark.parametrize("header", [False, True])
    def test_output_has_no_documented_violation(self, tmp_path, n_points, header):
        x, y = self.groove(n_points)
        path = write_ggp(tmp_path / "g.ggp", t=x, y=y, header=header)
        text = path.read_text(encoding="utf-8")
        assert documented_violations(text, control_strings=2 if header else 0) == []

    def test_a_single_space_separates_the_two_numbers(self, tmp_path):
        """10.2.2: "to separate X and Y values in strings use the space key"."""
        x, y = self.groove(50)
        text = write_ggp(tmp_path / "g.ggp", t=x, y=y).read_text(encoding="utf-8")
        assert all(line.count(" ") == 1 for line in text.splitlines())
        assert "\t" not in text

    def test_at_most_sixteen_digits_a_piece(self, tmp_path):
        """10.2.2 caps a number at 16 digits including the decimal point."""
        x, y = self.groove(50)
        text = write_ggp(tmp_path / "g.ggp", t=x, y=y).read_text(encoding="utf-8")
        assert max(len(field) for line in text.splitlines()
                   for field in line.split()) <= 16

    def test_the_point_count_is_bounded_by_the_setting(self):
        """`ggp_n_points` is validated against the same 5000 (10.2.3)."""
        from gratinglab.metrology.settings import MAX_GGP_POINTS, AnalysisSettings

        assert MAX_GGP_POINTS == 5000
        base = AnalysisSettings.from_config()
        assert 'ggp_n_points' not in [
            field for field, _ in base.with_(ggp_n_points=5000).validate()
        ]
        assert 'ggp_n_points' in [
            field for field, _ in base.with_(ggp_n_points=5001).validate()
        ]
