"""
Reading raw Nanoscope .spm files.

The fixture pair is the point of this file: `TASTE_ALS_A205_Ti_Pt.0_00003.spm`
and `TASTE_ALS_A205_Ti_Pt_flatten.txt` are the same scan by two routes - the
instrument's own file, and the Gwyddion export that has been this project's only
input until now. They must agree.
"""
import contextlib
import io
import os
import sys

import numpy as np
import pytest

from scans import real_scan

# src/ is placed on the path by tests/conftest.py; repeated here so the file
# also runs directly as a script, not only under pytest.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src'))

import matplotlib
matplotlib.use('Agg')

from gratinglab.metrology.core.processing import load_afm_data          # noqa: E402
from gratinglab.metrology.io.spm import (                               # noqa: E402
    is_nanoscope_file, list_channels, read_spm)

# A real Nanoscope binary and its Gwyddion export of the same scan. The
# whole point of this module is the binary format and the Sens. Zsens trap,
# so there is nothing to synthesise: these skip without the group's data.
SPM_NAME = 'TASTE_ALS_A205_Ti_Pt.0_00003.spm'
TXT_NAME = 'TASTE_ALS_A205_Ti_Pt_flatten.txt'

_cache = {}


def _txt():
    if 'txt' not in _cache:
        _cache['txt'] = np.genfromtxt(real_scan(TXT_NAME), skip_header=4)
    return _cache['txt']


# ── Header parsing ───────────────────────────────────────────────────────────

def test_recognises_a_nanoscope_file_by_content():
    """
    By content, not by extension.

    Nanoscope also writes companion files with no extension at all
    (`sample_flatten.0_00003`), which are the same format.
    """
    assert is_nanoscope_file(real_scan(SPM_NAME))
    assert not is_nanoscope_file(real_scan(TXT_NAME))


def test_lists_every_plane():
    channels = list_channels(real_scan(SPM_NAME))
    assert len(channels) == 4, [c.describe() for c in channels]
    names = {(c.channel, c.direction) for c in channels}
    assert names == {("Height Sensor", "Retrace"), ("Height Sensor", "Trace"),
                     ("Peak Force Error", "Retrace"), ("Peak Force Error", "Trace")}


def test_dimensions_are_not_assumed_square():
    """
    This scan is 512 wide by 401 lines - it stopped early.

    Assuming Samps/line applies to both axes would read past the plane and
    silently reshape another channel's bytes into the height map.
    """
    height = next(c for c in list_channels(real_scan(SPM_NAME)) if c.channel == "Height Sensor")
    assert (height.rows, height.cols) == (401, 512)
    data, _ = read_spm(real_scan(SPM_NAME))
    assert data.shape == (401, 512)


def test_scan_size_is_read_in_microns():
    """The header writes '2.000 2.000 ~m', where ~m is a mangled µm"""
    _, scan_x_size = read_spm(real_scan(SPM_NAME))
    assert np.isclose(scan_x_size, 2.0)


# ── The height scale: the trap this format sets ──────────────────────────────

def test_height_scale_matches_the_known_export():
    """
    The single most dangerous detail in the format.

    '@2:Z scale' names its sensitivity parameter in brackets - here
    '[Sens. ZsensSens]'. The file also defines '\\@Sens. Zsens', which looks like
    the same name with a redundant suffix and is 5.13x smaller (32.46 vs 166.63
    nm/V). Using it produces an image of the right shape, with grooves in the
    right places, and every blaze angle wrong by a factor of five.

    Peak-to-peak against a known-good export is what catches that.
    """
    data, _ = read_spm(real_scan(SPM_NAME))
    ratio = np.ptp(data) / np.ptp(_txt())
    assert np.isclose(ratio, 1.0, atol=0.001), (
        f"height scale is off by {ratio:.3f}x. If this is ~0.195 or ~5.13, the "
        f"reader resolved '@2:Z scale' to '\\@Sens. Zsens' instead of the "
        f"bracketed '\\@Sens. ZsensSens' - the bracket name must be used verbatim.")


def test_heights_are_returned_in_metres():
    """Text exports are in metres and everything downstream assumes it"""
    data, _ = read_spm(real_scan(SPM_NAME))
    assert np.ptp(data) < 1e-6, "peak-to-peak over a micron suggests nanometres"
    assert np.ptp(data) > 1e-9


# ── Which plane the export came from ─────────────────────────────────────────

def test_retrace_is_the_plane_the_export_came_from():
    """
    Retrace is the default because it is what was exported, not by convention.

    Both directions image the same grooves, so both correlate; Retrace correlates
    better, and that is the evidence for the default.
    """
    txt = _txt()
    retrace, _ = read_spm(real_scan(SPM_NAME), direction="Retrace")
    trace, _ = read_spm(real_scan(SPM_NAME), direction="Trace")

    corr_retrace = np.corrcoef(retrace.ravel(), txt.ravel())[0, 1]
    corr_trace = np.corrcoef(trace.ravel(), txt.ravel())[0, 1]

    assert corr_retrace > 0.95, f"Retrace correlation only {corr_retrace:.3f}"
    assert corr_retrace > corr_trace, (
        f"Trace ({corr_trace:.3f}) matched the export better than Retrace "
        f"({corr_retrace:.3f}); the default direction may be wrong")


def test_correlation_is_below_one_because_the_export_was_flattened():
    """
    Not a defect - a fact worth pinning.

    The .txt was flattened by Gwyddion on export; the .spm is raw. A correlation
    of exactly 1.0 would mean the export applied nothing, which would make the
    'flatten' in its filename a lie.
    """
    retrace, _ = read_spm(real_scan(SPM_NAME))
    corr = np.corrcoef(retrace.ravel(), _txt().ravel())[0, 1]
    assert 0.95 < corr < 0.9999, f"correlation {corr:.5f}"


# ── Errors ───────────────────────────────────────────────────────────────────

def test_missing_channel_names_what_the_file_does_contain():
    try:
        read_spm(real_scan(SPM_NAME), channel="Nonexistent Channel")
    except ValueError as exc:
        assert "Height Sensor" in str(exc), "the error should list real channels"
        return
    raise AssertionError("expected a ValueError for a missing channel")


def test_reading_a_text_file_as_spm_is_refused():
    try:
        read_spm(real_scan(TXT_NAME))
    except ValueError as exc:
        assert "Nanoscope" in str(exc)
        return
    raise AssertionError("expected a ValueError for a non-Nanoscope file")


# ── Dispatch and end-to-end agreement ────────────────────────────────────────

def test_load_afm_data_dispatches_on_content():
    """One loader, six call sites, neither of which should care about format"""
    from_spm, scan_spm = load_afm_data(real_scan(SPM_NAME), default_scan_size=2.0)
    from_txt, scan_txt = load_afm_data(real_scan(TXT_NAME), default_scan_size=2.0)
    assert from_spm.shape == from_txt.shape
    assert np.isclose(scan_spm, scan_txt)


def test_both_routes_give_the_same_blaze_angle():
    """
    The measurement that matters: one scan, two routes.

    They are not bit-identical - the export was flattened and the .spm is raw -
    but the software flattens each row group itself, which removes what Gwyddion
    removed. Measured difference is 0.002 degrees against a sigma of 2.13.
    """
    from gratinglab.metrology.analyzer import analyze_single_file
    from gratinglab.metrology.settings import AnalysisSettings

    settings = AnalysisSettings.from_config()
    with contextlib.redirect_stdout(io.StringIO()):
        via_spm = analyze_single_file(real_scan(SPM_NAME), show_plots=False, settings=settings)
        via_txt = analyze_single_file(real_scan(TXT_NAME), show_plots=False, settings=settings)

    assert via_spm['n_grooves'] == via_txt['n_grooves']
    difference = abs(via_spm['mean_angle'] - via_txt['mean_angle'])
    assert difference < 0.05, (
        f"the two routes disagree by {difference:.4f} deg. Measured at 0.002; "
        f"anything approaching the tolerance means the reader or the flattening "
        f"has changed.")


def test_direction_choice_reaches_the_analysis():
    """A setting that does not change the answer is not wired up"""
    from gratinglab.metrology.analyzer import analyze_single_file
    from gratinglab.metrology.settings import AnalysisSettings

    base = AnalysisSettings.from_config()
    with contextlib.redirect_stdout(io.StringIO()):
        retrace = analyze_single_file(real_scan(SPM_NAME), show_plots=False,
                                      settings=base.with_(spm_direction='Retrace'))
        trace = analyze_single_file(real_scan(SPM_NAME), show_plots=False,
                                    settings=base.with_(spm_direction='Trace'))
    assert retrace['mean_angle'] != trace['mean_angle']


if __name__ == '__main__':
    failures = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith('test_'):
            continue
        try:
            fn()
            print(f"PASS  {name}")
        except AssertionError as exc:
            failures += 1
            print(f"FAIL  {name}: {exc}")
    print(f"\n{'all tests passed' if not failures else f'{failures} failure(s)'}")
    sys.exit(1 if failures else 0)


class TestTheSlowAxisSize:
    """
    `\\Aspect Ratio` is how Nanoscope records a rectangular scan.

    `Scan Size` carries the same number twice even when the scan is not
    square, so a reader that takes the second value at face value makes the
    slow axis too long by the aspect ratio. That was true of 1013 of the 1245
    scans in this project's archive, 911 of them at 4:1, and since a groove
    rotation is an arctangent of a ratio of the two axes, every one of those
    reported a rotation too small by the same factor.

    Nothing in the committed fixtures could catch it: the TASTE scan and the
    synthetic one are both 1:1. So these build headers instead, with the field
    values taken from real files, and the last test checks the real archive
    when it is reachable.
    """

    def _header(self, fast_um, slow_um, aspect, lines, plane_aspect=None):
        """A header with the fields `_slow_size_um` reads, and nothing else."""
        text = (f"\\Scan Size: {fast_um} {slow_um} ~m\n"
                f"\\Aspect Ratio: {aspect}\n"
                f"\\Lines: {lines}\n")
        section = "\\Samps/line: 512\n"
        if plane_aspect is not None:
            # Each plane section carries a field of the same name that is not
            # a scan shape. See _aspect_ratio.
            section += f"\\Aspect Ratio: {plane_aspect}\n"
        return section, text

    def test_a_square_scan_is_unchanged(self):
        from gratinglab.metrology.io.spm import _slow_size_um

        section, header = self._header(2, 2, "1:1", 512)
        assert _slow_size_um(section, header, 512) == 2.0

    def test_a_rectangular_scan_is_divided_by_its_aspect_ratio(self):
        """
        The case Gwyddion settles: 1.5 µm at 4:1 is 0.375 µm of slow axis.

        Its export of such a file says `# Width: 1.500 µm` against
        `# Height: 0.375 µm`.
        """
        from gratinglab.metrology.io.spm import _slow_size_um

        section, header = self._header(1.5, 1.5, "4:1", 128)
        assert _slow_size_um(section, header, 128) == pytest.approx(0.375)

    @pytest.mark.parametrize('aspect,expected', [
        ("1:1", 2.0), ("2:1", 1.0), ("4:1", 0.5), ("8:1", 0.25),
    ])
    def test_every_aspect_ratio_in_the_archive(self, aspect, expected):
        from gratinglab.metrology.io.spm import _slow_size_um

        section, header = self._header(2, 2, aspect, 128)
        assert _slow_size_um(section, header, 128) == pytest.approx(expected)

    def test_the_plane_sections_field_of_the_same_name_is_ignored(self):
        """
        The trap, pinned.

        Every other field here is read from the plane section in preference to
        the file-level header, so writing this one that way is the obvious
        mistake. The section's `\\Aspect Ratio` is a bare number -- 0.352542
        on the TASTE file -- and is not a scan shape. Reading it turns that
        file's 1.57 µm slow axis into 4.44 µm.
        """
        from gratinglab.metrology.io.spm import _aspect_ratio, _slow_size_um

        section, header = self._header(2, 2, "1:1", 512, plane_aspect="0.352542")
        assert _aspect_ratio(header) == 1.0
        assert _slow_size_um(section, header, 512) == 2.0

    def test_a_bare_number_is_not_taken_as_a_ratio(self):
        from gratinglab.metrology.io.spm import _aspect_ratio

        assert _aspect_ratio("\\Aspect Ratio: 1.75489\n") == 1.0

    def test_an_absent_or_unusable_field_assumes_square(self):
        """1.0 is what this reader assumed unconditionally before."""
        from gratinglab.metrology.io.spm import _aspect_ratio

        assert _aspect_ratio("\\Lines: 512\n") == 1.0
        assert _aspect_ratio("\\Aspect Ratio: 0:1\n") == 1.0
        assert _aspect_ratio("\\Aspect Ratio: 4:0\n") == 1.0

    def test_the_aspect_ratio_composes_with_a_scan_stopped_early(self):
        """
        Both corrections apply, and they are independent.

        A 4:1 scan set to 512 lines that recorded 401 covers 2/4 µm of slow
        axis over the full scan and 401/512 of that in the file.
        """
        from gratinglab.metrology.io.spm import _slow_size_um

        section, header = self._header(2, 2, "4:1", 512)
        assert _slow_size_um(section, header, 401) == pytest.approx(
            0.5 * 401 / 512)

    def test_against_gwyddions_own_numbers(self):
        """
        Every .spm in the archive with a Gwyddion export beside it.

        Gwyddion states the slow extent as `# Height:` in its export, which is
        an independent reading of the same instrument file. Before this fix 3
        of 70 pairs agreed -- exactly the 1:1 ones; after it, 70 of 70.

        Skips unless the archive is reachable. Point `GRATINGLAB_AFM_DIR` at a
        tree of dated run directories to run it over all of them.
        """
        import re

        from gratinglab.metrology.io.spm import scan_y_size_um
        from scans import scan_dir

        directory = scan_dir()
        if directory is None:
            pytest.skip("AFM scans not found; set GRATINGLAB_AFM_DIR")

        checked = []
        for spm in sorted(directory.rglob('*.spm')):
            export = spm.with_suffix('.txt')
            if not export.exists():
                continue
            stated = re.search(r"#\s*Height:\s*([\d.]+)\s*(µm|um|nm)",
                               export.read_text(errors='replace')[:400])
            if stated is None:
                continue
            gwyddion = float(stated.group(1))
            if stated.group(2) == 'nm':
                gwyddion /= 1000.0
            ours = scan_y_size_um(spm)
            assert ours is not None, spm.name
            checked.append((spm.name, ours, gwyddion))

        if not checked:
            pytest.skip(f"no .spm/.txt pairs under {directory}")

        wrong = [(n, o, g) for n, o, g in checked if abs(o - g) / g >= 0.02]
        assert not wrong, (
            f"{len(wrong)} of {len(checked)} pairs disagree with Gwyddion's "
            f"stated height, e.g. {wrong[:3]}")
