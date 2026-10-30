"""
Scar repair and the region mask: the two steps that used to send a scan to Gwyddion.

The property that matters is recovery. The synthetic scan is a known grating, so
the tests damage it in the two ways these stages exist for, seven scarred scan
lines and a particle, and assert two things each time: that the damage really
does bias the answer when left alone, and that the repaired or masked scan gives
the clean scan's answer. A test of only the second half would pass on damage too
small to matter, and so would prove nothing.
"""
import contextlib
import io

import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg')

from scans import SYNTHETIC, real_scan                                  # noqa: E402

from gratinglab.metrology.analyzer import analyze_single_file           # noqa: E402
from gratinglab.metrology.boundary import build_boundary_profile       # noqa: E402
from gratinglab.metrology.core.image_flatten import flatten_image      # noqa: E402
from gratinglab.metrology.core.processing import (                     # noqa: E402
    load_afm_data, raw_data, raw_data_multi_group)
from gratinglab.metrology.core.region_mask import (                    # noqa: E402
    apply_region_mask, bridge_nan, describe_mask, pixels_from_drag,
    region_mask, reject_grooves_near_nan)
from gratinglab.metrology.core.scars import (                          # noqa: E402
    mark_scars, repair_scars, repair_scars_for)
from gratinglab.metrology.io.file_io import save_results_to_file       # noqa: E402
from gratinglab.metrology.io.metrics import write_profile_metrics      # noqa: E402
from gratinglab.metrology.settings import AnalysisSettings             # noqa: E402

#: Seven scarred lines: two pairs and a triple, so bands of every width up to
#: three, plus two single lines. Each is raised by 60 nm over a third of its
#: length: over half, which part of the line is the scar is a matter of opinion.
SCAR_LINES = (10, 11, 30, 45, 46, 47, 58)
SCAR_HEIGHT_M = 60e-9

#: A 150 nm Gaussian particle, and the rectangle that covers it.
PARTICLE = (20, 44, 100, 160)


def _quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


@pytest.fixture(scope='module')
def clean():
    data, size = _quiet(load_afm_data, str(SYNTHETIC), default_scan_size=2.0)
    return data, size


@pytest.fixture(scope='module')
def scarred(clean):
    data, size = clean
    rng = np.random.default_rng(1)
    out = data.copy()
    truth = np.zeros(data.shape, dtype=bool)
    for line in SCAR_LINES:
        start = rng.integers(0, data.shape[1] // 3)
        stop = start + data.shape[1] // 3
        out[line, start:stop] += SCAR_HEIGHT_M
        truth[line, start:stop] = True
    return out, truth


@pytest.fixture(scope='module')
def particle(clean):
    data, _ = clean
    r0, r1, c0, c1 = PARTICLE
    out = data.copy()
    yy, xx = np.mgrid[r0:r1, c0:c1]
    out[r0:r1, c0:c1] += 150e-9 * np.exp(
        -(((yy - 32) / 8) ** 2 + ((xx - 130) / 18) ** 2))
    return out


def _write(path, data):
    """Four header lines then tab-separated metres, as the fixture is written."""
    lines = ["# Channel: ZSensor", "# Width: 2.000 µm", "# Height: 0.333 µm",
             "# Value units: m"]
    lines += ["\t".join(f"{v:.5e}" for v in row) for row in data]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


@pytest.fixture(scope='module')
def files(tmp_path_factory, clean, scarred, particle):
    root = tmp_path_factory.mktemp('damaged')
    # The clean scan is written through the same formatter as the damaged ones,
    # so the comparison sees only the damage and not a rounding difference.
    return {'clean': _write(root / 'clean.txt', clean[0]),
            'scarred': _write(root / 'scarred.txt', scarred[0]),
            'particle': _write(root / 'particle.txt', particle)}


BASE = AnalysisSettings()


def _angle(path, settings=BASE):
    result = _quiet(analyze_single_file, path, show_plots=False,
                    settings=settings)
    assert result is not None
    return result


def _groove(data, size, settings=BASE):
    """The averaged boundary groove, prepared the way the export prepares it."""
    data, _ = apply_region_mask(data, settings)
    data = flatten_image(data, settings.image_flatten_method)
    data, _ = repair_scars_for(data, settings)
    return build_boundary_profile(data, size, settings)


# ── Scar detection ───────────────────────────────────────────────────────────

class TestMarking:
    def test_the_clean_scan_has_no_scars(self, clean):
        """Noise alone must never be marked: the length rule is what ensures it."""
        assert not mark_scars(flatten_image(clean[0], 'align_rows')).any()

    def test_every_scarred_line_is_found_and_no_other(self, scarred):
        marks = mark_scars(flatten_image(scarred[0], 'align_rows'))
        assert tuple(np.flatnonzero(marks.any(axis=1))) == SCAR_LINES

    def test_every_injected_pixel_is_covered(self, scarred):
        """
        Including where two adjacent scars overlap for only a few pixels, the
        case a length test per band width used to cut out of the run.
        """
        data, truth = scarred
        marks = mark_scars(flatten_image(data, 'align_rows'))
        assert not (truth & ~marks).any()

    def test_a_negative_scar_is_found_too(self, clean):
        data = clean[0].copy()
        data[25, 50:250] -= SCAR_HEIGHT_M
        marks = mark_scars(flatten_image(data, 'align_rows'))
        assert tuple(np.flatnonzero(marks.any(axis=1))) == (25,)

    def test_a_short_excursion_is_not_a_scar(self, clean):
        """Added after levelling, to test the length rule alone."""
        data = flatten_image(clean[0], 'align_rows')
        data[25, 50:60] += SCAR_HEIGHT_M          # ten pixels, under min_length
        assert not mark_scars(data).any()

    def test_what_levelling_did_to_a_line_is_not_a_scar(self, clean):
        """
        Ten high pixels move their line's median by about 2 nm (a grating's
        heights spread over its whole depth), so align_rows shifts the whole
        line. The detector levels lines by their median step instead, which
        the spike cannot move, so only the spike is left, and it is short.
        """
        data = clean[0].copy()
        data[25, 50:60] += SCAR_HEIGHT_M
        assert not mark_scars(flatten_image(data, 'align_rows')).any()

    def test_a_line_wider_than_max_width_is_not_a_scar(self, clean):
        data = clean[0].copy()
        data[20:26, 50:250] += SCAR_HEIGHT_M      # six lines
        assert not mark_scars(flatten_image(data, 'align_rows'),
                              max_width=4).any()


class TestRepair:
    def test_interpolation_restores_the_clean_image(self, clean, scarred):
        """To within the noise of interpolating between the lines either side."""
        reference = flatten_image(clean[0], 'align_rows')
        repaired, report = repair_scars(flatten_image(scarred[0], 'align_rows'),
                                        'interpolate')
        assert report.n_lines == len(SCAR_LINES)
        assert report.lines == SCAR_LINES
        assert np.max(np.abs(repaired - reference)) < 1e-9

    def test_lines_without_scars_are_untouched(self, clean, scarred):
        before = flatten_image(scarred[0], 'align_rows')
        repaired, _ = repair_scars(before, 'interpolate')
        others = [r for r in range(before.shape[0]) if r not in SCAR_LINES]
        assert np.array_equal(repaired[others], before[others])

    def test_drop_blanks_exactly_the_scarred_lines(self, scarred):
        dropped, report = repair_scars(flatten_image(scarred[0], 'align_rows'),
                                       'drop')
        blank = np.flatnonzero(np.all(np.isnan(dropped), axis=1))
        assert tuple(blank) == SCAR_LINES == report.lines
        assert np.isfinite(np.delete(dropped, SCAR_LINES, axis=0)).all()

    def test_none_is_the_identity(self, scarred):
        out, report = repair_scars(scarred[0], 'none')
        assert np.array_equal(out, scarred[0])
        assert report.n_lines == 0

    def test_the_input_is_never_modified(self, scarred):
        data = flatten_image(scarred[0], 'align_rows')
        copy = data.copy()
        repair_scars(data, 'interpolate')
        repair_scars(data, 'drop')
        assert np.array_equal(data, copy)

    def test_the_summary_says_what_was_done(self, scarred):
        _, report = repair_scars(flatten_image(scarred[0], 'align_rows'),
                                 'interpolate')
        assert report.summary.startswith("interpolate: 7 line(s) interpolated")

    def test_an_unknown_method_is_refused(self, clean):
        with pytest.raises(ValueError, match="unknown scar repair"):
            repair_scars(clean[0], 'gwyddion')


# ── The region mask ──────────────────────────────────────────────────────────

class TestMaskShapes:
    def test_the_three_shapes_combine(self):
        settings = BASE.with_(mask_rect=(1, 3, 2, 4), mask_rows=((6, 7),),
                              mask_cols=((9, 10),))
        mask = region_mask((8, 12), settings)
        assert mask[1:3, 2:4].all() and mask[6].all() and mask[:, 9].all()
        assert mask.sum() == 4 + 12 + 8 - 1

    def test_ranges_past_the_image_are_clipped(self):
        mask = region_mask((5, 5), BASE.with_(mask_rect=(3, 99, 3, 99)))
        assert mask.sum() == 4

    def test_no_mask_returns_the_same_array(self, clean):
        """So an unmasked run is bit-identical to the pipeline before masks."""
        out, n = apply_region_mask(clean[0], BASE)
        assert out is clean[0] and n == 0

    def test_description(self):
        settings = BASE.with_(mask_rect=(1, 3, 2, 4), mask_cols=((9, 10),))
        assert describe_mask(settings) == "rect rows 1-3 cols 2-4; cols 9-10"
        assert describe_mask(BASE) == ""

    def test_a_drag_maps_to_the_pixels_it_touched(self):
        # 10 lines, 20 columns over 2 um: a pixel is 0.1 um wide.
        assert pixels_from_drag(0.25, 2.5, 0.55, 4.2, (10, 20), 2.0) == (2, 5, 2, 6)
        # Either direction of drag, and clipped at the far edges.
        assert pixels_from_drag(5.0, 99, -1, -3, (10, 20), 2.0) == (0, 10, 0, 20)


class TestRowGrouping:
    def test_masked_lines_join_no_group(self, clean):
        data, size = clean
        masked, _ = apply_region_mask(data, BASE.with_(mask_rows=((0, 24),)))
        _, profiles, info = _quiet(raw_data_multi_group, masked, size, 10)
        assert info['rows_excluded'] == 24
        assert info['group_ranges'][0][0] == 24
        assert all(np.isfinite(p).all() for p in profiles)

    def test_groups_are_cut_from_what_is_left(self, clean):
        data, size = clean
        masked, _ = apply_region_mask(data, BASE.with_(mask_rows=((0, 24),)))
        _, kept, _ = _quiet(raw_data_multi_group, data[24:], size, 10)
        _, profiles, _ = _quiet(raw_data_multi_group, masked, size, 10)
        for a, b in zip(profiles, kept):
            np.testing.assert_array_equal(a, b)

    def test_the_first_band_starts_below_the_lines_left_out_at_the_top(self, clean):
        """The bands are centred, so the top band is not the top of the image.

        64 lines in 10 bands of 6 leaves 4 over; two come off each end.
        """
        data, size = clean
        _, profiles, info = _quiet(raw_data_multi_group, data, size, 10)
        assert info['rows_excluded'] == 0
        assert info['rows_unbanded'] == 4
        rows = info['rows_per_group']
        first = data[2:2 + rows]
        np.testing.assert_array_equal(
            profiles[0], 1e9 * (first.mean(axis=0) - first.mean(axis=0).min()))

    def test_the_bands_are_the_same_set_read_either_way(self, clean):
        """Reversing the slow axis must not re-cut the bands.

        It is not a hypothetical: a Nanoscope `.spm` and the Gwyddion export of
        the same scan disagree about which end of the slow axis comes first, so
        the same image reaches this function both ways up. While the leftover
        lines went to the last band, that moved every boundary.
        """
        data, size = clean
        _, upright, _ = _quiet(raw_data_multi_group, data, size, 20)
        _, flipped, _ = _quiet(raw_data_multi_group, data[::-1], size, 20)
        assert len(upright) == len(flipped) == 20
        for i, (a, b) in enumerate(zip(upright, reversed(flipped))):
            np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12,
                                       err_msg=f"band {i}")

    def test_a_fully_masked_column_is_a_gap(self, clean):
        data, size = clean
        masked, _ = apply_region_mask(data, BASE.with_(mask_cols=((100, 160),)))
        _, profile = raw_data(masked, size)
        assert np.isnan(profile[100:160]).all()
        assert np.isfinite(np.delete(profile, range(100, 160))).all()

    def test_grooves_reaching_into_a_gap_are_rejected(self):
        profile = np.zeros(100)
        profile[40:50] = np.nan
        kept, n = reject_grooves_near_nan([10, 35, 45, 55, 80], profile, 6)
        assert list(kept) == [10, 80] and n == 3

    def test_bridging_is_a_straight_line(self):
        y = np.array([0.0, np.nan, np.nan, 3.0])
        np.testing.assert_allclose(bridge_nan(y), [0, 1, 2, 3])


# ── Recovery: the point of both stages ───────────────────────────────────────

class TestScarRecovery:
    def test_scars_bias_the_blaze_angle_and_repair_removes_it(self, files):
        clean = _angle(files['clean'])['mean_angle']
        raw = _angle(files['scarred'])['mean_angle']
        assert abs(raw - clean) > 0.02, "the scars are too weak to test anything"
        for method in ('interpolate', 'drop'):
            repaired = _angle(files['scarred'], BASE.with_(scar_repair=method))
            assert abs(repaired['mean_angle'] - clean) < 0.005, method
            assert repaired['scar_lines'] == len(SCAR_LINES)

    def test_scars_bias_the_averaged_groove_and_repair_removes_it(self, clean,
                                                                  scarred):
        reference = _groove(*clean).y_avg_nm
        raw = _groove(scarred[0], clean[1]).y_avg_nm
        assert np.max(np.abs(raw - reference)) > 0.5
        for method in ('interpolate', 'drop'):
            repaired = _groove(scarred[0], clean[1],
                               BASE.with_(scar_repair=method)).y_avg_nm
            assert np.max(np.abs(repaired - reference)) < 0.1, method

    def test_repair_on_a_clean_scan_changes_nothing(self, files):
        a = _angle(files['clean'])
        b = _angle(files['clean'], BASE.with_(scar_repair='interpolate'))
        assert a['mean_angle'] == b['mean_angle'] and b['scar_lines'] == 0

    def test_dropped_lines_are_counted_out_of_row_grouping(self, files):
        result = _angle(files['scarred'], BASE.with_(scar_repair='drop'))
        assert result['rows_excluded'] == len(SCAR_LINES)


class TestMaskRecovery:
    def test_a_particle_biases_the_blaze_angle(self, files):
        clean = _angle(files['clean'])['mean_angle']
        # Unmasked, the particle is fitted as a facet, and the fit says so.
        with pytest.warns(UserWarning, match='unusually large'):
            damaged = _angle(files['particle'])['mean_angle']
        assert abs(damaged - clean) > 0.3

    @pytest.mark.parametrize('mask', [
        {'mask_rect': PARTICLE},
        {'mask_rows': ((PARTICLE[0], PARTICLE[1]),)},
        {'mask_cols': ((PARTICLE[2], PARTICLE[3]),)},
    ], ids=['rect', 'rows', 'cols'])
    def test_masked_the_particle_has_no_influence_left(self, files, mask):
        """Exactly the clean scan's answer under the same mask."""
        settings = BASE.with_(**mask)
        damaged = _angle(files['particle'], settings)
        reference = _angle(files['clean'], settings)
        assert damaged['mean_angle'] == pytest.approx(reference['mean_angle'],
                                                      abs=1e-9)
        assert damaged['masked_pixels'] > 0
        assert damaged['region_mask'] == describe_mask(settings)

    @pytest.mark.parametrize('mask', [
        {'mask_rect': PARTICLE},
        {'mask_rows': ((PARTICLE[0], PARTICLE[1]),)},
    ], ids=['rect', 'rows'])
    def test_masked_the_answer_is_the_clean_one(self, files, mask):
        """
        Against the *unmasked* clean scan, within the sampling difference of
        measuring fewer grooves. Columns are left out of this one: masking 60 of
        384 columns keeps 40 measurements of 100, too few for a tight bound.
        """
        clean = _angle(files['clean'])['mean_angle']
        masked = _angle(files['particle'], BASE.with_(**mask))['mean_angle']
        assert abs(masked - clean) < 0.02

    def test_masked_columns_contribute_no_grooves(self, files):
        r0, r1, c0, c1 = PARTICLE
        result = _angle(files['particle'],
                        BASE.with_(mask_cols=((c0, c1),),
                                   use_row_groups=False))
        x = np.asarray(result['raw_x'])
        positions = x[np.asarray(result['groove_centers'])]
        assert not ((positions >= x[c0]) & (positions < x[c1])).any()
        assert result['mask_rejected_grooves'] > 0

    def test_the_particle_biases_the_averaged_groove_and_the_mask_removes_it(
            self, clean, particle):
        reference = _groove(*clean)
        raw = _groove(particle, clean[1])
        masked = _groove(particle, clean[1], BASE.with_(mask_rect=PARTICLE))
        assert np.max(np.abs(raw.y_avg_nm - reference.y_avg_nm)) > 2.0
        assert np.max(np.abs(masked.y_avg_nm - reference.y_avg_nm)) < 0.2

    def test_masked_columns_leave_the_averaged_groove_untouched(self, clean,
                                                                particle):
        settings = BASE.with_(mask_cols=((PARTICLE[2], PARTICLE[3]),))
        masked = _groove(particle, clean[1], settings)
        reference = _groove(clean[0], clean[1], settings)
        np.testing.assert_allclose(masked.y_avg_nm, reference.y_avg_nm,
                                   atol=1e-9)
        assert masked.n_mask_rejected > 0
        assert masked.metrics['mask_rejected_grooves'] == masked.n_mask_rejected


# ── Settings and provenance ──────────────────────────────────────────────────

class TestSettings:
    def test_config_lists_become_tuples(self, monkeypatch):
        from gratinglab.metrology import config
        monkeypatch.setattr(config, 'MASK_RECT', [1, 2, 3, 4], raising=False)
        monkeypatch.setattr(config, 'MASK_ROWS', [[5, 6]], raising=False)
        monkeypatch.setattr(config, 'SCAR_REPAIR', 'drop', raising=False)
        settings = AnalysisSettings.from_config()
        assert settings.mask_rect == (1, 2, 3, 4)
        assert settings.mask_rows == ((5, 6),)
        assert settings.scar_repair == 'drop'
        hash(settings)        # frozen and hashable, like every other setting

    @pytest.mark.parametrize('change, field', [
        ({'scar_repair': 'gwyddion'}, 'scar_repair'),
        ({'scar_repair': 'interpolate', 'scar_threshold': 0}, 'scar_threshold'),
        ({'mask_rect': (5, 5, 0, 3)}, 'mask_rect'),
        ({'mask_rect': (-1, 5, 0, 3)}, 'mask_rect'),
        ({'mask_rows': ((4, 2),)}, 'mask_rows'),
        ({'mask_cols': ((0, 0),)}, 'mask_cols'),
    ])
    def test_bad_values_are_reported(self, change, field):
        fields = [f for f, _ in BASE.with_(**change).validate()]
        assert field in fields

    def test_the_defaults_are_valid(self):
        assert BASE.validate() == ()


class TestProvenance:
    def test_the_summary_file_records_both(self, files, tmp_path):
        settings = BASE.with_(scar_repair='interpolate', mask_rect=PARTICLE)
        result = _angle(files['scarred'], settings)
        _quiet(save_results_to_file, [result], output_dir=str(tmp_path))
        text = next(tmp_path.glob('analysis_summary_*.txt')).read_text(
            encoding='utf-8')
        assert "Scar repair: interpolate, 7 line(s)" in text
        assert "Region mask: rect rows 20-44 cols 100-160" in text

    def test_an_untouched_run_records_neither(self, files, tmp_path):
        _quiet(save_results_to_file, [_angle(files['clean'])],
               output_dir=str(tmp_path))
        text = next(tmp_path.glob('analysis_summary_*.txt')).read_text(
            encoding='utf-8')
        assert "Scar repair" not in text and "Region mask" not in text

    def test_the_boundary_metrics_record_the_mask(self, clean, particle,
                                                  tmp_path):
        profile = _groove(particle, clean[1], BASE.with_(mask_rect=PARTICLE))
        path = tmp_path / 'metrics.txt'
        write_profile_metrics(path, profile.metrics)
        assert "Region mask: rect rows 20-44 cols 100-160" in path.read_text(
            encoding='utf-8')

    def test_the_export_workflow_records_scar_repair(self, files, tmp_path,
                                                     monkeypatch):
        from gratinglab.metrology import workflows
        monkeypatch.setattr(workflows, 'RESULTS_DIR', str(tmp_path))
        out = _quiet(workflows.run_boundary_profile_export,
                     BASE.with_(scar_repair='interpolate'), files['scarred'])
        assert out['metrics']['scar_lines'] == len(SCAR_LINES)
        text = next(tmp_path.glob('groove_analysis_metrics_*.txt')).read_text(
            encoding='utf-8')
        assert "Scar repair: interpolate (7 line(s))" in text


# ── A real scan ──────────────────────────────────────────────────────────────

def test_a_real_scan_is_not_riddled_with_false_scars():
    """
    The scan with the worst row-offset spread in the dataset. Whatever scars it
    really has, the default threshold must not mark a large share of its lines,
    or the repair would be rewriting the measurement rather than mending it.
    """
    data, _ = _quiet(load_afm_data, real_scan('20250820_280C_00004.txt'),
                     default_scan_size=2.0)
    marks = mark_scars(flatten_image(data, 'align_rows'))
    assert marks.any(axis=1).mean() < 0.1
