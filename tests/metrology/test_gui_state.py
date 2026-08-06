"""
The pure half of the GUI: form validation and result formatting.

No Qt import anywhere in this file, deliberately - these rules are the ones worth
testing, and they should be testable without a display.
"""
import math

from afm_analysis.gui.state import (FormState, build, validate, summarize_result)
from afm_analysis.settings import AnalysisSettings, MAX_FACET_TRIM


def test_defaults_are_valid():
    form = FormState.from_settings(AnalysisSettings.from_config())
    settings, errors = build(form)
    assert errors == ()
    assert settings is not None


def test_form_overrides_reach_the_settings():
    settings, errors = build(FormState(period_est=400.0, facet_trim=0.2,
                                       blaze_side='positive_slope',
                                       edge_exclusion_periods=1.0,
                                       use_row_groups=False, n_row_groups=8,
                                       scan_x_size=1.5))
    assert errors == ()
    assert settings.period_est == 400.0
    assert settings.facet_trim == 0.2
    assert settings.blaze_side == 'positive_slope'
    assert settings.use_row_groups is False
    assert settings.n_row_groups == 8


def test_excessive_facet_trim_is_rejected_with_a_useful_message():
    """
    The failure this prevents is silent: above ~0.286 the trim removes the whole
    facet, every fit fails, and the analysis returns None. A user sees "no
    measurements" with nothing pointing at the control responsible.
    """
    settings, errors = build(FormState(facet_trim=0.35))
    assert settings is None
    fields = [e.field for e in errors]
    assert 'facet_trim' in fields
    message = next(e.message for e in errors if e.field == 'facet_trim')
    assert str(MAX_FACET_TRIM) in message
    assert 'no angles' in message.lower() or 'whole facet' in message.lower()


def test_invalid_blaze_side_is_rejected():
    settings, errors = build(FormState(blaze_side='left'))
    assert settings is None
    assert 'blaze_side' in [e.field for e in errors]


def test_build_returns_no_settings_whenever_there_are_errors():
    """A caller that ignores the errors must not get usable settings"""
    for bad in (FormState(facet_trim=0.9), FormState(period_est=0.0),
                FormState(blaze_side='nonsense'), FormState(n_row_groups=1)):
        settings, errors = build(bad)
        assert errors, f"expected errors for {bad}"
        assert settings is None


def test_summarize_handles_a_missing_result():
    assert summarize_result(None) == "No measurements."


def test_summarize_survives_absent_optional_keys():
    """quality/sem/period_std are optional; a blank panel is a bad failure mode"""
    minimal = {
        'mean_angle': 30.0, 'std_angle': 1.0, 'min_angle': 28.0,
        'max_angle': 32.0, 'n_grooves': 10, 'period_nm': 315.0,
    }
    text = summarize_result(minimal)
    assert '30.000' in text
    assert 'n/a' in text          # no quality entries -> no R2 to report
    assert 'N = 10' in text
