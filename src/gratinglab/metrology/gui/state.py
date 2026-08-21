"""
Form state, and the rules for turning it into settings.

Pure: no Qt import, no file access, no plotting. Everything a widget needs to
decide *what to compute* is decided here, so it can be tested without a display.

The window reads values off its controls into a `FormState`, calls `build`, and
either gets an `AnalysisSettings` or a list of field errors it can attach to the
offending inputs.
"""
from __future__ import annotations

from dataclasses import dataclass

from ..settings import (AnalysisSettings, MAX_FACET_TRIM, VALID_BLAZE_SIDES)


@dataclass(frozen=True, slots=True)
class FieldError:
    """A problem with one named input"""
    field: str
    message: str


@dataclass(frozen=True, slots=True)
class FormState:
    """
    Raw values as the controls hold them.

    Numbers are already typed here rather than being strings, because Qt spin
    boxes hand back floats and ints. What this layer adds is the *semantic*
    checking the widgets cannot express - a spin box can cap a range, but it
    cannot know that facet trim above 0.28 leaves no facet to fit.
    """
    period_est: float = 315.0
    facet_trim: float = 0.1
    blaze_side: str = 'negative_slope'
    edge_exclusion_periods: float = 0.6
    use_row_groups: bool = True
    n_row_groups: int = 20
    scan_x_size: float = 2.0
    # Only meaningful for raw Nanoscope input; a text export has one plane.
    spm_direction: str = 'Retrace'

    # Owned by the Import tab. Image flattening acts on the 2-D scan; profile
    # flattening acts on each averaged row group and is the one that moves the
    # measured angle.
    image_flatten_method: str = 'align_rows'
    flatten_method: str = 'level_grooves'
    flatten_poly_order: int = 2
    flatten_feature: str = 'peaks'
    flatten_exclude_edges: float = 0.05

    @classmethod
    def from_settings(cls, settings: AnalysisSettings) -> "FormState":
        """Initial control values, taken from config.py defaults"""
        return cls(
            period_est=settings.period_est,
            facet_trim=settings.facet_trim,
            blaze_side=settings.blaze_side,
            edge_exclusion_periods=settings.edge_exclusion_periods,
            use_row_groups=settings.use_row_groups,
            n_row_groups=settings.n_row_groups,
            scan_x_size=settings.scan_x_size,
            spm_direction=settings.spm_direction,
            image_flatten_method=settings.image_flatten_method,
            flatten_method=settings.flatten_method,
            flatten_poly_order=settings.flatten_poly_order,
            flatten_feature=settings.flatten_feature,
            flatten_exclude_edges=settings.flatten_exclude_edges,
        )


def validate(form: FormState) -> tuple[FieldError, ...]:
    """
    Problems that would make the analysis fail or mislead.

    Delegates to AnalysisSettings.validate so the GUI and any other caller apply
    one set of rules, then adds the messages worth phrasing for someone looking
    at a window rather than a traceback.
    """
    settings = _to_settings(form)
    errors = [FieldError(field, message) for field, message in settings.validate()]

    # Phrased for the window: the analysis would return None and the user would
    # see "no measurements" with no hint that one control caused it.
    if form.facet_trim > MAX_FACET_TRIM:
        errors = [e for e in errors if e.field != 'facet_trim']
        errors.append(FieldError(
            'facet_trim',
            f"Above {MAX_FACET_TRIM} the trim removes the whole facet and no "
            f"angles can be measured at all."))

    if form.blaze_side not in VALID_BLAZE_SIDES:
        pass  # already reported by settings.validate, with the valid list

    return tuple(errors)


def _to_settings(form: FormState) -> AnalysisSettings:
    """Overlay the form onto config.py defaults, without validating"""
    return AnalysisSettings.from_config().with_(
        period_est=form.period_est,
        facet_trim=form.facet_trim,
        blaze_side=form.blaze_side,
        edge_exclusion_periods=form.edge_exclusion_periods,
        use_row_groups=form.use_row_groups,
        n_row_groups=form.n_row_groups,
        scan_x_size=form.scan_x_size,
        spm_direction=form.spm_direction,
        image_flatten_method=form.image_flatten_method,
        flatten_method=form.flatten_method,
        flatten_poly_order=form.flatten_poly_order,
        flatten_feature=form.flatten_feature,
        flatten_exclude_edges=form.flatten_exclude_edges,
    )


def build(form: FormState) -> tuple[AnalysisSettings | None, tuple[FieldError, ...]]:
    """
    Turn a form into settings.

    Returns (settings, errors). Settings is None when errors is non-empty, so a
    caller cannot accidentally run an analysis on invalid input by ignoring the
    second element.
    """
    errors = validate(form)
    if errors:
        return None, errors
    return _to_settings(form), ()


def summarize_result(result) -> str:
    """
    The results panel text.

    Pure so its wording is testable, and because formatting numbers is exactly
    the kind of thing that quietly goes wrong (a missing key, a nan) in a place
    where the only symptom is a blank widget.
    """
    if result is None:
        return "No measurements."

    r2 = [q.get('blaze_r2') for q in result.get('quality', [])]
    r2 = [v for v in r2 if v is not None]
    worst_r2 = f"{min(r2):.4f}" if r2 else "n/a"

    n_groups = result.get('n_groups')
    mode = f"row groups x{n_groups}" if n_groups else "averaged profile"

    icc = result.get('icc')
    sem = result.get('sem_corrected', result.get('sem', float('nan')))
    icc_note = (f"   ICC {icc:.2f}, N_eff {result.get('n_effective', float('nan')):.0f}"
                if icc is not None and icc == icc else "")

    return (
        f"Mean blaze angle : {result['mean_angle']:.3f}deg  "
        f"+/- {sem:.3f}deg (SEM){icc_note}\n"
        f"Spread           : sigma = {result['std_angle']:.3f}deg   "
        f"range {result['min_angle']:.2f}-{result['max_angle']:.2f}deg\n"
        f"Measurements     : N = {result['n_grooves']}   ({mode})\n"
        f"Period           : {result['period_nm']:.2f} "
        f"+/- {result.get('period_std', 0):.2f} nm   worst fit R2 = {worst_r2}"
    )
