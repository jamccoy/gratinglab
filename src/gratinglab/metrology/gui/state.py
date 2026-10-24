"""
Form state, and the rules for turning it into settings.

The module is pure, with no Qt import, no file access and no plotting.
Everything a widget needs to decide what to compute is decided here, so it can
be tested without a display.

The window reads values off its controls into a `FormState`, calls `build`, and
either gets an `AnalysisSettings` or a list of field errors it can attach to the
offending inputs.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from endstation.forms import FieldError as _FieldError

from ..settings import (AnalysisSettings, MAX_FACET_TRIM, VALID_BLAZE_SIDES)


#: A problem with one named input. The same type the main window uses -- this
#: package carried an identical definition of its own until both became
#: :class:`endstation.forms.FieldError`.
#:
#: The two `build()` contracts still differ: the main window's
#: raises, this one returns `(settings, errors)`. Only the duplicated type went.
FieldError = _FieldError


@dataclass(frozen=True, slots=True)
class FormState:
    """
    Raw values as the controls hold them.

    Numbers are already typed here instead of being strings, because Qt spin
    boxes hand back floats and ints. What this layer adds is the semantic
    checking the widgets cannot express: a spin box can cap a range, but it
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
    # Groove orientation, also owned by the Import tab: measured on the 2-D
    # scan after image flattening and corrected when large enough to matter.
    groove_orientation: str = 'auto'
    flatten_method: str = 'level_grooves'
    flatten_poly_order: int = 2
    flatten_feature: str = 'peaks'
    flatten_exclude_edges: float = 0.05

    # Also owned by the Import tab: scar repair beside image flattening, and
    # the region mask drawn on the image, in pixels (see core/region_mask.py).
    scar_repair: str = 'none'
    scar_threshold: float = 4.0
    mask_rect: Optional[tuple] = None
    mask_rows: tuple = ()
    mask_cols: tuple = ()

    # Tip deconvolution, applied to the 2-D scan before either the blaze-angle
    # analysis or the boundary profile sees it. Off by default: a corrected and
    # an uncorrected depth are different measurements.
    tip_correction: str = 'none'
    tip_radius_nm: float = 1.0
    tip_half_angle_deg: float = 18.0
    tip_half_angle_right_deg: Optional[float] = None

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
            groove_orientation=settings.groove_orientation,
            flatten_method=settings.flatten_method,
            flatten_poly_order=settings.flatten_poly_order,
            flatten_feature=settings.flatten_feature,
            flatten_exclude_edges=settings.flatten_exclude_edges,
            scar_repair=settings.scar_repair,
            scar_threshold=settings.scar_threshold,
            mask_rect=settings.mask_rect,
            mask_rows=settings.mask_rows,
            mask_cols=settings.mask_cols,
            tip_correction=settings.tip_correction,
            tip_radius_nm=settings.tip_radius_nm,
            tip_half_angle_deg=settings.tip_half_angle_deg,
            tip_half_angle_right_deg=settings.tip_half_angle_right_deg,
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
        groove_orientation=form.groove_orientation,
        flatten_method=form.flatten_method,
        flatten_poly_order=form.flatten_poly_order,
        flatten_feature=form.flatten_feature,
        flatten_exclude_edges=form.flatten_exclude_edges,
        scar_repair=form.scar_repair,
        scar_threshold=form.scar_threshold,
        mask_rect=form.mask_rect,
        mask_rows=tuple(form.mask_rows),
        mask_cols=tuple(form.mask_cols),
        tip_correction=form.tip_correction,
        tip_radius_nm=form.tip_radius_nm,
        tip_half_angle_deg=form.tip_half_angle_deg,
        tip_half_angle_right_deg=form.tip_half_angle_right_deg,
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

    The function is pure so that its wording is testable, and because number
    formatting can fail quietly (a missing key, a nan) in a place where the
    only symptom is a blank widget.
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

    # What was taken out of the scan before it was measured, only when
    # something was: a repaired or masked scan is a different measurement.
    cleanup = ""
    if result.get('scar_repair', 'none') != 'none':
        cleanup += (f"\nScar repair      : {result['scar_repair']}, "
                    f"{result.get('scar_lines', 0)} line(s)")
    if result.get('region_mask'):
        cleanup += (f"\nRegion mask      : {result['region_mask']}   "
                    f"{result.get('mask_rejected_grooves', 0)} groove(s) "
                    f"rejected")

    return (
        f"Mean blaze angle : {result['mean_angle']:.3f}deg  "
        f"+/- {sem:.3f}deg (SEM){icc_note}\n"
        f"Spread           : sigma = {result['std_angle']:.3f}deg   "
        f"range {result['min_angle']:.2f}-{result['max_angle']:.2f}deg\n"
        f"Measurements     : N = {result['n_grooves']}   ({mode})\n"
        f"Period           : {result['period_nm']:.2f} "
        f"+/- {result.get('period_std', 0):.2f} nm   worst fit R2 = {worst_r2}"
        + orientation_line(result.get('groove_orientation'))
        + cleanup
    )


def orientation_line(report) -> str:
    """
    The groove-rotation line of the results panel, or nothing.

    Takes the report as the plain dict a result carries. Stated whether or not
    the correction was applied: a period read along a rotated scan line is
    long by 1/cos(theta), and the panel is where a reader would notice.
    """
    if not report:
        return ""
    angle = report.get('angle_deg')
    if angle is None:
        return f"\nGroove rotation  : not measured ({report.get('reason', '')})"
    done = ("corrected" if report.get('applied')
            else f"not corrected, {report.get('reason', '')}")
    return (f"\nGroove rotation  : {angle:+.3f} +/- "
            f"{report.get('stderr_deg', float('nan')):.3f}deg   ({done})")
