"""
Analysis settings as a value, rather than as module state.

`config.py` remains the file a user edits, and is still the default source of
every value here. What this module adds is the ability to *pass* a different
configuration without mutating anything.

Before this existed, the only way to run an analysis with different parameters
was to rebind names inside `analyzer` itself:

    analyzer.FACET_TRIM = 0.2

which worked because `analyzer` binds its configuration at import time. The GUI
did exactly that. It is also why a CLI, a test, or two analyses with different
settings in one process were all impossible. Settings now travel as an argument.
"""
from dataclasses import dataclass, replace

from .core.image_flatten import VALID_IMAGE_FLATTEN_METHODS
from .core.tip import VALID_TIP_CORRECTIONS

# The blaze facet is trimmed 2.5x harder on the trough side than the land side
# (see core/analysis.py trim_facet), so the total removed is FACET_TRIM * 3.5.
# At 1/3.5 = 0.2857 nothing is left, every fit fails, and the analysis returns no
# measurements rather than degrading. Kept here so the GUI, any validator and the
# documentation all cite one number.
MAX_FACET_TRIM = 0.28

VALID_BLAZE_SIDES = ('negative_slope', 'positive_slope', 'longer')

#: Scan directions a Nanoscope file records. Retrace is the default because it is
#: the plane this project's existing Gwyddion exports were taken from, so a .spm
#: and its .txt export agree by construction rather than by luck.
VALID_SPM_DIRECTIONS = ('Retrace', 'Trace')
VALID_FLATTEN_METHODS = ('linear', 'polynomial', 'groove_peaks', 'level_grooves')
VALID_FLATTEN_FEATURES = ('peaks', 'troughs', 'both')


@dataclass(frozen=True)
class AnalysisSettings:
    """Every knob the blaze-angle analysis reads. Immutable by design."""

    # Scan
    scan_x_size: float = 2.0
    period_est: float = 315.0

    # Raw Nanoscope (.spm) input. Ignored for text exports, which carry one
    # channel and no direction.
    spm_channel: str = 'Height Sensor'
    spm_direction: str = 'Retrace'

    # Groove detection
    prominence_factor: float = 0.01
    distance_factor: float = 0.3
    edge_exclusion_periods: float = 0.6

    # Facet fitting
    facet_trim: float = 0.1
    blaze_side: str = 'negative_slope'

    # Image flattening: the 2-D stage, before rows are averaged. Affine methods
    # provably cannot change a blaze angle here (see core/image_flatten.py), so
    # the default costs nothing - measured 0.0000 deg on all eight samples.
    image_flatten_method: str = 'align_rows'

    # Tip correction: undoing what the tip shape did to the image, where that
    # is possible, on the 2-D array right after image flattening. Off by
    # default because it changes measured numbers (depth especially) and must
    # be asked for. The tip is a cone with a spherical apex cap; the radius is
    # the *apex radius* (a "2 nm wide" tip is radius 1), and the half angle is
    # measured from the tip axis. See core/tip.py for what can and cannot be
    # recovered.
    tip_correction: str = 'none'
    tip_radius_nm: float = 1.0
    tip_half_angle_deg: float = 18.0

    # Profile flattening: the 1-D stage, after each row group is averaged. This
    # one does move the answer - about 0.49 deg across the four methods.
    flatten_method: str = 'level_grooves'
    flatten_poly_order: int = 2
    flatten_exclude_edges: float = 0.05
    flatten_feature: str = 'peaks'

    # PCGrate boundary export. Only read by the boundary path; the blaze-angle
    # analysis ignores them.
    ggp_n_points: int = 2000
    ggp_apply_smoothing: bool = True
    ggp_smoothing_window: int = 5
    ggp_min_half_width: int = 10

    # Row groups
    use_row_groups: bool = True
    n_row_groups: int = 20

    # Diagnostic plots
    show_2d_image: bool = False
    show_individual_grooves: bool = False
    show_full_profile: bool = True
    show_flattening_diagnostic: bool = False
    show_local_angle_distribution: bool = True
    show_analyzed_regions: bool = False

    @classmethod
    def from_config(cls):
        """Build from the values currently in config.py.

        Read at call time rather than import time, so editing config.py and
        re-running in the same session behaves as a user would expect.
        """
        from . import config
        return cls(
            scan_x_size=config.SCAN_X_SIZE,
            period_est=config.PERIOD_EST,
            prominence_factor=config.PROMINENCE_FACTOR,
            distance_factor=config.DISTANCE_FACTOR,
            edge_exclusion_periods=config.EDGE_EXCLUSION_PERIODS,
            facet_trim=config.FACET_TRIM,
            blaze_side=config.BLAZE_SIDE,
            image_flatten_method=getattr(config, 'IMAGE_FLATTEN_METHOD',
                                         'align_rows'),
            tip_correction=getattr(config, 'TIP_CORRECTION', 'none'),
            tip_radius_nm=getattr(config, 'TIP_RADIUS_NM', 1.0),
            tip_half_angle_deg=getattr(config, 'TIP_HALF_ANGLE_DEG', 18.0),
            flatten_method=config.FLATTEN_METHOD,
            flatten_poly_order=config.FLATTEN_POLY_ORDER,
            flatten_exclude_edges=config.FLATTEN_EXCLUDE_EDGES,
            flatten_feature=config.FLATTEN_FEATURE,
            spm_channel=getattr(config, 'SPM_CHANNEL', 'Height Sensor'),
            spm_direction=getattr(config, 'SPM_DIRECTION', 'Retrace'),
            ggp_n_points=getattr(config, 'GGP_N_POINTS', 2000),
            ggp_apply_smoothing=getattr(config, 'GGP_APPLY_SMOOTHING', True),
            ggp_smoothing_window=getattr(config, 'GGP_SMOOTHING_WINDOW', 5),
            ggp_min_half_width=getattr(config, 'GGP_MIN_HALF_WIDTH', 10),
            use_row_groups=config.USE_ROW_GROUPS,
            n_row_groups=config.N_ROW_GROUPS,
            show_2d_image=config.SHOW_2D_IMAGE,
            show_individual_grooves=config.SHOW_INDIVIDUAL_GROOVES,
            show_full_profile=config.SHOW_FULL_PROFILE,
            show_flattening_diagnostic=config.SHOW_FLATTENING_DIAGNOSTIC,
            show_local_angle_distribution=config.SHOW_LOCAL_ANGLE_DISTRIBUTION,
            show_analyzed_regions=config.SHOW_ANALYZED_REGIONS,
        )

    def with_(self, **changes):
        """A copy with some fields replaced - the frozen equivalent of setattr"""
        return replace(self, **changes)

    def validate(self):
        """
        Problems that would make an analysis fail or mislead.

        Returns a tuple of (field_name, message). Empty means usable. Reported
        rather than raised so a GUI can mark individual fields; callers wanting
        an exception can raise on a non-empty result.
        """
        errors = []

        if self.period_est <= 0:
            errors.append(('period_est', "must be greater than 0 nm"))
        if self.scan_x_size <= 0:
            errors.append(('scan_x_size', "must be greater than 0 um"))

        if not 0 <= self.facet_trim <= MAX_FACET_TRIM:
            errors.append(('facet_trim',
                           f"must be between 0 and {MAX_FACET_TRIM}; above that "
                           f"the trim consumes the whole facet and no angle can "
                           f"be fitted"))

        if self.blaze_side not in VALID_BLAZE_SIDES:
            errors.append(('blaze_side',
                           f"must be one of {', '.join(VALID_BLAZE_SIDES)}"))
        if self.image_flatten_method not in VALID_IMAGE_FLATTEN_METHODS:
            errors.append(('image_flatten_method',
                           f"must be one of "
                           f"{', '.join(VALID_IMAGE_FLATTEN_METHODS)}"))

        if self.tip_correction not in VALID_TIP_CORRECTIONS:
            errors.append(('tip_correction',
                           f"must be one of {', '.join(VALID_TIP_CORRECTIONS)}"))
        if self.tip_correction != 'none':
            if self.tip_radius_nm <= 0:
                errors.append(('tip_radius_nm', "must be greater than 0 nm"))
            if not 0 < self.tip_half_angle_deg < 90:
                errors.append(('tip_half_angle_deg',
                               "must lie strictly between 0 and 90 degrees; "
                               "90 would be a flat punch and 0 an infinitely "
                               "thin needle"))

        if self.flatten_method not in VALID_FLATTEN_METHODS:
            errors.append(('flatten_method',
                           f"must be one of {', '.join(VALID_FLATTEN_METHODS)}"))
        if self.spm_direction not in VALID_SPM_DIRECTIONS:
            errors.append(('spm_direction',
                           f"must be one of {', '.join(VALID_SPM_DIRECTIONS)}"))

        if self.flatten_feature not in VALID_FLATTEN_FEATURES:
            errors.append(('flatten_feature',
                           f"must be one of {', '.join(VALID_FLATTEN_FEATURES)}"))

        if self.use_row_groups and self.n_row_groups < 2:
            errors.append(('n_row_groups', "needs at least 2 groups"))
        if self.ggp_n_points < 10:
            errors.append(('ggp_n_points', "needs at least 10 points"))
        if self.ggp_smoothing_window < 1:
            errors.append(('ggp_smoothing_window', "must be at least 1"))
        if self.ggp_min_half_width < 1:
            errors.append(('ggp_min_half_width', "must be at least 1"))

        if self.edge_exclusion_periods < 0:
            errors.append(('edge_exclusion_periods', "cannot be negative"))
        if not 0 <= self.flatten_exclude_edges < 0.5:
            errors.append(('flatten_exclude_edges', "must be between 0 and 0.5"))
        if self.flatten_poly_order < 1:
            errors.append(('flatten_poly_order', "must be at least 1"))

        return tuple(errors)
