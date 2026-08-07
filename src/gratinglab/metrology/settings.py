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

    # Flattening
    flatten_method: str = 'level_grooves'
    flatten_poly_order: int = 2
    flatten_exclude_edges: float = 0.05
    flatten_feature: str = 'peaks'

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
            flatten_method=config.FLATTEN_METHOD,
            flatten_poly_order=config.FLATTEN_POLY_ORDER,
            flatten_exclude_edges=config.FLATTEN_EXCLUDE_EDGES,
            flatten_feature=config.FLATTEN_FEATURE,
            spm_channel=getattr(config, 'SPM_CHANNEL', 'Height Sensor'),
            spm_direction=getattr(config, 'SPM_DIRECTION', 'Retrace'),
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
        if self.edge_exclusion_periods < 0:
            errors.append(('edge_exclusion_periods', "cannot be negative"))
        if not 0 <= self.flatten_exclude_edges < 0.5:
            errors.append(('flatten_exclude_edges', "must be between 0 and 0.5"))
        if self.flatten_poly_order < 1:
            errors.append(('flatten_poly_order', "must be at least 1"))

        return tuple(errors)
