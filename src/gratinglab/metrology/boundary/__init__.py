"""
Boundary profile generation for PCGrate

This package averages the detected grooves of an AFM scan into a single
representative groove, normalises it to one period, and exports it as a
PCGrate .ggp boundary profile. `trace_cross_section` does the same for an SEM
micrograph of a cleaved cross-section, and fits a trapezoid to it.

It shares the profile front-end (loading, groove detection, groove windowing)
with the blaze-angle analysis in gratinglab.metrology.core.
"""

from .pipeline import BoundaryProfile, build_boundary_profile
from .swaths import SwathStatistics, swath_statistics
from .cross_section import CrossSection, TrapezoidFit, trace_cross_section
from .average import (
    flatten_endpoints,
    average_grooves,
    groove_residuals,
    normalize_profile,
    profile_metrics,
    select_grooves,
)

__all__ = [
    'BoundaryProfile',
    'build_boundary_profile',
    'SwathStatistics',
    'swath_statistics',
    'CrossSection',
    'TrapezoidFit',
    'trace_cross_section',
    'flatten_endpoints',
    'average_grooves',
    'groove_residuals',
    'normalize_profile',
    'profile_metrics',
    'select_grooves',
]
