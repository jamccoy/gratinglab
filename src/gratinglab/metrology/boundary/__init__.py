"""
Boundary profile generation for PCGrate

This package averages the detected grooves of an AFM scan into a single
representative groove, normalises it to one period, and exports it as a
PCGrate .ggp boundary profile.

It shares the profile front-end (loading, groove detection, groove windowing)
with the blaze-angle analysis in gratinglab.metrology.core.
"""

from .pipeline import BoundaryProfile, build_boundary_profile
from .swaths import SwathStatistics, swath_statistics
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
    'flatten_endpoints',
    'average_grooves',
    'groove_residuals',
    'normalize_profile',
    'profile_metrics',
    'select_grooves',
]
