"""
Boundary profile generation for PCGrate

Averages the detected grooves of an AFM scan into a single representative
groove, normalises it to one period, and exports it as a PCGrate .ggp
boundary profile.

Shares the profile front-end (loading, groove detection, groove windowing)
with the blaze-angle analysis in afm_analysis.core.
"""

from .average import (
    flatten_endpoints,
    average_grooves,
    normalize_profile,
    profile_metrics,
)

__all__ = [
    'flatten_endpoints',
    'average_grooves',
    'normalize_profile',
    'profile_metrics',
]
