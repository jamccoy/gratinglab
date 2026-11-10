"""Groove profiles, in the four representations the solver methods need.

Profiles are normalised. Both the coordinate ``t`` and the returned height
are expressed in units of the period, so a profile is a pure shape and the
period is stored once, on :class:`~gratinglab.problem.Problem`. A solver that needs
physical height multiplies by ``problem.period``. This also matches PCGrate's
``.ggp`` format, which stores ``Period: 1`` and normalised points.

The four methods require different geometry from the same profile:

=================  ==========================================================
Method             Representation
=================  ==========================================================
Scalar             :meth:`Profile.height` -- the groove phase function
RCWA               :meth:`Profile.slice_layers` -- piecewise-constant slabs
C-method           :meth:`Profile.slope` -- smooth dy/dt; no vertical facets
Integral           :meth:`Profile.boundary` -- parametrised curve with
                   outward normals and arc-length quadrature weights
=================  ==========================================================

All four were designed together deliberately, because retrofitting normals and
arc lengths onto a shipped API once the integral method arrived would have
required changes to every profile type.

A profile that cannot supply a representation raises an error. ``Blazed`` with
a vertical anti-blaze facet has no finite slope there, so :meth:`Blazed.slope`
refuses instead of returning a large number that would silently corrupt a
C-method solve.
"""

from __future__ import annotations

# Importers use the names on this module. The parts are how the source is
# arranged, not an interface. The names missing from __all__ are the ones the
# solvers reach for, which is why they are imported without being used here.
from ._polyline import (  # noqa: F401
    CORNER_TURN_DEGREES,
    _polyline_boundary,
    _runs_to_intervals,
    polyline_corner_count,
    polyline_has_bends,
    polyline_has_corners,
    polyline_rounded_turn,
    smoothed_polygon_sides,
    smoothed_polyline,
)
from ._shapes import Blazed, FromProfileData, Lamellar, Sinusoidal
from ._trapezoid import Trapezoidal
from ._types import BoundaryCurve, Layer, Profile, ProfileRepresentationError

__all__ = [
    "Profile",
    "Layer",
    "BoundaryCurve",
    "ProfileRepresentationError",
    "Blazed",
    "Lamellar",
    "Trapezoidal",
    "Sinusoidal",
    "FromProfileData",
]
