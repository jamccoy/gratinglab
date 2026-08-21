"""
Groove metrology: measured surfaces, not modelled ones.

Reads AFM scans of a diffraction grating and produces two things: blaze angles
per groove with row-group statistics, and one averaged groove normalised to a
period -- the boundary profile that
:class:`gratinglab.profiles.FromProfileData` consumes and every solver can be
driven from.

This is the half of the project that measures. The rest of :mod:`gratinglab`
computes what a given geometry would do; here the geometry arrives from an
instrument, with an uncertainty attached. Conventions are shared and normative:
see ``docs/conventions.md``, particularly §1 on units -- metres at the
instrument boundary, microns laterally, nanometres in profile space, and a
dimensionless fraction of the period at the boundary-file boundary.

Requires the ``metrology`` extra::

    pip install -e '.[metrology]'
"""
from __future__ import annotations

__all__ = [
    "require_matplotlib",
    "MATPLOTLIB_MISSING_MESSAGE",
    "run_single_file_analysis",
    "run_multiple_file_analysis",
    "run_comparison_analysis",
    "analyze_single_file",
]

MATPLOTLIB_MISSING_MESSAGE = (
    "gratinglab.metrology needs matplotlib, which is an optional extra.\n"
    "Install it with:\n\n"
    "    pip install -e '.[metrology]'\n\n"
    "The solvers, the comparison harness and the .ggp reader all work without it."
)


def require_matplotlib() -> None:
    """Raise something useful rather than a bare ImportError.

    Mirrors :func:`gratinglab.metrology.gui.require_qt`, and for the same
    reason: ``No module named 'matplotlib'`` tells a user nothing about which
    extra they are missing, and this is the most likely first-run failure for
    anyone who installed the core distribution.

    ``find_spec`` rather than an import, so checking availability does not cost
    a full matplotlib import in the case where it *is* present.
    """
    from importlib.util import find_spec

    if find_spec("matplotlib") is None:  # pragma: no cover - depends on environment
        raise ModuleNotFoundError(MATPLOTLIB_MISSING_MESSAGE)


# Before the re-exports below, which reach plotting code transitively. Checking
# here is what turns an unhelpful traceback from four modules down into one line
# naming the extra.
require_matplotlib()

from .analyzer import analyze_single_file  # noqa: E402
from .workflows import (  # noqa: E402
    run_comparison_analysis,
    run_multiple_file_analysis,
    run_single_file_analysis,
)
