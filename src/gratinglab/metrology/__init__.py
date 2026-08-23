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

from importlib import import_module
from typing import TYPE_CHECKING

from endstation.optional import missing_extra_message, require_extra

__all__ = [
    "require_matplotlib",
    "MATPLOTLIB_MISSING_MESSAGE",
    "run_single_file_analysis",
    "run_multiple_file_analysis",
    "run_comparison_analysis",
    "analyze_single_file",
]

#: Kept as a module constant because the metrology tests assert on it -- both
#: that it names matplotlib and that it names the extra.
MATPLOTLIB_MISSING_MESSAGE = missing_extra_message(
    "matplotlib",
    "metrology",
    still_works=(
        "The solvers, the comparison harness and the .ggp reader all work "
        "without it."
    ),
)


def require_matplotlib() -> None:
    """Raise something useful rather than a bare ImportError.

    Mirrors :func:`endstation.qt.require_qt`, and for the same reason:
    ``No module named 'matplotlib'`` tells a user nothing about which extra they
    are missing, and this is the most likely first-run failure for anyone who
    installed the core distribution.

    The check itself is :func:`endstation.optional.require_extra`, which uses
    ``find_spec`` rather than an import so asking the question costs nothing in
    the case where the answer is yes.
    """
    require_extra(
        "matplotlib",
        "metrology",
        still_works=(
            "The solvers, the comparison harness and the .ggp reader all work "
            "without it."
        ),
    )


# Eager, and cheap enough to be: `require_extra` is a `find_spec` call, not an
# import. What it buys is that *any* attempt to use this package -- including
# importing a subpackage -- fails on one line naming the extra, rather than on
# an unhelpful traceback from four modules down the first time an analysis
# function is touched.
require_matplotlib()

#: The public names, and the submodule each is fetched from on first access.
#: They used to be imported here at module scope, which meant `import
#: gratinglab.metrology.gui` paid for matplotlib -- roughly a second -- before
#: `require_qt()` could say a word about a missing PySide6. The console script
#: `gratinglab-metrology-gui` exists partly to make that message fast, so the
#: parent package now hands these out lazily (PEP 562).
#:
#: This does not weaken the matplotlib check above: that still runs on import.
#: Only the *plotting code* is deferred, to the first caller who wants it.
_LAZY_EXPORTS = {
    "analyze_single_file": ".analyzer",
    "run_comparison_analysis": ".workflows",
    "run_multiple_file_analysis": ".workflows",
    "run_single_file_analysis": ".workflows",
}

if TYPE_CHECKING:  # so type checkers and IDEs still resolve the re-exports
    from .analyzer import analyze_single_file
    from .workflows import (
        run_comparison_analysis,
        run_multiple_file_analysis,
        run_single_file_analysis,
    )


def __getattr__(name: str):
    """Import the submodule behind ``name`` the first time it is asked for."""
    try:
        module = _LAZY_EXPORTS[name]
    except KeyError:
        raise AttributeError(
            f"module {__name__!r} has no attribute {name!r}"
        ) from None

    value = getattr(import_module(module, __name__), name)
    # Bind it here so this runs once per name; `__getattr__` fires only for
    # attributes that are not found normally.
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    """Keep ``dir()`` and tab-completion showing names not yet imported."""
    return sorted({*globals(), *__all__})
