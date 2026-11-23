"""gratinglab: an open comparison platform for grating efficiency analysis.

The organising idea is that the physical problem and the numerical method are
separate objects, and the problem spec is serializable. One ``Problem`` goes to
every solver, so disagreement between methods becomes visible instead of hidden
behind incompatible conventions.

``docs/conventions.md`` is normative. In short: time dependence ``exp(-iwt)``,
lossy ``n = n' + ik``, grating equation
``sin(alpha) + sin(beta_m) = m*lambda/(p*sin(gamma))``, efficiencies absolute,
lengths in nm, angles in degrees at the API boundary.

>>> from gratinglab import Illumination
>>> Illumination.offplane(graze=1.5, azimuth=0.97).gamma_deg
1.5

What a script needs is importable from here::

    from gratinglab import Problem, Blazed, Illumination, solve, materials, Mirror

The grating, the profiles, the light and the result types are imported with
the package. The solvers, the comparison and the mirror are looked up on first
use (:func:`__getattr__`), so ``import gratinglab`` does not pay for scipy's
linear algebra until something is solved. ``examples/`` shows them in use.
"""

from __future__ import annotations

from importlib import import_module
from typing import TYPE_CHECKING

from . import materials
from .illumination import Illumination, Polarization
from .measured import MeasuredEfficiency, MeasurementProvenance
from .problem import ConformalLayer, Problem
from .profiles import Blazed, Lamellar, Sinusoidal, Trapezoidal
from .resolution import ResolvingPower, resolving_power
from .result import EfficiencyScan, OrderEfficiency, Provenance

if TYPE_CHECKING:  # for editors and type checkers; at run time these are lazy
    from .compare import align, sweep
    from .mirror import Mirror
    from .solvers import available_solvers, get_solver, solve

#: Names resolved on first use, and the module each lives in.
_LAZY = {
    "solve": "gratinglab.solvers",
    "get_solver": "gratinglab.solvers",
    "available_solvers": "gratinglab.solvers",
    "sweep": "gratinglab.compare",
    "align": "gratinglab.compare",
    "Mirror": "gratinglab.mirror",
}


def __getattr__(name: str):
    """Import a lazy name's module the first time the name is asked for."""
    if name in _LAZY:
        value = getattr(import_module(_LAZY[name]), name)
        globals()[name] = value  # so the next lookup is an ordinary one
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted({*globals(), *_LAZY})

__version__ = "0.5.2"

#: The day this version was released, ISO 8601. It sits beside the version
#: rather than being read from CHANGELOG.md at run time because the changelog
#: is not in the wheel -- only `docs/` is -- so an installed copy would have
#: nothing to read and the About box would be blank for everyone but a
#: developer. `tests/test_changelog.py` asserts that it matches the date on
#: the `## <version>` heading, so the release commit that forgets to move it
#: fails at the point where it is cheapest to fix, the same arrangement as
#: the version's own section.
__version_date__ = "2026-11-23"

__all__ = [
    "Problem",
    "ConformalLayer",
    "Blazed",
    "Lamellar",
    "Sinusoidal",
    "Trapezoidal",
    "Illumination",
    "Polarization",
    "solve",
    "get_solver",
    "available_solvers",
    "sweep",
    "align",
    "Mirror",
    "materials",
    "EfficiencyScan",
    "OrderEfficiency",
    "Provenance",
    "MeasuredEfficiency",
    "MeasurementProvenance",
    "ResolvingPower",
    "resolving_power",
    "__version__",
    "__version_date__",
]
