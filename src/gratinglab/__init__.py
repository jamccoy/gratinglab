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
"""

from __future__ import annotations

from .illumination import Illumination, Polarization
from .measured import MeasuredEfficiency, MeasurementProvenance
from .resolution import ResolvingPower, resolving_power
from .result import EfficiencyScan, OrderEfficiency, Provenance

__version__ = "0.3.6"

#: The day this version was released, ISO 8601. It sits beside the version
#: rather than being read from CHANGELOG.md at run time because the changelog
#: is not in the wheel -- only `docs/` is -- so an installed copy would have
#: nothing to read and the About box would be blank for everyone but a
#: developer. `tests/test_changelog.py` asserts that it matches the date on
#: the `## <version>` heading, so the release commit that forgets to move it
#: fails at the point where it is cheapest to fix, the same arrangement as
#: the version's own section.
__version_date__ = "2026-11-13"

__all__ = [
    "Illumination",
    "Polarization",
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
