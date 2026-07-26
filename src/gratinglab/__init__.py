"""gratinglab -- an open comparison platform for grating efficiency analysis.

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
from .result import EfficiencyScan, OrderEfficiency, Provenance

__version__ = "0.0.1"

__all__ = [
    "Illumination",
    "Polarization",
    "EfficiencyScan",
    "OrderEfficiency",
    "Provenance",
    "__version__",
]
