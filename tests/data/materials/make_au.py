"""Write ``Au.txt``, the suite's stand-in for an X-ray gold table.

It is invented, not measured: two power laws with roughly gold's magnitudes,
written in CXRO's export format on the grid a CXRO request for 200 eV to 2 keV
in 500 points returns. The repository carries no third-party optical constants
(see ``src/gratinglab/materials/data/SOURCES.md``), and the solver tests need a
lossy X-ray table to exercise absorption, range guards and the readers, not
gold's true values. It has no absorption edges, so nothing here is evidence
about real gold.

``conftest.py`` makes this directory the bundled library for the suite, so
``coating="Au"`` resolves to this file in tests and nowhere else. Run it from
anywhere to rewrite the table; the output is deterministic.
"""

from pathlib import Path

import numpy as np

#: Decrement and absorption at 200 eV, and how fast each falls with energy.
#: Chosen by hand to be gold-sized, not fitted to any table.
DELTA_200, DELTA_SLOPE = 0.0113, 1.2
BETA_200, BETA_SLOPE = 0.0095, 1.7

energy = 200.0 * 10.0 ** (np.arange(500) / 499.0)
delta = DELTA_200 * (energy / 200.0) ** -DELTA_SLOPE
beta = BETA_200 * (energy / 200.0) ** -BETA_SLOPE

lines = [" Au synthetic, see tests/data/materials/make_au.py", " Energy(eV), Delta, Beta"]
lines += [f"  {e:.9g}  {d:.9g}  {b:.9g}" for e, d, b in zip(energy, delta, beta)]
out = Path(__file__).with_name("Au.txt")
out.write_bytes(("\n".join(lines) + "\n").encode("ascii"))
print(f"wrote {out}: {len(energy)} rows, {energy[0]:g} to {energy[-1]:g} eV")
