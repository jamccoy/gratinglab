"""Diffraction efficiency of a blazed grating, from scratch.

The three objects every calculation starts from:

- a ``Problem``: the grating itself (period, groove profile, material),
- an ``Illumination``: how the light arrives,
- a solver, asked for by name.

Run it:  python examples/01_efficiency_scan.py

It needs nothing installed beyond gratinglab. The coating is ``typical:Au``,
the built-in approximate X-ray model, so no optical-constant table is needed;
with a table of your own in your library, write ``coating="Au"`` instead.
"""

import numpy as np

from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed
from gratinglab.solvers import get_solver

# A 160 nm period sawtooth with a 30 degree blaze, in gold. Lengths are nm.
problem = Problem(
    period=160.0,
    profile=Blazed(blaze_angle=30.0),
    coating="typical:Au",
)

# Extreme off-plane mount, as used for soft X-rays: light grazes the grating
# at 1.5 degrees, turned 30 degrees about the normal (near Littrow for this
# blaze). Illumination.classical(alpha) is the ordinary in-plane mount.
light = Illumination.offplane(graze=1.5, azimuth=30.0)

wavelengths = np.linspace(1.0, 5.0, 41)  # nm

# "scalar" is the fast approximate solver. "rcwa", "integral" and "cmethod"
# are the rigorous ones; see 02_compare_solvers.py.
scan = get_solver("scalar").solve(problem, light, wavelengths)

# A scan is a table: one row per wavelength, one column per order.
print(f"solved {len(scan)} wavelengths for orders {scan.orders.min()} to {scan.orders.max()}")
print(f"coating: {scan.provenance.notes['coating']}\n")

# One order across the scan...
first = scan.order(1)
best = int(np.argmax(first))
print(f"order 1 peaks at {wavelengths[best]:.2f} nm with efficiency {first[best]:.3f}")

# ...or every order at one wavelength.
row = scan.at(2.0)
print("\nat 2.0 nm:")
for order in row.propagating_orders:
    print(f"  order {order:+d}: {row[int(order)]:.4f}")
print(f"  total reflected into orders: {row.total:.4f}")

# to_records() gives plain dicts, ready for pandas.DataFrame or csv.DictWriter.
print(f"\nfirst record: {scan.to_records()[0]}")
