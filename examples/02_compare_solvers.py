"""The same grating through two solvers, and how far apart they are.

The scalar solver is fast and approximate. The RCWA is rigorous and slower.
``compare.sweep`` runs several solvers on one wavelength grid and
``compare.align`` lines their orders up, so the difference is one call.

Run it:  python examples/02_compare_solvers.py     (takes a few seconds)
"""

import numpy as np

from gratinglab import Blazed, Illumination, Problem, align, sweep

problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0), coating="typical:Au")
light = Illumination.offplane(graze=1.5, azimuth=30.0)
wavelengths = np.linspace(1.0, 5.0, 9)

# Solver settings go in by name. 20 is the RCWA's default truncation; it is
# written out here to show where an accuracy setting goes.
scans = sweep(
    problem, light, wavelengths,
    methods=["scalar", "rcwa"],
    options={"rcwa": {"truncation": 20}},
)
comparison = align(scans)

print("order 1 efficiency")
print("  wavelength   scalar     rcwa   difference")
both = comparison.order(1)
for wl, s, r in zip(wavelengths, both["scalar"], both["rcwa"]):
    print(f"  {wl:7.2f} nm  {s:7.4f}  {r:7.4f}  {s - r:+9.4f}")

summary = comparison.summary("scalar", "rcwa")
print(f"\nlargest difference over all orders: {summary['max_abs_difference']:.3f}, "
      f"in order {summary['at_order']} at {summary['at_wavelength']:.2f} nm")
print(f"RMS difference: {summary['rms_difference']:.3f}")

# The rigorous solver also says where the light that was not diffracted went.
rcwa = scans[1]
print(f"\nrcwa at {wavelengths[0]:.1f} nm: reflected {rcwa.total[0]:.3f}, "
      f"absorbed {rcwa.absorption[0]:.3f}")

# A difference this size is the scalar approximation, not the RCWA: doubling
# its truncation moves these numbers by about 0.001. To check that for your
# own grating, see gratinglab.convergence.check_convergence.
