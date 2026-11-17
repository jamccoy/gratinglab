"""Something of your own: sweep a design parameter the GUI has no slider for.

Which blaze angle puts the most light into first order at 2 nm? A script can
ask that directly, because a ``Problem`` is just data: build one per value,
solve, and keep the number you care about.

Run it:  python examples/05_parameter_sweep.py
"""

import numpy as np

from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed
from gratinglab.solvers import get_solver

solver = get_solver("scalar")
wavelength = np.array([2.0])  # nm; a scan of one wavelength
graze = 1.5  # degrees


def first_order_efficiency(blaze_angle: float) -> float:
    """Order +1 at 2 nm for a 160 nm gold sawtooth with this blaze angle,
    in the Littrow-like mount where the azimuth equals the blaze."""
    problem = Problem(
        period=160.0,
        profile=Blazed(blaze_angle=blaze_angle),
        coating="typical:Au",
    )
    light = Illumination.offplane(graze=graze, azimuth=blaze_angle)
    return float(solver.solve(problem, light, wavelength).order(1)[0])


angles = np.arange(10.0, 46.0, 2.5)
efficiency = np.array([first_order_efficiency(a) for a in angles])

print("blaze angle   order +1 at 2 nm")
for angle, value in zip(angles, efficiency):
    bar = "#" * int(round(40 * value))
    print(f"   {angle:5.1f} deg    {value:.4f}  {bar}")

best = angles[int(np.argmax(efficiency))]
print(f"\nbest of these: {best:.1f} deg, efficiency {efficiency.max():.4f}")

# A Problem is immutable; model_copy makes a changed one. The same pattern
# sweeps the period, the coating, or a layer's thickness.
base = Problem(period=160.0, profile=Blazed(blaze_angle=best), coating="typical:Au")
for coating in ("typical:Au", "typical:Pt", "typical:Ni"):
    variant = base.model_copy(update={"coating": coating})
    light = Illumination.offplane(graze=graze, azimuth=best)
    value = float(solver.solve(variant, light, wavelength).order(1)[0])
    print(f"  {coating:11s} order +1: {value:.4f}")

# The scalar solver is the survey tool. Check a design it likes with a
# rigorous solver before trusting it: see 02_compare_solvers.py.
