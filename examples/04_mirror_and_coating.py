"""A coating as a flat mirror: reflectivity, penetration depth, and when a
layer is thick enough to treat as bulk.

No grating here. ``Mirror`` answers what a coating does before grooves are
cut into it, which decides how the coating should be modelled on a grating.

Run it:  python examples/04_mirror_and_coating.py
"""

import numpy as np

from gratinglab.mirror import Mirror
from gratinglab.problem import ConformalLayer

# 20 nm of gold on 5 nm of nickel on silicon. Layers are listed innermost
# first, so the last one is the surface the light meets.
mirror = Mirror(
    substrate="typical:Si",
    layers=(
        ConformalLayer(material="typical:Ni", thickness=5.0),
        ConformalLayer(material="typical:Au", thickness=20.0),
    ),
    roughness=0.3,  # nm RMS, at every interface
)

graze = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0])  # degrees from the surface
wavelength = 1.0  # nm

# Reflectivity broadcasts over angle and wavelength.
reflectivity = mirror.reflectivity(graze, wavelength)
bare_gold = Mirror("typical:Au", roughness=0.3).reflectivity(graze, wavelength)

# The 1/e intensity depth in each material, taken as bulk.
depth = mirror.penetration_depth(graze, wavelength)

print(f"at {wavelength:g} nm")
print("  graze    R (stack)   R (bulk Au)   1/e depth in Au")
for g, r, b, d in zip(graze, reflectivity, bare_gold, depth["typical:Au"]):
    print(f"  {g:4.1f} deg {r:9.4f}   {b:11.4f}   {d:8.2f} nm")

# How thick must the gold be before the nickel and silicon under it change
# the reflectivity by less than one part in a thousand?
print("\ngold needed to read as bulk (within 0.1% of R):")
for g in (1.0, 1.5, 3.0):
    needed = mirror.bulk_thickness(g, wavelength, tolerance=1e-3)
    verdict = "20 nm is enough" if needed <= 20.0 else "20 nm is not enough"
    print(f"  at {g:.1f} deg: {needed:7.1f} nm   ({verdict})")

# So at 1.5 degrees this gold may be modelled as the grating material alone,
# Problem(coating="typical:Au"), and at 3 degrees the layers must be modelled:
# Problem(coating="typical:Si", layers=mirror.layers), which the "rcwa" solver takes.

# Amplitudes, with phase, for anything that needs more than the intensity.
r_s, r_p = mirror.amplitude(1.5, wavelength)
print(f"\nat 1.5 deg: r_s = {complex(r_s):.4f}, r_p = {complex(r_p):.4f}")
