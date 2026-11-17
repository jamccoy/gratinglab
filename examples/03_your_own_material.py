"""Using optical constants of your own, and the built-in ways to name a material.

A material reaches a calculation by name. This shows the four kinds of name:

- a table in your library (``"Glass"``),
- the built-in X-ray model (``"typical:Au"``),
- tabulated scattering factors with a formula and density (``"factors:SiO2@2.2"``),
- a library table at another density (``"Au@17.5"``).

Run it:  python examples/03_your_own_material.py

It writes one small made-up table to a temporary folder, so it changes
nothing in your real library.
"""

import os
import tempfile
from pathlib import Path

from gratinglab import materials

# -- 1. A table of your own -------------------------------------------------------
# Columns of wavelength, n and k, such as an ellipsometry export. The header
# must name the wavelength unit (nm, um, A or eV); a file that does not is
# refused rather than guessed at. These numbers are invented.
folder = Path(tempfile.mkdtemp(prefix="gratinglab-example-"))
table = folder / "Glass.nk"
table.write_text(
    "# an invented glass, for the example\n"
    "wavelength (nm), n, k\n"
    "400, 1.530, 0.0\n"
    "550, 1.518, 0.0\n"
    "700, 1.513, 0.0\n",
    encoding="utf-8",
)

# Read one file directly, wherever it is...
glass = materials.from_nk_file(table)
low, high = glass.range_nm
print(f"{glass.name}: covers {low:g} to {high:g} nm")
print(f"  index at 500 nm: {complex(glass.n(500.0)):.4f}")

# ...or keep it in a library, a folder of such files, and ask for it by name.
# Your library is ~/.gratinglab/materials; GRATINGLAB_MATERIALS names another.
os.environ["GRATINGLAB_MATERIALS"] = str(folder)
print(f"  in the library: {materials.available()}")
print(f"  looked up by name: {materials.lookup('Glass').source}")

# Asking outside a table's range is an error, never an extrapolation.
try:
    glass.n(900.0)
except ValueError as refused:
    print(f"  at 900 nm: refused ({str(refused)[:60]}...)")

# -- 2. The built-in X-ray model ---------------------------------------------------
# No table needed: Au, Pt, Ni, C, Si and O, and compounds of them, 30 eV to
# 30 keV. It is approximate, and says so in its source and its caveats.
gold = materials.lookup("typical:Au")
decrement, absorption = gold.at(1.0)
print(f"\ntypical:Au at 1 nm: 1 - n = {float(decrement):.3e}, k = {float(absorption):.3e}")
film = materials.lookup("typical:C@2.0")  # carbon at 2.0 g/cm^3
print(f"typical:C@2.0 is: {film.source[:52]}...")
for caveat in gold.caveats([0.55, 0.57]):  # across gold's M edges
    print(f"  caveat near 0.56 nm: {caveat[:70]}...")

# -- 3. Scattering factors, if you have them ---------------------------------------
# With CXRO's per-element files (Si.nff, O.nff) in your library,
# "factors:SiO2@2.2" computes any compound at any density from them.
elements = materials.factor_elements()
if elements:
    print(f"\nscattering factors installed for: {', '.join(elements)}")
else:
    print("\nno scattering-factor files installed, so factors: names are refused:")
    try:
        materials.lookup("factors:SiO2@2.2")
    except materials.UnknownMaterial as missing:
        print(f"  {str(missing)[:96]}...")

# -- 4. A table at another density ---------------------------------------------------
# "Au@17.5" rescales a CXRO export of gold to 17.5 g/cm^3, for a sputtered
# film. It needs a CXRO table named Au in the library, so it is only named here.
