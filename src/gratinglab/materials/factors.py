r"""A material's X-ray index from tabulated atomic scattering factors.

```python
from gratinglab import materials

Problem(..., coating="factors:SiO2@2.2")   # needs Si.nff and O.nff
Problem(..., coating="factors:Au")         # a default density, see below
pmma = materials.factors("C5H8O2", density=1.19)
```

Each atom scatters forward with :math:`f_1 + i f_2` electrons, and in the
independent-atom approximation a material is the sum of its atoms. With
:math:`n` formula units per volume, the complex index
:math:`1 - \delta + i\beta` has

.. math:: \delta + i\beta = \frac{r_e \lambda^2}{2\pi}\, n
          \sum_{\text{atoms}} (f_1 + i f_2),

which is :func:`index` below, and is the relation CXRO's own calculator
uses. :mod:`.typical` uses the same function with modelled :math:`f_1, f_2`;
this module takes them from per-element tables you supply instead, so it
keeps the edge structure the model lacks. Both ignore chemical bonding, so
near an edge a measured table of the actual film is still better.

The absorption is photoabsorption only, because that is all :math:`f_2` holds.
CXRO's own index calculator gives a larger :math:`\beta` at high energy for
light elements, consistent with it also counting Compton scattering as loss.
Checked against it on 2026-11-17 with CXRO's ``au.nff``, ``si.nff`` and
``o.nff``: the decrement agrees to 1e-4 everywhere from 50 eV to 20 keV, and
so does Au's absorption (median 9e-5, worst 0.7% at its M edges), while
Si's absorption reads 3.4% low at 20 keV and O's 21% low, both under 0.5%
below 7 keV. For a mirror the decrement dominates there: SiO2's grazing
reflectivity (0.02 to 0.6 deg) differs from CXRO's by at most 0.23% at
20 keV. A transmission or attenuation calculation for a light material above
about 10 keV should not use this.

The tables are read from your material library as ``<Symbol>.nff``, the
three-column ``E(eV) f1 f2`` format CXRO publishes one element at a time
(``https://henke.lbl.gov/optical_constants/asf.html``). A row whose
:math:`f_1` is -9999 is one where the compilation gives none, and is skipped.
None is bundled, for the reason ``data/SOURCES.md`` gives for every table.

Read by :func:`gratinglab.materials.lookup` for a ``factors:`` name, which
finds the files; nothing here touches a directory.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .elements import weight
from .optical import HC_EV_NM, OpticalConstants, _numeric_row

__all__ = ["PREFIX", "SUFFIX", "ScatteringFactors", "read_nff", "parse_formula", "index", "factors"]

#: The prefix that sends a material name here instead of to a table library.
PREFIX = "factors:"
#: The file suffix a library holds one element's factors under.
SUFFIX = ".nff"

#: Classical electron radius, m (CODATA 2018).
R_E = 2.8179403262e-15
#: Avogadro constant, 1/mol (exact, SI 2019).
N_A = 6.02214076e23

#: Henke's placeholder for an energy with no tabulated f1.
_NO_F1 = -9000.0

_TERM = re.compile(r"([A-Z][a-z]?)(\d+(?:\.\d+)?|\.\d+)?")


@dataclass(frozen=True)
class ScatteringFactors:
    """One element's :math:`f_1, f_2` against ascending energy in eV."""

    symbol: str
    energy_ev: NDArray[np.float64]
    f1: NDArray[np.float64]
    f2: NDArray[np.float64]
    source: str


def read_nff(path: "str | Path", symbol: str | None = None) -> ScatteringFactors:
    """Read an ``E(eV) f1 f2`` table; ``symbol`` defaults to the file stem."""
    path = Path(path)
    symbol = symbol or path.stem.capitalize()
    weight(symbol)  # refuses a stem that is not an element
    rows = [
        row
        for line in path.read_text(encoding="utf-8").splitlines()
        if (row := _numeric_row(line)) is not None and row[1] > _NO_F1
    ]
    if len(rows) < 2:
        raise ValueError(
            f"{path}: expected rows of three numbers (E in eV, f1, f2) with f1 "
            f"given, found {len(rows)}"
        )
    data = np.asarray(rows, dtype=np.float64)
    energy = data[:, 0]
    if (energy <= 0).any() or (np.diff(energy) < 0).any():
        raise ValueError(f"{path}: energies must be positive and ascending")
    return ScatteringFactors(symbol, energy, data[:, 1], data[:, 2], f"{path.name}")


def parse_formula(formula: str) -> tuple[tuple[str, float], ...]:
    """``"SiO2"`` as ``(("Si", 1.0), ("O", 2.0))``; counts may be fractional."""
    if not formula or "".join(m.group(0) for m in _TERM.finditer(formula)) != formula:
        raise ValueError(f"{formula!r} is not a chemical formula such as Au, SiO2 or C5H8O2")
    counts: dict[str, float] = {}
    for symbol, count in _TERM.findall(formula):
        weight(symbol)
        n = float(count) if count else 1.0
        if n <= 0:
            raise ValueError(f"{formula!r}: {symbol} needs a positive count")
        counts[symbol] = counts.get(symbol, 0.0) + n
    return tuple(counts.items())


def molar_mass(formula: str) -> float:
    """Grams per mole of formula units."""
    return sum(weight(s) * n for s, n in parse_formula(formula))


def atoms_per_m3(formula: str, density: float) -> float:
    """Atoms per cubic metre at ``density`` g/cm^3 (the :math:`N_a` of a table
    of number densities, counting every atom rather than formula units)."""
    atoms = sum(n for _, n in parse_formula(formula))
    return atoms * density * 1e6 / molar_mass(formula) * N_A


def index(
    f1: ArrayLike, f2: ArrayLike, mass: float, density: float, wavelength_nm: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """``(decrement, absorption)`` from :math:`f_1, f_2` summed over one formula
    unit of molar mass ``mass`` g/mol, at ``density`` g/cm^3."""
    formula_units = density * 1e6 / mass * N_A  # per m^3
    scale = R_E * (np.asarray(wavelength_nm) * 1e-9) ** 2 * formula_units / (2 * np.pi)
    return scale * np.asarray(f1), scale * np.asarray(f2)


def factors(
    formula: str,
    density: float,
    tables: Mapping[str, ScatteringFactors],
    name: str | None = None,
) -> OpticalConstants:
    """The index of ``formula`` at ``density`` g/cm^3 from per-element ``tables``.

    The result is sampled on every energy any of its elements is tabulated at,
    so each element's edges keep their own sampling, over the range all of
    them cover. Outside it a lookup is refused, as a table's is.
    """
    if not density > 0:
        raise ValueError(f"density must be positive, got {density!r}")
    composition = parse_formula(formula)
    missing = [s for s, _ in composition if s not in tables]
    if missing:
        raise KeyError(f"no scattering factors for {', '.join(missing)}")
    used = [tables[s] for s, _ in composition]
    low = max(t.energy_ev[0] for t in used)
    high = min(t.energy_ev[-1] for t in used)
    if not low < high:
        raise ValueError(f"{formula}: its elements' tables share no energy range")
    grid = np.unique(np.concatenate([t.energy_ev for t in used]))
    energy = grid[(grid >= low) & (grid <= high)]
    f1 = np.zeros_like(energy)
    f2 = np.zeros_like(energy)
    for (symbol, n), table in zip(composition, used):
        f1 += n * np.interp(energy, table.energy_ev, table.f1)
        f2 += n * np.interp(energy, table.energy_ev, table.f2)
    wavelength = HC_EV_NM / energy[::-1]
    decrement, absorption = index(f1[::-1], f2[::-1], molar_mass(formula), density, wavelength)
    files = ", ".join(t.source for t in used)
    return OpticalConstants(
        name=name or f"{PREFIX}{formula}@{density:g}",
        wavelength_nm=wavelength,
        decrement=decrement,
        absorption=absorption,
        source=(
            f"tabulated scattering factors ({files}), {formula} at {density:g} "
            "g/cm^3, independent atoms; see gratinglab.materials.factors"
        ),
    )
