r"""A typical material's X-ray index, from physical facts instead of a table.

```python
from gratinglab import materials

au = materials.typical("Au")                  # bulk density
c = materials.typical("C", density=2.0)       # an evaporated film, say
Problem(..., coating="typical:SiO2")          # the same, by name
Problem(..., coating="typical:C@2.0")         # with a density, by name
```

The repository commits no tabulated optical constants (``data/SOURCES.md``),
so this builds an index from a handful of physical facts per element instead:
atomic number, atomic weight, the core-level binding energies and how many
electrons sit in each level, and a density per material. It is a model, and
it says so: every result carries ``source`` naming it, and a scan that comes
near an absorption edge carries a provenance warning (:meth:`caveats`).

The model
=========

Each atom scatters with :math:`f_1 + i f_2` electrons, and a material is the
sum of its atoms (the independent-atom approximation, which ignores chemical
bonding and so fails where bonding matters: near edges and at low energy).
Then, with :math:`n_a` atoms per volume,

.. math:: \text{decrement} + i\,\text{absorption}
          = \frac{r_e \lambda^2}{2\pi} \sum_{\text{atoms}} n_a (f_1 + i f_2).

The free-electron limit takes :math:`f_1 = Z`. That holds only above every
edge of the atom, because an electron bound more tightly than the photon's
energy cannot follow the field freely: far below its edge a shell adds almost
nothing to :math:`f_1`, and just below it adds a negative amount (anomalous
dispersion). For gold at 1 nm, where the K, L and M shells (28 electrons) are
all still bound, the free-electron value is about 1.5 times too large.

So each core shell :math:`s` (binding energy :math:`E_s`, occupancy
:math:`N_s`) absorbs by Kramers' law, :math:`\sigma \propto E^{-3}` above its
edge and nothing below it, which is :math:`f_2 \propto E^{-2}` since
:math:`f_2 \propto E\sigma`. Its strength is fixed by the shell's own sum rule,
:math:`(2/\pi)\int_0^\infty f_2(E)\,dE/E = N_s`, which says the shell's
:math:`N_s` electrons stop contributing to :math:`f_1` at zero frequency:

.. math:: f_{2,s}(E) = \pi N_s (E_s/E)^2 \quad (E > E_s).

The Kramers-Kronig relation then gives :math:`f_1` in closed form, with
:math:`x = E/E_s`:

.. math:: f_1(E) = Z^* + \sum_s N_s\,x^{-2} \ln\lvert x^2 - 1\rvert,

which tends to :math:`Z^*` far above every edge, to the valence count far
below, and dips logarithmically at each edge. A real edge is broadened by the
core hole's lifetime and has fine structure, so the dip is deepest exactly
where the model is least to be trusted; it is clamped at
:data:`_EDGE_CLAMP` so an energy exactly on an edge returns a number.
:math:`Z^* = Z - (Z/82.5)^{2.37}` is the relativistic correction of Henke,
Gullikson and Davis (1993).

The valence electrons (those with no listed core level) are given the same
treatment, with the atom's first ionization energy as their edge. Far above it
they are nearly free, adding about one each to :math:`f_1`, and the tail of
their Kramers absorption keeps :math:`f_2` above zero below the first core
edge. Without it carbon between 100 and 284 eV would have no absorption at
all, and a solver would treat it as a lossless dielectric that transmits. The
ionization energy is an atomic value, and in a solid the valence oscillator
strength sits higher, around the plasmon; that is why the model stops at
:data:`LOW_EV`: below about 30 eV the valence response is the whole story
(plasmons, interband transitions), and nothing this simple describes it.
It also stops at :data:`HIGH_EV`, past which Compton scattering, which is not
in :math:`f_2`, matters to absorption.

Where it is approximate
=======================

Kramers' law puts a shell's whole absorption at its edge. For d and f shells
in heavy atoms that is wrong in a known way: the centrifugal barrier holds the
outgoing electron back, so absorption rises over several times the binding
energy instead (the delayed onset), and gold's and platinum's 4d and 4f shells
(335-353 eV and 84-88 eV) are where the model is worst; up to
:data:`DELAYED_REACH` times the binding energy of a d or f level is reported
by :meth:`TypicalConstants.caveats` as well, for the 4d and 4f levels. The
light elements have no such shells.

The per-shell sum rule is approximate in a many-electron atom (oscillator
strength moves between shells), Kramers' law is a semiclassical estimate, and
real edges have structure the step lacks. Expect :math:`f_1`, and so the
decrement, to be good between edges and well above them, and :math:`f_2`, the
absorption, to be cruder. Within :data:`EDGE_FRACTION` of an edge (or
:data:`EDGE_MIN_EV`, whichever is wider) the index is qualitative and
:meth:`TypicalConstants.caveats` says so; so it is below :data:`VALENCE_EV`.
How good it is elsewhere is to be measured against real tables, not claimed
here.

The facts
=========

Binding energies: X-ray Data Booklet (Lawrence Berkeley National Laboratory),
Table 1-1, G. P. Williams, compiled from J. A. Bearden and A. F. Burr, Rev.
Mod. Phys. 39, 125 (1967), M. Cardona and L. Ley (1978), and J. C. Fuggle and
N. Martensson, J. Electron Spectrosc. 21, 275 (1980). Occupancies follow from
the ground-state configuration. First ionization energies: CRC Handbook of
Chemistry and Physics, "Ionization Energies of Atoms". Atomic weights: IUPAC standard atomic weights.
Densities: bulk values as stated in the CRC Handbook of Chemistry and Physics,
except carbon, see :data:`DENSITIES`. Each is a measured physical quantity
stated with its source, not a reproduced table.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .optical import HC_EV_NM, OpticalConstants

__all__ = [
    "ELEMENTS",
    "DENSITIES",
    "TypicalConstants",
    "typical",
    "from_name",
    "PREFIX",
]

#: Classical electron radius, m (CODATA 2018).
R_E = 2.8179403262e-15
#: Avogadro constant, 1/mol (exact, SI 2019).
N_A = 6.02214076e23

#: The model's range, eV. Requests outside it are refused, as a table's are.
LOW_EV, HIGH_EV = 30.0, 30000.0
#: Below this the valence response matters and the index is qualitative.
VALENCE_EV = 100.0
#: Half-width of the band around an edge that is reported as near it.
EDGE_FRACTION, EDGE_MIN_EV = 0.10, 30.0
#: The 4d and 4f levels, whose absorption is delayed to several times their
#: binding energy, and how far above the edge that is reported (as a multiple).
#: Deeper d levels (gold's M4,5 at 2.2 keV) bind too tightly for the barrier to
#: matter much, and measured against CXRO they show no such error; 4x is what
#: covers platinum's 4f, whose excess absorption still reads 180% at 3.3x.
DELAYED_LEVELS = frozenset({"N4", "N5", "N6", "N7"})
DELAYED_REACH = 4.0
#: Smallest |x^2 - 1| the edge logarithm is evaluated at.
_EDGE_CLAMP = 1e-6

#: The prefix that sends a material name here instead of to a table library.
PREFIX = "typical:"


@dataclass(frozen=True)
class Element:
    """One element's facts: ``shells`` is ``(level, binding eV, electrons)``
    for the core levels, and ``ionization`` (eV) is the valence electrons' edge."""

    symbol: str
    z: int
    weight: float
    ionization: float
    shells: tuple[tuple[str, float, int], ...]

    @property
    def valence(self) -> int:
        return self.z - sum(n for _, _, n in self.shells)

    @property
    def absorbers(self) -> tuple[tuple[str, float, int], ...]:
        """The core shells and then the valence electrons, as the model uses them."""
        return (*self.shells, ("valence", self.ionization, self.valence))


def _gold_like(k, l1, l2, l3, m, n, o):
    """The 68 core electrons of Pt and Au, which differ only in energies."""
    return (
        ("K", k, 2), ("L1", l1, 2), ("L2", l2, 2), ("L3", l3, 4),
        *zip(("M1", "M2", "M3", "M4", "M5"), m, (2, 2, 4, 4, 6)),
        *zip(("N1", "N2", "N3", "N4", "N5", "N6", "N7"), n, (2, 2, 4, 4, 6, 6, 8)),
        *zip(("O1", "O2", "O3"), o, (2, 2, 4)),
    )


#: Atomic weight, first ionization energy (eV) and core binding energies (eV,
#: X-ray Data Booklet Table 1-1) with each level's occupancy.
ELEMENTS = {
    e.symbol: e
    for e in (
        Element("C", 6, 12.011, 11.2603, (("K", 284.2, 2),)),
        Element("O", 8, 15.999, 13.6181, (("K", 543.1, 2), ("L1", 41.6, 2))),
        Element("Si", 14, 28.085, 8.1517, (
            ("K", 1839.0, 2), ("L1", 149.7, 2), ("L2", 99.82, 2), ("L3", 99.42, 4),
        )),
        Element("Ni", 28, 58.693, 7.6399, (
            ("K", 8333.0, 2), ("L1", 1008.6, 2), ("L2", 870.0, 2), ("L3", 852.7, 4),
            ("M1", 110.8, 2), ("M2", 68.0, 2), ("M3", 66.2, 4),
        )),
        Element("Pt", 78, 195.084, 8.9588, _gold_like(
            78395.0, 13880.0, 13273.0, 11564.0,
            (3296.0, 3027.0, 2645.0, 2202.0, 2122.0),
            (725.4, 609.1, 519.4, 331.6, 314.6, 74.5, 71.2),
            (101.7, 65.3, 51.7),
        )),
        Element("Au", 79, 196.96657, 9.2256, _gold_like(
            80725.0, 14353.0, 13734.0, 11919.0,
            (3425.0, 3148.0, 2743.0, 2291.0, 2206.0),
            (762.1, 642.7, 546.3, 353.2, 335.1, 87.6, 84.0),
            (107.2, 74.2, 57.2),
        )),
    )
}

#: Default densities, g/cm^3. Thin films often run a few percent below bulk,
#: which moves the decrement in proportion; pass ``density`` when it is known.
#: Carbon has no single bulk density (graphite 2.26, diamond 3.51); 2.2 is a
#: typical amorphous coating, and is a choice, not a fact.
DENSITIES = {
    "Au": 19.32, "Pt": 21.45, "Ni": 8.908, "Si": 2.329, "C": 2.2, "SiO2": 2.20,
}

_TERM = re.compile(r"([A-Z][a-z]?)(\d*)")


def _parse(formula: str) -> tuple[tuple[Element, int], ...]:
    """``"SiO2"`` as ``((Si, 1), (O, 2))``, refusing an element with no facts."""
    if not formula or "".join(m.group(0) for m in _TERM.finditer(formula)) != formula:
        raise ValueError(f"{formula!r} is not a chemical formula such as Au or SiO2")
    counts: dict[str, int] = {}
    for symbol, count in _TERM.findall(formula):
        if symbol not in ELEMENTS:
            raise ValueError(
                f"no facts for {symbol!r} in the typical-material model; it "
                f"knows {', '.join(ELEMENTS)}"
            )
        counts[symbol] = counts.get(symbol, 0) + int(count or 1)
    return tuple((ELEMENTS[s], n) for s, n in counts.items())


def scattering_factors(
    element: Element, energy_ev: ArrayLike
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """``(f1, f2)`` for one atom, in electrons, by the closed forms above."""
    energy = np.asarray(energy_ev, dtype=np.float64)
    z_star = element.z - (element.z / 82.5) ** 2.37
    f1 = np.full_like(energy, z_star)
    f2 = np.zeros_like(energy)
    for _, binding, electrons in element.absorbers:
        x2 = (energy / binding) ** 2
        # log1p below the edge, where x^2 can be too small for 1 - x^2 to
        # differ from 1 (gold's K shell seen from 30 eV).
        gap = np.maximum(np.abs(x2 - 1.0), _EDGE_CLAMP)
        log_gap = np.where(x2 < 1.0, np.log1p(-np.minimum(x2, 1.0 - _EDGE_CLAMP)), np.log(gap))
        f1 += electrons / x2 * log_gap
        f2 += np.where(x2 > 1.0, np.pi * electrons / x2, 0.0)
    return f1, f2


@dataclass(frozen=True, slots=True)
class TypicalConstants(OpticalConstants):
    """An :class:`OpticalConstants` evaluated from the model, not interpolated.

    The inherited arrays hold the model sampled across its range, for anything
    that plots or inspects a table; :meth:`at` evaluates the closed forms
    directly, so there is no interpolation error near an edge.
    """

    formula: str = ""
    density: float = 0.0
    composition: tuple = field(default=(), repr=False)

    def at(self, wavelength_nm):
        values = np.asarray(wavelength_nm, dtype=np.float64)
        self._require_in_range(values)
        return _index(self.composition, self.density, values)

    def edges(self) -> tuple[tuple[str, float], ...]:
        """Every core level in the material, as ``("Au M5", eV)``, by energy."""
        found = {
            (f"{el.symbol} {level}", binding)
            for el, _ in self.composition
            for level, binding, _ in el.shells
        }
        return tuple(sorted(found, key=lambda item: item[1]))

    def caveats(self, wavelength_nm) -> tuple[str, ...]:
        """Where these wavelengths reach the model's qualitative bands."""
        energy = HC_EV_NM / np.atleast_1d(np.asarray(wavelength_nm, dtype=np.float64))
        found = []
        for label, binding in self.edges():
            width = max(EDGE_FRACTION * binding, EDGE_MIN_EV)
            near = np.abs(energy - binding) <= width
            if near.any():
                found.append(f"{label} at {binding:g} eV ({HC_EV_NM / binding:.4g} nm)")
        delayed = []
        for el, _ in self.composition:
            for level, binding, _ in el.shells:
                band = (energy > binding) & (energy < DELAYED_REACH * binding)
                if level in DELAYED_LEVELS and band.any():
                    delayed.append(f"{el.symbol} {level} at {binding:g} eV")
        out = []
        if found:
            out.append(
                f"{self.name} is the typical-material model, and this scan comes "
                f"within {100 * EDGE_FRACTION:g}% (or {EDGE_MIN_EV:g} eV) of "
                f"absorption edge(s): {'; '.join(found)}. The model has no edge "
                "structure, so its index there is qualitative; supply a "
                "measured table for this range"
            )
        if delayed:
            out.append(
                f"{self.name} is the typical-material model, and this scan lies "
                f"within {DELAYED_REACH:g} times the binding energy above 4d or 4f "
                f"core level(s): {'; '.join(delayed)}. Their absorption rises "
                "over that band instead of at the edge (the delayed onset), which "
                "the model does not have, so its absorption can be out by a "
                "factor of two or more there, and its decrement by tens of percent"
            )
        if (energy < VALENCE_EV).any():
            out.append(
                f"{self.name} is the typical-material model, and this scan goes "
                f"below {VALENCE_EV:g} eV ({HC_EV_NM / VALENCE_EV:.4g} nm), where "
                "the valence electrons the model treats as free decide the "
                "index; read it as qualitative there"
            )
        return tuple(out)


def _index(composition, density, wavelength_nm):
    """``(decrement, absorption)`` for a composition at a density in g/cm^3."""
    energy = HC_EV_NM / wavelength_nm
    molar_mass = sum(el.weight * n for el, n in composition)
    formula_units = density * 1e3 / (molar_mass * 1e-3) * N_A  # per m^3
    f1 = np.zeros_like(energy)
    f2 = np.zeros_like(energy)
    for el, n in composition:
        a, b = scattering_factors(el, energy)
        f1 += n * a
        f2 += n * b
    scale = R_E * (wavelength_nm * 1e-9) ** 2 * formula_units / (2 * np.pi)
    return scale * f1, scale * f2


def typical(formula: str, density: float | None = None) -> TypicalConstants:
    """The model for ``formula`` (``"Au"``, ``"SiO2"``) at ``density`` g/cm^3.

    ``density`` defaults to :data:`DENSITIES` and is required for anything not
    listed there, since a compound's density is not derivable from its atoms.
    """
    composition = _parse(formula)
    if density is None:
        if formula not in DENSITIES:
            raise ValueError(
                f"no default density for {formula!r}; pass density in g/cm^3 "
                f"(defaults exist for {', '.join(DENSITIES)})"
            )
        density = DENSITIES[formula]
    if not density > 0:
        raise ValueError(f"density must be positive, got {density!r}")
    # The ends set exactly, so the range is the documented one to the last bit.
    energy = np.geomspace(LOW_EV, HIGH_EV, 1201)
    energy[[0, -1]] = LOW_EV, HIGH_EV
    wavelength = HC_EV_NM / energy[::-1]
    decrement, absorption = _index(composition, density, wavelength)
    name = f"{PREFIX}{formula}" + ("" if density == DENSITIES.get(formula) else f"@{density:g}")
    return TypicalConstants(
        name=name,
        wavelength_nm=wavelength,
        decrement=decrement,
        absorption=absorption,
        source=(
            f"typical-material model, {formula} at {density:g} g/cm^3 "
            "(Kramers' law per core shell, Kramers-Kronig, free valence "
            "electrons); approximate, see gratinglab.materials.typical"
        ),
        formula=formula,
        density=float(density),
        composition=composition,
    )


def from_name(name: str) -> TypicalConstants:
    """A ``coating`` string, ``typical:<formula>[@<density>]``, as a model."""
    body = name[len(PREFIX):]
    formula, _, rho = body.partition("@")
    try:
        density = float(rho) if rho else None
    except ValueError:
        raise ValueError(f"{name!r}: the density after @ must be a number in g/cm^3")
    return typical(formula, density)
