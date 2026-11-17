"""Atomic numbers and weights for hydrogen to uranium, the range the
tabulated scattering factors in :mod:`.factors` cover.

Read by :mod:`.factors`, which needs a weight for every element of a formula
to turn a mass density into atoms per volume. :mod:`.typical` keeps its own
weights beside each element's other facts; a test holds the two equal.

Weights are the IUPAC standard atomic weights (CIAAW), given as the
conventional single value where IUPAC states an interval (H, Li, B, C, N, O,
Mg, Si, S, Cl, Ar, Br, Tl). An element with no stable isotope (Tc, Pm, and
Po to Ac) has none, and carries the mass number of its longest-lived isotope
instead, which is what CXRO's calculator does too. Each is a measured
physical constant, not a reproduced table.
"""

from __future__ import annotations

__all__ = ["WEIGHTS", "ATOMIC_NUMBER", "weight"]

#: Symbol to standard atomic weight in g/mol, in order of atomic number.
WEIGHTS: dict[str, float] = {
    "H": 1.008, "He": 4.002602, "Li": 6.94, "Be": 9.0121831,
    "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998403163, "Ne": 20.1797, "Na": 22.98976928, "Mg": 24.305,
    "Al": 26.9815384, "Si": 28.085, "P": 30.973761998, "S": 32.06,
    "Cl": 35.45, "Ar": 39.95, "K": 39.0983, "Ca": 40.078,
    "Sc": 44.955907, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961,
    "Mn": 54.938043, "Fe": 55.845, "Co": 58.933194, "Ni": 58.6934,
    "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630,
    "As": 74.921595, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.4678, "Sr": 87.62, "Y": 88.905838, "Zr": 91.222,
    "Nb": 92.90637, "Mo": 95.95, "Tc": 98.0, "Ru": 101.07,
    "Rh": 102.90549, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.414,
    "In": 114.818, "Sn": 118.710, "Sb": 121.760, "Te": 127.60,
    "I": 126.90447, "Xe": 131.293, "Cs": 132.90545196, "Ba": 137.327,
    "La": 138.90547, "Ce": 140.116, "Pr": 140.90766, "Nd": 144.242,
    "Pm": 145.0, "Sm": 150.36, "Eu": 151.964, "Gd": 157.249,
    "Tb": 158.925354, "Dy": 162.500, "Ho": 164.930329, "Er": 167.259,
    "Tm": 168.934219, "Yb": 173.045, "Lu": 174.9668, "Hf": 178.486,
    "Ta": 180.94788, "W": 183.84, "Re": 186.207, "Os": 190.23,
    "Ir": 192.217, "Pt": 195.084, "Au": 196.966570, "Hg": 200.592,
    "Tl": 204.38, "Pb": 207.2, "Bi": 208.98040, "Po": 209.0,
    "At": 210.0, "Rn": 222.0, "Fr": 223.0, "Ra": 226.0,
    "Ac": 227.0, "Th": 232.0377, "Pa": 231.03588, "U": 238.02891,
}

#: Symbol to atomic number, from the order above.
ATOMIC_NUMBER: dict[str, int] = {s: z for z, s in enumerate(WEIGHTS, start=1)}


def weight(symbol: str) -> float:
    """``symbol``'s atomic weight, refusing one past uranium or misspelt."""
    try:
        return WEIGHTS[symbol]
    except KeyError:
        raise ValueError(
            f"{symbol!r} is not an element from H to U, the range tabulated "
            "scattering factors cover"
        ) from None
