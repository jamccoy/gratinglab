r"""Optical constants, and where they came from.

In the soft X-ray the refractive index is written

.. math:: n(\lambda) = 1 - \text{decrement} + i\,\text{absorption}

with the ``exp(-iwt)`` time convention fixed by ``docs/conventions.md`` §2, which
is what puts the ``+i`` there. Both parts are small: for Au at 1 nm the
decrement is ~1e-3, so ``n`` differs from 1 in the fourth decimal place and
total external reflection below the critical angle is the only reason grazing
optics work at all.

Why ``decrement`` and ``absorption`` and not the usual symbols
=============================================================

The literature calls these :math:`\delta` and :math:`\beta`. This codebase
already uses :math:`\delta` for the **blaze angle** and :math:`\beta` for the
**diffracted-order angle**, everywhere, and ``conventions.md`` §8 exists
because exactly this kind of collision has bitten the two primary references
already (they use ``d`` for opposite things).

`fresnel.reflectivity` will be called within a few lines of
`geometry.facet_graze(gamma, blaze_angle, alpha)`. Two meanings of
:math:`\delta` in one call stack is how a sign error gets written and not
noticed, so the names here are unambiguous words rather than the conventional
letters. §8 carries the mapping.

Range is enforced, not extrapolated
===================================

A table covers what it covers -- the shipped Au data runs 200 eV to 2 keV. Ask
for 20 nm and this raises, because a linearly extrapolated optical constant is
a plausible number with nothing behind it, which is the failure mode this
project is arranged against. `EfficiencyScan` refuses to carry efficiency on a
non-propagating order for the same reason.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "HC_EV_NM",
    "OpticalConstants",
    "read_cxro",
    "read_ari",
    "write_ari",
]

#: Planck constant times the speed of light, in eV nm, so that
#: ``wavelength_nm = HC_EV_NM / energy_eV``. The value the prototype used
#: (``panter1.py``), kept identical so a re-conversion of the corpus tables
#: reproduces them.
HC_EV_NM = 1239.84193


@dataclass(frozen=True, slots=True)
class OpticalConstants:
    """A tabulated index, and a record of where it came from.

    Attributes
    ----------
    name
        Material identifier, e.g. ``"Au"``.
    wavelength_nm
        Ascending. The table's own sampling; no resampling is done on load.
    decrement, absorption
        :math:`\\delta` and :math:`\\beta` in the usual notation --  see the
        module docstring for why they are not called that here.
    source
        Free text naming the file and its provenance. Ends up on
        ``Provenance.notes`` so a result can say which table produced it.
    """

    name: str
    wavelength_nm: NDArray[np.float64]
    decrement: NDArray[np.float64]
    absorption: NDArray[np.float64]
    source: str = "unknown"

    def __post_init__(self) -> None:
        for field in ("wavelength_nm", "decrement", "absorption"):
            values = np.asarray(getattr(self, field), dtype=np.float64)
            values.setflags(write=False)
            object.__setattr__(self, field, values)

        if self.wavelength_nm.ndim != 1:
            raise ValueError("wavelength_nm must be one-dimensional")
        if len(self.wavelength_nm) < 2:
            raise ValueError(
                f"{self.name}: a table needs at least two points to interpolate "
                f"between, got {len(self.wavelength_nm)}"
            )
        if len(self.decrement) != len(self.wavelength_nm) or len(
            self.absorption
        ) != len(self.wavelength_nm):
            raise ValueError(
                f"{self.name}: decrement and absorption must match "
                f"wavelength_nm in length"
            )
        if not np.all(np.diff(self.wavelength_nm) > 0):
            raise ValueError(
                f"{self.name}: wavelength_nm must be strictly ascending. CXRO "
                "tables are in ascending *energy*, which is descending "
                "wavelength -- reverse on load, as `read_cxro` does."
            )
        if (self.wavelength_nm <= 0).any():
            raise ValueError(f"{self.name}: wavelengths must be positive")

    # -- the range, stated -------------------------------------------------

    @property
    def range_nm(self) -> tuple[float, float]:
        return float(self.wavelength_nm[0]), float(self.wavelength_nm[-1])

    def covers(self, wavelength_nm: ArrayLike) -> NDArray[np.bool_]:
        """Which of these wavelengths the table can answer for."""
        low, high = self.range_nm
        values = np.asarray(wavelength_nm, dtype=np.float64)
        return (values >= low) & (values <= high)

    # -- lookup ------------------------------------------------------------

    def n(self, wavelength_nm: ArrayLike) -> NDArray[np.complex128]:
        r"""Complex index :math:`1 - \text{decrement} + i\,\text{absorption}`."""
        d, b = self.at(wavelength_nm)
        return (1.0 - d) + 1j * b

    def at(
        self, wavelength_nm: ArrayLike
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """``(decrement, absorption)``, linearly interpolated.

        Linear rather than spline: the shipped tables are sampled ~500 points
        per decade, so the interpolation error is far below the tabulated
        precision, and a spline can overshoot near an absorption edge -- which
        would invent structure the data does not have.
        """
        values = np.asarray(wavelength_nm, dtype=np.float64)
        self._require_in_range(values)
        return (
            np.interp(values, self.wavelength_nm, self.decrement),
            np.interp(values, self.wavelength_nm, self.absorption),
        )

    def critical_angle(self, wavelength_nm: ArrayLike) -> NDArray[np.float64]:
        r"""Total-external-reflection critical graze angle, radians.

        :math:`\theta_c \approx \sqrt{2\,\text{decrement}}`, measured from the
        **surface** -- the same sense as `geometry.facet_graze`, so the two are
        directly comparable. That comparison is the validity guard
        ``docs/theory/scalar.md`` §7 has carried as "needs a materials layer to
        evaluate".
        """
        decrement, _ = self.at(wavelength_nm)
        return np.sqrt(2.0 * decrement)

    def _require_in_range(self, values: NDArray[np.float64]) -> None:
        outside = ~self.covers(values)
        if not outside.any():
            return
        low, high = self.range_nm
        worst = values[outside]
        raise ValueError(
            f"{self.name} is tabulated over {low:.4g}-{high:.4g} nm; asked for "
            f"{worst.min():.4g}"
            + (f" to {worst.max():.4g}" if worst.min() != worst.max() else "")
            + " nm. Extrapolating an optical constant would produce a "
            "plausible number with nothing behind it -- supply a table that "
            "covers the scan, or shorten the scan."
        )

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        low, high = self.range_nm
        return (
            f"OpticalConstants({self.name!r}, {len(self.wavelength_nm)} points, "
            f"{low:.4g}-{high:.4g} nm)"
        )


# -- readers ---------------------------------------------------------------


def read_cxro(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read a CXRO/Henke ``(energy, delta, beta)`` export.

    The format, as downloaded from ``henke.lbl.gov``::

         Au Density=19.32
         Energy(eV), Delta, Beta
          200.  0.0112813823  0.00951793324

    **Two** header lines, then whitespace-separated columns in ascending
    energy -- so descending wavelength, which is reversed here.

    A port of ``CXRO_to_n_k`` from the prototype
    (``~/Documents/diffraction_efficiency/panter1.py``), with one correction.
    The prototype hardcodes ``skip_header=3`` against a file with two header
    lines, so it silently discards the first data row -- the lowest energy,
    which is the *longest* wavelength. On the Au table that truncates the range
    from 6.199 nm to 6.171 nm, and the ``.ari`` files in the corpus carry the
    truncation because they were written by that code.

    So the header is detected rather than counted: a line is data when its
    first three whitespace-separated fields all parse as floats. That is
    robust to the two- and three-line variants both, and to an export format
    that grows a line -- and no row is dropped for being at an end of the
    table, which is where a range guard will later be at its most sensitive.
    """
    path = Path(path)
    rows = [
        parsed
        for line in path.read_text().splitlines()
        if (parsed := _numeric_row(line)) is not None
    ]
    if len(rows) < 2:
        raise ValueError(
            f"{path}: expected rows of three numbers "
            f"(energy_eV, delta, beta), found {len(rows)}"
        )
    data = np.asarray(rows, dtype=np.float64)

    energy_ev, decrement, absorption = data[:, 0], data[:, 1], data[:, 2]
    if (energy_ev <= 0).any():
        raise ValueError(f"{path}: energies must be positive")

    order = np.argsort(HC_EV_NM / energy_ev)
    return OpticalConstants(
        name=name or _name_from(path),
        wavelength_nm=(HC_EV_NM / energy_ev)[order],
        decrement=decrement[order],
        absorption=absorption[order],
        source=f"CXRO/Henke export {path.name}",
    )


def read_ari(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read a PCGrate ``.ari`` table.

    Three whitespace-separated columns, **already in ascending wavelength**::

        0.619921 0.000510 0.000106

    Worth stating because the name misleads: the prototype's writer is called
    ``CXRO_to_n_k`` and its output files are named ``*_optical_constants_n_k``,
    but the columns are ``(wavelength_nm, decrement, absorption)`` -- *not*
    ``n`` and ``k``. The conversion ``n = 1 - delta`` appears only in a comment
    there. Reading these as ``n`` and ``k`` would put the index at ~5e-4
    instead of ~1.
    """
    path = Path(path)
    data = np.genfromtxt(path)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(
            f"{path}: expected three columns "
            f"(wavelength_nm, decrement, absorption), got shape {data.shape}"
        )
    return OpticalConstants(
        name=name or _name_from(path),
        wavelength_nm=data[:, 0],
        decrement=data[:, 1],
        absorption=data[:, 2],
        source=f"PCGrate .ari table {path.name}",
    )


def write_ari(constants: OpticalConstants, path: "str | Path") -> Path:
    """Write a ``.ari`` PCGrate can consume.

    ``%f`` and a plain space, matching what the prototype produced, so a file
    written here is byte-comparable with the corpus ones. Note that ``%f`` is
    six decimal places: absorption at short wavelengths rounds to zero, which
    is a property of the format rather than of this writer.
    """
    path = Path(path)
    np.savetxt(
        path,
        np.column_stack(
            [constants.wavelength_nm, constants.decrement, constants.absorption]
        ),
        delimiter=" ",
        fmt="%f",
    )
    return path


def _numeric_row(line: str) -> "list[float] | None":
    """The first three fields as floats, or ``None`` if this is not a data row.

    Deliberately not a regex: "does it parse as a float" is exactly the
    question, and `float` answers it for every spelling a table might use --
    ``200.``, ``1e-3``, and the ``3.55e-006`` three-digit exponents the
    Windows-written ``.ari`` files carry.
    """
    fields = line.split()
    if len(fields) < 3:
        return None
    try:
        return [float(x) for x in fields[:3]]
    except ValueError:
        return None


def _name_from(path: Path) -> str:
    """``Au_CXRO_SXR.txt`` -> ``Au``. The element is what precedes the first
    underscore in every corpus file; a caller who disagrees passes ``name``."""
    return path.stem.split("_")[0]
