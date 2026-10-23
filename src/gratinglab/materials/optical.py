r"""Optical constants, and where they came from.

In the soft X-ray the refractive index is written

.. math:: n(\lambda) = 1 - \text{decrement} + i\,\text{absorption}

with the ``exp(-iwt)`` time convention fixed by ``docs/conventions.md`` §2, which
is what puts the ``+i`` there. Both parts are small: for Au at 1 nm the
decrement is ~1e-3, so ``n`` differs from 1 in the fourth decimal place and
total external reflection below the critical angle is the only reason grazing
optics work at all.

Why the names ``decrement`` and ``absorption`` are used
=======================================================

The literature calls these :math:`\delta` and :math:`\beta`. This codebase
already uses :math:`\delta` for the blaze angle and :math:`\beta` for the
diffracted-order angle, everywhere, and ``conventions.md`` §8 exists because
this kind of collision has already caused trouble in the two primary
references (they use ``d`` for opposite things).

`fresnel.reflectivity` will be called within a few lines of
`geometry.facet_graze(gamma, blaze_angle, alpha)`. Two meanings of
:math:`\delta` in one call stack make it easy to write a sign error without
noticing it, so the names here are unambiguous words in place of the
conventional letters. §8 carries the mapping.

The tabulated range is enforced
===============================

A table is valid only over the range it covers; the shipped Au data runs
200 eV to 2 keV. A request for 20 nm raises an error, because a linearly
extrapolated optical constant is a plausible number with nothing behind it,
which is the failure mode this project is arranged to prevent. `EfficiencyScan`
refuses to carry efficiency on a non-propagating order for the same reason.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = [
    "HC_EV_NM",
    "NEAR_UNITY",
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

#: How far from 1 the index may sit, as :math:`|n - 1|`, and still be in the
#: regime where a critical angle means something. Every X-ray and EUV table
#: this project has met sits well inside it (the shipped Au peaks at 0.015, and
#: Mo or Ru at 13.5 nm are near 0.1); ordinary UV-to-IR materials sit well
#: outside it (a glass at 0.3 to 1, a visible metal at several). It is also the
#: regime where reflectivity has nothing else to stand on: at
#: :math:`|n - 1| = 0.2` normal-incidence reflectivity is about 1%, so above
#: the critical angle it genuinely collapses, which is what makes losing total
#: external reflection worth a warning at all.
NEAR_UNITY = 0.2


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
        :math:`\\delta` and :math:`\\beta` in the usual notation. See the
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

        Linear interpolation is used in preference to a spline. The shipped
        tables are sampled ~500 points per decade, so the interpolation error
        is far below the tabulated precision, and a spline can overshoot near
        an absorption edge, which would invent structure the data does not
        have.
        """
        values = np.asarray(wavelength_nm, dtype=np.float64)
        self._require_in_range(values)
        return (
            np.interp(values, self.wavelength_nm, self.decrement),
            np.interp(values, self.wavelength_nm, self.absorption),
        )

    def critical_angle(self, wavelength_nm: ArrayLike) -> NDArray[np.float64]:
        r"""Total-external-reflection critical graze angle, radians.

        :math:`\theta_c = \arccos(1 - \text{decrement})`, measured from the
        surface, in the same sense as `geometry.facet_graze`, so the two are
        directly comparable. That comparison is the validity guard that
        ``docs/theory/scalar.md`` §7 has carried as "needs a materials layer to
        evaluate".

        This is the exact lossless edge, in place of the small-angle one.
        Snell's law from the surface, :math:`\cos\theta = n'\cos\theta_t`, has
        no real transmitted angle once :math:`\cos\theta > n'`, which gives the
        arccos. :math:`\sqrt{2\,\text{decrement}}` is its leading term, and for
        the :math:`10^{-3}` decrements of the soft X-ray the two differ by a
        relative :math:`\text{decrement}/12`, about 1e-4. The exact form costs
        nothing and remains correct at EUV, where the decrement reaches 0.1.
        Absorption softens the edge without moving it, which is why the
        imaginary part takes no part in the value.

        Zero means there is no critical angle, and this can arise in two
        ways:

        - the real part of n is not below 1. Just below an X-ray absorption
          edge anomalous dispersion drives the decrement negative (Cr near
          574 eV, Ni near 853 eV), and every glass in the visible is there
          anyway. A bare square root returned NaN for the first, and
          ``graze > nan`` is False, so the scalar solver's guard passed
          silently at exactly the wavelengths it should warn about most.
        - the index is not near 1 (see :meth:`near_unity`), which is nearly
          all of the UV to IR. Au at 600 nm has :math:`n \approx 0.25 + 3.1i`,
          real part below 1, and the arccos would say 76 degrees; its Fresnel
          reflectivity stays between 0.90 and 0.93 from 10 degrees to normal
          incidence, so there is no edge to be on either side of.

        Zero is safe to compare against in either case. Whether a graze above
        it warrants a warning is a separate question, which :meth:`near_unity`
        answers, and the scalar solver asks both.
        """
        n = self.n(wavelength_nm)
        has_edge = (n.real < 1.0) & self._near_unity(n)
        # The clip only keeps arccos in its domain where has_edge is False and
        # the value is discarded; near unity it never changes a kept one.
        return np.where(has_edge, np.arccos(np.clip(n.real, -1.0, 1.0)), 0.0)

    def near_unity(self, wavelength_nm: ArrayLike) -> NDArray[np.bool_]:
        r"""Where :math:`|n - 1| \le` ``NEAR_UNITY``: the X-ray-like regime.

        Here reflectivity is only appreciable at grazing incidence, so a facet
        past the critical angle, or a wavelength with no critical angle at all,
        loses its reflectivity and warrants a warning. Outside it, for a glass
        or a metal from the UV to the IR, reflectivity does not depend on
        grazing, a near-normal design is ordinary, and :meth:`critical_angle`
        is zero because there is no edge, without anything being wrong.

        The test is on :math:`|n - 1|` and not on absorption against
        decrement, because the shipped Au table has absorption above its
        decrement over a third of its range, where the edge is soft but
        reflectivity still collapses above it.
        """
        return self._near_unity(self.n(wavelength_nm))

    @staticmethod
    def _near_unity(n: NDArray[np.complex128]) -> NDArray[np.bool_]:
        return np.abs(n - 1.0) <= NEAR_UNITY

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

    The file has two header lines, then whitespace-separated columns in
    ascending energy, and therefore descending wavelength, which is reversed
    here.

    This is a port of ``CXRO_to_n_k`` from the prototype
    (``~/Documents/diffraction_efficiency/panter1.py``), with one correction.
    The prototype hardcodes ``skip_header=3`` against a file with two header
    lines, so it silently discards the first data row: the lowest energy,
    which is the longest wavelength. On the Au table that truncates the range
    from 6.199 nm to 6.171 nm, and the ``.ari`` files in the corpus carry the
    truncation because they were written by that code.

    The header is therefore detected instead of counted: a line is data when
    its first three whitespace-separated fields all parse as floats. This is
    robust to both the two- and three-line variants and to an export format
    that gains a line, and no row is dropped for being at an end of the
    table, which is where a range guard will later be at its most sensitive.
    """
    path = Path(path)
    rows = [
        parsed
        for line in path.read_text(encoding="utf-8").splitlines()
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

    Three whitespace-separated columns, already in ascending wavelength::

        0.619921 0.000510 0.000106

    The column meaning is stated because the name is misleading: the
    prototype's writer is called ``CXRO_to_n_k`` and its output files are named
    ``*_optical_constants_n_k``, but the columns are
    ``(wavelength_nm, decrement, absorption)`` and are not ``n`` and ``k``.
    The conversion ``n = 1 - delta`` appears only in a comment there. Reading
    these as ``n`` and ``k`` would put the index at ~5e-4 instead of ~1.
    """
    path = Path(path)
    data = np.genfromtxt(path, encoding="utf-8")
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

    No regex is used, because "does it parse as a float" is the question
    being asked, and `float` answers it for every spelling a table might use:
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
