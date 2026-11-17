r"""Reading a plain table of ``wavelength, n, k``, such as your own measurement.

```text
# Ni film, ellipsometry, 2026-11-16
wavelength (nm), n, k
250, 1.62, 2.10
260, 1.63, 2.19
```

For an index from anywhere that is not CXRO, refractiveindex.info or PCGrate:
an ellipsometry export, a column pair from a paper, a spreadsheet. Saved in a
library as ``<Name>.nk`` or ``<Name>.csv`` it is ``coating="<Name>"``.

The format
==========

Rows of two or three numbers, separated by commas, semicolons, tabs or spaces.
Lines starting ``#`` are comments, and any other line that is not numbers is a
header. The first column is wavelength or photon energy, in the unit the
header line just above the data names in its first field, so one of:

========================  =========================
header names              unit
========================  =========================
``nm``                    nanometres
``um``, ``µm``, micron    micrometres
``A``, ``Å``, angstrom    ångströms
``eV``                    photon energy, electronvolts
========================  =========================

``wavelength (nm)``, ``lambda_um`` and ``Energy [eV]`` all qualify. A file
with no such header is refused rather than read as nanometres, because a
table in micrometres read as nanometres is a plausible material a thousand
times too blue, and nothing downstream could tell.

With two columns ``k`` is zero, and the ``source`` says it was not given
rather than measured as zero. Rows are sorted, so an energy scan reads as
well as a wavelength scan. A wavelength given twice is refused, naming both
lines. A negative ``k``, which ellipsometry fits produce as noise around a
transparent region, is clipped to zero and the ``source`` records how many
and how low, as :mod:`.refractiveindex` does.

The index is :math:`n' + ik` with :math:`k \ge 0` (``docs/conventions.md``),
stored as ``decrement = 1 - n'`` and ``absorption = k``.

Read by :func:`gratinglab.materials.lookup` through its suffix table, and by
:func:`gratinglab.materials.from_nk_file` for a file outside any library.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .optical import HC_EV_NM, OpticalConstants

__all__ = ["SUFFIXES", "read_nk"]

#: File suffixes this reader takes, for the library search.
SUFFIXES = (".nk", ".csv")

_SPLIT = re.compile(r"[,;\s]+")

#: Header word to (factor to nm, is energy). Matched as whole words in the
#: header's first field, after brackets and underscores become spaces.
_UNITS = {
    "nm": (1.0, False),
    "um": (1000.0, False),
    "µm": (1000.0, False),
    "μm": (1000.0, False),  # Greek mu, which is a different code point
    "micron": (1000.0, False),
    "microns": (1000.0, False),
    "a": (0.1, False),
    "å": (0.1, False),
    "angstrom": (0.1, False),
    "angstroms": (0.1, False),
    "ev": (None, True),
}


def _fields(line: str) -> list[str]:
    return [f for f in _SPLIT.split(line.strip()) if f]


def _numbers(line: str) -> "list[float] | None":
    """Two or three floats, or ``None`` when this is not a data row."""
    fields = _fields(line)
    if not 2 <= len(fields) <= 3:
        return None
    try:
        return [float(f) for f in fields]
    except ValueError:
        return None


def _unit(header: str, path: Path) -> "tuple[float | None, bool]":
    """The first column's unit, from the first field of the header line."""
    first = re.split(r"[,;\t]|\s{2,}", header.strip().lstrip("#").strip())[0]
    words = re.split(r"[\s()\[\]{}_/]+", first.lower())
    found = [_UNITS[w] for w in words if w in _UNITS]
    if len(found) != 1:
        raise ValueError(
            f"{path}: the header line above the data must name the first "
            "column's unit, one of nm, um, A or eV (for example "
            f"'wavelength (nm), n, k'); found {header.strip()!r}"
        )
    return found[0]


def read_nk(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read a ``wavelength, n[, k]`` table; ``name`` defaults to the file stem."""
    path = Path(path)
    header = None
    rows: list[tuple[int, list[float]]] = []
    for number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = _numbers(line)
        if values is None:
            if rows:
                raise ValueError(f"{path}, line {number}: not a row of numbers: {line.strip()!r}")
            header = line
        else:
            if rows and len(values) != len(rows[0][1]):
                raise ValueError(f"{path}, line {number}: {len(values)} columns where the first row has {len(rows[0][1])}")
            rows.append((number, values))
    if header is None:
        raise ValueError(f"{path}: no header line naming the first column's unit (nm, um, A or eV)")
    if len(rows) < 2:
        raise ValueError(f"{path}: expected at least two rows of numbers, found {len(rows)}")
    factor, is_energy = _unit(header, path)

    lines = np.array([n for n, _ in rows])
    data = np.array([v for _, v in rows], dtype=np.float64)
    x = data[:, 0]
    if (x <= 0).any():
        raise ValueError(f"{path}: the first column must be positive")
    wavelength = HC_EV_NM / x if is_energy else x * factor
    order = np.argsort(wavelength, kind="stable")
    wavelength, data, lines = wavelength[order], data[order], lines[order]
    same = np.flatnonzero(np.diff(wavelength) == 0)
    if same.size:
        a, b = lines[same[0]], lines[same[0] + 1]
        raise ValueError(f"{path}: lines {a} and {b} give the same wavelength; keep one")

    n = data[:, 1]
    notes = []
    if data.shape[1] == 3:
        k = data[:, 2].copy()
        negative = k < 0
        if negative.any():
            notes.append(f"{int(negative.sum())} negative k clipped to 0 (lowest {k.min():.3g})")
            k[negative] = 0.0
    else:
        k = np.zeros_like(n)
        notes.append("no k column, so k = 0: not measured, not shown lossless")
    unit = "eV" if is_energy else {1.0: "nm", 1000.0: "um", 0.1: "A"}[factor]
    source = f"n,k table {path.name} (first column in {unit})"
    return OpticalConstants(
        name=name or path.stem,
        wavelength_nm=wavelength,
        decrement=1.0 - n,
        absorption=k,
        source=source + "".join(f"; {note}" for note in notes),
    )
