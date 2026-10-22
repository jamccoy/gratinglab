r"""Reading refractiveindex.info's YAML tables, ultraviolet to infrared.

The refractiveindex.info database (Polyanskiy, *Sci. Data* **11**, 94 (2024))
is released under CC0, so unlike CXRO's tables there is no question about
redistributing it. One file is one dataset: its ``DATA`` list holds a
tabulated ``n`` and ``k`` together, or each separately, or ``n`` as one of nine
dispersion formulas, often with a tabulated ``k`` beside it. Wavelengths are
micrometres; :class:`~.optical.OpticalConstants` is nanometres, and the
conversion happens here and nowhere else.

What becomes of each entry:

``tabulated nk``, ``tabulated n``, ``tabulated k``
    Rows of ``wavelength n k``. Sorted on read, since 17 blocks in the
    database are not. A wavelength given twice, which happens where two
    measurement ranges were joined and the overlap row differs in the fourth
    digit, is averaged, and the table's ``source`` says how many were.
``formula 1`` to ``formula 9``
    Evaluated exactly as the database documents them (``doc/Dispersion
    formulas.pdf``), then sampled across the entry's ``wavelength_range``,
    because a table is what every solver here consumes. The sampling is
    adaptive: an interval is halved until linear interpolation across it misses
    the formula by under :data:`FORMULA_TOLERANCE` in ``n`` at its midpoint. A
    fixed grid cannot do this, since a range often ends beside a resonance:
    2,001 log-spaced points left a 0.27 error in Xe's index next to its UV
    pole while the median dataset was good to 5e-8.
``tabulated n2``
    **Refused.** That is the nonlinear (Kerr) index :math:`n_2`, from the
    database's separate n2 catalogue, in the same file shape. Read as ``n`` it
    would be a plausible-looking material with nothing right about it.

``n`` with no ``k`` means the dataset is transparent where it is defined, and
``k`` is taken as zero; the ``source`` says so, because "lossless" and "loss
not measured" are different claims and only the file knows which it meant. A
file with ``k`` and no ``n`` is refused, as is one giving ``n`` twice. Where
``n`` and ``k`` cover different ranges, the table covers their overlap only.
A negative ``k`` is clipped to zero and the ``source`` records how many and how
low, and a formula whose stated range crosses its own pole is refused.

Measured on the whole database on 2026-10-21 (3,720 files): 3,435 read; 240
refused as ``n2``; 19 for having no ``n``; 5 for formulas crossing a pole
(CS2/Chemnitz and four GaSe sets); 3 for ``n`` and ``k`` ranges that do not
overlap; 2 single-point tables; 1 giving ``n`` twice.

Optical constants here are :math:`n = n' + ik` with :math:`k \ge 0`, the
``exp(-i\omega t)`` convention of ``conventions.md``, which is refractiveindex.info's
own. They are stored as ``decrement = 1 - n'`` and ``absorption = k``, so a
glass in the visible has a *negative* decrement. That is the same quantity the
X-ray tables carry, merely far from zero.

PyYAML is needed to read these, and only these: ``pip install
"gratinglab[refractiveindex]"``. It is imported when a file is read, never at
package import, so the core install stays numpy, scipy and pydantic.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from .optical import OpticalConstants

__all__ = ["FORMULA_TOLERANCE", "SUFFIXES", "read_refractiveindex"]

#: File suffixes this reader takes, for the library search.
SUFFIXES = (".yml", ".yaml")

#: Largest error in ``n`` that linear interpolation between a formula's samples
#: may leave, at an interval's midpoint. Below the sixth significant figure the
#: coefficients are typically published to.
FORMULA_TOLERANCE = 1e-7

#: Where sampling starts (log-spaced, since ranges span decades) and the most
#: it may grow to before giving up rather than looping on a singularity.
_FORMULA_START = 257
_FORMULA_MAX = 200_001

_UM_TO_NM = 1000.0


def read_refractiveindex(path: "str | Path", name: str | None = None) -> OpticalConstants:
    """Read one refractiveindex.info YAML file.

    ``name`` defaults to the file's stem. In the database a file is named for
    its author (``main/Au/nk/Johnson.yml``), so save it into a library under
    the material's name, ``Au.yml``, or pass ``name``.
    """
    from endstation.optional import require_extra

    require_extra(
        "yaml",
        "refractiveindex",
        package="gratinglab",
        still_works=(
            "CXRO and PCGrate tables, and everything else in gratinglab, work "
            "without it."
        ),
    )
    import yaml

    path = Path(path)
    document = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(document, dict) or not isinstance(document.get("DATA"), list):
        raise ValueError(f"{path}: not a refractiveindex.info file (no DATA list)")

    n_part: tuple[NDArray[np.float64], NDArray[np.float64], str] | None = None
    k_part: tuple[NDArray[np.float64], NDArray[np.float64]] | None = None
    averaged = 0

    def give_n(wl, values, how):
        nonlocal n_part
        if n_part is not None:
            raise ValueError(
                f"{path}: gives n twice ({n_part[2]} and {how}); which one is "
                "meant is not something a reader should guess"
            )
        n_part = (wl, values, how)

    def give_k(wl, values):
        nonlocal k_part
        if k_part is not None:
            raise ValueError(f"{path}: gives k twice")
        k_part = (wl, values)

    for entry in document["DATA"]:
        kind = str(entry.get("type", "")).strip()
        if kind == "tabulated n2":
            raise ValueError(
                f"{path}: 'tabulated n2' is the nonlinear (Kerr) index n2, not "
                "a refractive index; this file is from refractiveindex.info's "
                "n2 catalogue and cannot be a coating"
            )
        if kind in ("tabulated nk", "tabulated n", "tabulated k"):
            columns = {"tabulated nk": 3, "tabulated n": 2, "tabulated k": 2}[kind]
            table, repeats = _rows(entry.get("data", ""), columns, path, kind)
            averaged += repeats
            if kind == "tabulated nk":
                give_n(table[:, 0], table[:, 1], kind)
                give_k(table[:, 0], table[:, 2])
            elif kind == "tabulated n":
                give_n(table[:, 0], table[:, 1], kind)
            else:
                give_k(table[:, 0], table[:, 1])
        elif kind.startswith("formula "):
            number = kind.split()[1]
            low, high = (float(v) for v in str(entry["wavelength_range"]).split())
            coefficients = [float(v) for v in str(entry["coefficients"]).split()]
            wl, values = _sample(number, coefficients, low, high, path)
            give_n(wl, values, kind)
        else:
            raise ValueError(f"{path}: unrecognised DATA type {kind!r}")

    if n_part is None:
        raise ValueError(
            f"{path}: gives k but no n, so there is no index to compute with"
        )
    wl_n, n, how = n_part
    notes = [how]

    if k_part is None:
        wavelengths, real, k = wl_n, n, np.zeros_like(n)
        notes.append("k not given, taken as 0")
    else:
        wl_k, k_values = k_part
        low = max(wl_n[0], wl_k[0])
        high = min(wl_n[-1], wl_k[-1])
        if low >= high:
            raise ValueError(
                f"{path}: n covers {wl_n[0]:g}-{wl_n[-1]:g} um and k covers "
                f"{wl_k[0]:g}-{wl_k[-1]:g} um, which do not overlap"
            )
        grid = np.union1d(wl_n, wl_k)
        wavelengths = grid[(grid >= low) & (grid <= high)]
        real = np.interp(wavelengths, wl_n, n)
        k = np.interp(wavelengths, wl_k, k_values)
        if k_part is not None and how != "tabulated nk":
            notes.append("tabulated k")

    if averaged:
        notes.append(f"averaged {averaged} repeated wavelength(s)")
    if not (np.isfinite(real).all() and np.isfinite(k).all()):
        raise ValueError(f"{path}: the index is not finite everywhere in its range")
    if (k < 0).any():
        # 15 files in the database carry some. Most are round-off in a
        # transparent region (-1e-17 to -1e-5); a few are real excursions in
        # Kramers-Kronig-derived data, to -0.15. Either way, read literally it
        # is gain. Clipped rather than refused, since refusing would reject
        # e.g. Crystran's fused silica over a -2e-9, and recorded rather than
        # silent, since the large ones are a property of the data a reader of
        # the result should be able to see.
        notes.append(f"clipped {(k < 0).sum()} negative k to 0, the lowest {k.min():.3g}")
        k = np.clip(k, 0.0, None)

    reference = _plain(str(document.get("REFERENCES", "")))
    source = f"refractiveindex.info {path.name} ({', '.join(notes)})"
    if reference:
        source += f": {reference}"
    return OpticalConstants(
        name=name or path.stem,
        wavelength_nm=wavelengths * _UM_TO_NM,
        decrement=1.0 - real,
        absorption=k,
        source=source,
    )


def dispersion(
    number: str, c: list[float], wl: NDArray[np.float64], path: "str | Path" = "?"
) -> NDArray[np.float64]:
    """Real index from one of the database's nine formulas, wavelength in um.

    Written from ``doc/Dispersion formulas.pdf``. Missing trailing
    coefficients are zero, as the database's own tooling treats them.
    """
    C = [0.0] + list(c) + [0.0] * (17 - len(c))  # C[1]..C[17], 1-based
    w2 = wl**2
    if number == "1":  # Sellmeier
        n2 = 1 + C[1] + sum(C[i] * w2 / (w2 - C[i + 1] ** 2) for i in range(2, 17, 2))
    elif number == "2":  # Sellmeier-2
        n2 = 1 + C[1] + sum(C[i] * w2 / (w2 - C[i + 1]) for i in range(2, 17, 2))
    elif number == "3":  # polynomial
        n2 = C[1] + sum(C[i] * wl ** C[i + 1] for i in range(2, 17, 2))
    elif number == "4":  # RefractiveIndex.INFO
        n2 = (
            C[1]
            + C[2] * wl ** C[3] / (w2 - C[4] ** C[5])
            + C[6] * wl ** C[7] / (w2 - C[8] ** C[9])
            + sum(C[i] * wl ** C[i + 1] for i in range(10, 17, 2))
        )
    elif number == "5":  # Cauchy
        return C[1] + sum(C[i] * wl ** C[i + 1] for i in range(2, 11, 2))
    elif number == "6":  # gases
        return 1 + C[1] + sum(C[i] / (C[i + 1] - wl**-2.0) for i in range(2, 11, 2))
    elif number == "7":  # Herzberger
        L = 1 / (w2 - 0.028)
        return C[1] + C[2] * L + C[3] * L**2 + C[4] * w2 + C[5] * w2**2 + C[6] * w2**3
    elif number == "8":  # Retro: (n^2 - 1) / (n^2 + 2)
        r = C[1] + C[2] * w2 / (w2 - C[3]) + C[4] * w2
        n2 = (1 + 2 * r) / (1 - r)
    elif number == "9":  # Exotic
        n2 = C[1] + C[2] / (w2 - C[3]) + C[4] * (wl - C[5]) / ((wl - C[5]) ** 2 + C[6])
    else:
        raise ValueError(f"{path}: unknown dispersion formula {number!r}")
    if (n2 < 0).any():
        raise ValueError(
            f"{path}: formula {number} gives n^2 < 0 inside its own stated range"
        )
    return np.sqrt(n2)


def _sample(
    number: str, c: list[float], low: float, high: float, path: "str | Path"
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """A formula sampled so linear interpolation stays within tolerance."""
    wl = np.geomspace(low, high, _FORMULA_START)
    n = dispersion(number, c, wl, path)
    while True:
        mid = np.sqrt(wl[:-1] * wl[1:])  # geometric, matching the log spacing
        exact = dispersion(number, c, mid, path)
        # Linear in wavelength between neighbours, which is exactly what
        # `OpticalConstants.at` does with the finished table.
        t = (mid - wl[:-1]) / (wl[1:] - wl[:-1])
        linear = n[:-1] + t * (n[1:] - n[:-1])
        bad = np.abs(exact - linear) > FORMULA_TOLERANCE
        if not bad.any():
            return wl, n
        if len(wl) + bad.sum() > _FORMULA_MAX:
            raise ValueError(
                f"{path}: formula {number} still varies too fast to tabulate "
                f"within {FORMULA_TOLERANCE:g} after {_FORMULA_MAX} points; "
                "its stated range likely contains a resonance pole"
            )
        wl = np.insert(wl, np.flatnonzero(bad) + 1, mid[bad])
        n = np.insert(n, np.flatnonzero(bad) + 1, exact[bad])


def _rows(
    text: str, columns: int, path: Path, kind: str
) -> tuple[NDArray[np.float64], int]:
    """Parse a tabulated block: sorted, repeats averaged. Returns the repeat count."""
    rows = [line.split() for line in str(text).splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"{path}: '{kind}' has no rows")
    if any(len(r) < columns for r in rows):
        raise ValueError(f"{path}: '{kind}' rows need {columns} columns")
    table = np.array([[float(v) for v in r[:columns]] for r in rows])
    if (table[:, 0] <= 0).any():
        raise ValueError(f"{path}: wavelengths must be positive")
    table = table[np.argsort(table[:, 0], kind="stable")]
    unique, inverse, counts = np.unique(table[:, 0], return_inverse=True, return_counts=True)
    if len(unique) == len(table):
        return table, 0
    merged = np.zeros((len(unique), columns))
    np.add.at(merged, inverse, table)
    merged /= counts[:, None]
    return merged, int((counts > 1).sum())


def _plain(html: str) -> str:
    """A reference as one line of text: tags dropped, whitespace collapsed."""
    text = re.sub(r"<[^>]+>", "", html)
    text = re.sub(r"\s+", " ", text).strip()
    return text if len(text) <= 160 else text[:157].rstrip() + "..."
