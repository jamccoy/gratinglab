"""Reader and writer for PCGrate ``.ggp`` polygonal boundary profiles.

A ``.ggp`` holds a groove profile normalised to the period: ``x`` runs 0 to 1
across exactly one period and ``y`` is height as a fraction of that period.
The scale is the one :mod:`gratinglab.profiles` uses, so no rescaling is
needed. The direction is not: ``x`` advances along the dispersion direction
and a profile's ``t`` against it (``docs/conventions.md``, "Direction of the
profile parameter"), so :func:`to_profile` and :func:`write_ggp` reflect the
points through :meth:`~gratinglab.profiles.FromProfileData.along_dispersion`
and its inverse. The PCGrate corpus is what fixes the file's direction: the
TASTE run matches only with ``x`` along ``+d-hat``.

Canonical form::

    3 0 - Polygonal type
    Period: 1 PSC: 1
    0.000000 0.000000
    0.012346 0.000808
    ...

The shared convention allows a scalar calculation and a PCGrate run to be
driven from identical geometry, which isolates the method difference from any
difference in assumed groove shape.

Three header variants exist on disk, and only the first is accepted by PCGrate:

===================================  ================  ==============
Header                               PCGrate accepts   Seen in
===================================  ================  ==============
``3 0 - Polygonal type`` + ``Period`` yes              most files
``# 3 0 - Polygonal type`` (hashed)   **no**           ``*_fix6/7.ggp``
``# X(normalized...)`` only           **no**           ``*_fix{,2,3,4}.ggp``
===================================  ================  ==============

This module reads all three, since they are real files that need to be loaded,
and writes only the canonical form. ``np.savetxt(header=...)`` is what produces
the hashed variants, so it must not be used here.

This is now the only ``.ggp`` writer in the project. The metrology package
previously shipped a second copy, and the two were kept in agreement by a
comment and a byte-compatibility test. :mod:`gratinglab.metrology` now calls
this one.

The format is the reason the sidecar exists: a ``.ggp`` cannot carry the period
in nm. It holds the shape of a grating without its scale, so anything
downstream must obtain the period from another source.
``benchmarks/corpus.toml`` is that other source, recovered by hand.
Code holding a measured profile should reach for
:meth:`gratinglab.metrology.boundary.BoundaryProfile.to_problem` instead, which
carries the period it measured. See ``docs/roadmap.md`` on the native boundary
format for the file-based version of the same fix.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike

from ..profiles import FromProfileData

__all__ = [
    "GgpFile",
    "GGP_HEADER",
    "read_ggp",
    "read_profile",
    "write_ggp",
    "to_profile",
]

#: The only header PCGrate accepts. Never prefix these lines with ``#``.
GGP_HEADER = "3 0 - Polygonal type\nPeriod: 1 PSC: 1"

_CANONICAL_FIRST = "3 0 - Polygonal type"


@dataclass(frozen=True, slots=True)
class GgpFile:
    """A parsed ``.ggp``, with the header provenance kept.

    ``format_valid`` records whether the file as written would
    load in PCGrate. This is surfaced because a hashed header fails there
    while parsing without error here.
    """

    t: np.ndarray
    y: np.ndarray
    header_variant: str
    format_valid: bool
    source: Path

    @property
    def depth(self) -> float:
        """Peak-to-valley height, in units of the period."""
        return float(np.ptp(self.y))


def read_ggp(path: str | Path) -> GgpFile:
    """Parse a ``.ggp``, tolerating all three header variants.

    Raises
    ------
    ValueError
        On a malformed data line. The error is deliberate, because a truncated
        point that was silently dropped would shift the profile and change
        every efficiency computed from it.
    """
    path = Path(path)
    lines = _read_text(path).splitlines()

    header_lines: list[str] = []
    points: list[tuple[float, float]] = []
    first_data_line = None

    for number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped:
            continue

        # A line is data if and only if it is exactly two numbers. Anything
        # else before the first data line is header; anything else after it is
        # an error. Classifying by content rather than by position is what
        # stops a malformed point like "0.25 abc" from being quietly absorbed
        # as a header and dropped.
        pair = _numeric_pair(stripped)
        if pair is not None:
            points.append(pair)
            if first_data_line is None:
                first_data_line = number
            continue

        if first_data_line is not None:
            raise ValueError(
                f"{path}:{number}: expected an 'x y' pair of numbers, got "
                f"{stripped!r}. Header lines must precede all data, and a "
                "truncated point would silently shift the profile."
            )
        header_lines.append(stripped)

    if not points:
        raise ValueError(f"{path}: no data points found")
    if len(points) < 3:
        raise ValueError(f"{path}: only {len(points)} point(s); need at least 3")

    variant, compatible = _classify(header_lines)
    array = np.asarray(points, dtype=np.float64)

    return GgpFile(
        t=array[:, 0],
        y=array[:, 1],
        header_variant=variant,
        format_valid=compatible,
        source=path,
    )


def _read_text(path: Path) -> str:
    """UTF-8, then PCGrate's own cp1252, then latin-1, which cannot fail.

    `write_ggp` writes UTF-8 and PCGrate runs on Windows and writes cp1252, so
    both are real. UTF-8 goes first because valid UTF-8 is almost never valid
    by accident, while cp1252 accepts nearly anything and would silently misread
    a UTF-8 header. This is the same chain as `efficiency_table`, which reads the other
    file PCGrate writes. Before this the platform's locale decided, so a header
    that loaded on Linux could raise on Windows.
    """
    raw = path.read_bytes()
    for encoding in ("utf-8", "cp1252"):
        try:
            return raw.decode(encoding)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1")


def _numeric_pair(line: str) -> tuple[float, float] | None:
    """Return ``(x, y)`` if the line is exactly two numbers, else ``None``."""
    fields = line.split()
    if len(fields) != 2:
        return None
    try:
        return float(fields[0]), float(fields[1])
    except ValueError:
        return None


def _classify(header_lines: list[str]) -> tuple[str, bool]:
    """Name the header variant and say whether PCGrate would accept it."""
    if not header_lines:
        return "none", False
    if header_lines[0] == _CANONICAL_FIRST:
        return "canonical", True
    if header_lines[0].lstrip("# ").startswith("3 0"):
        return "commented-canonical", False
    return "commented-columns", False


def to_profile(ggp: GgpFile) -> FromProfileData:
    """Convert to a :class:`~gratinglab.profiles.FromProfileData`.

    Heights are referenced to the profile minimum so the valley sits at zero,
    matching every analytic profile in :mod:`gratinglab.profiles`. The file's
    ``x`` runs along ``+d-hat``, so the points are reflected into ``t``.
    """
    return FromProfileData.along_dispersion(ggp.t, ggp.y - ggp.y.min())


def read_profile(path: str | Path) -> FromProfileData:
    """Read a ``.ggp`` straight into a profile."""
    return to_profile(read_ggp(path))


def write_ggp(
    path: str | Path,
    profile: FromProfileData | None = None,
    *,
    t: ArrayLike | None = None,
    y: ArrayLike | None = None,
    fmt: str = "%f",
) -> Path:
    """Write a profile in the canonical, PCGrate-loadable form.

    Pass either a ``profile`` or explicit ``t``/``y`` arrays. A profile is
    reflected onto the file's ``+d-hat`` axis, the inverse of
    :func:`to_profile`. Explicit arrays are written as given, so they are
    file coordinates already, as the metrology pipeline's ``x_norm`` is.

    Values are written to six decimals, single-space separated and newline
    terminated. The format is
    pinned by ``tests/test_ggp.py``; PCGrate rejects variations on it.
    """
    if (profile is None) == (t is None):
        raise ValueError("pass either a profile or both t and y, not both or neither")

    if profile is not None:
        if not profile.is_single_valued():
            raise ValueError(
                "cannot write an undercut profile as .ggp: PCGrate's polygonal "
                "border expects a boundary that advances monotonically in x"
            )
        t_values, y_values = profile.dispersion_axis()
    else:
        t_values = np.asarray(t, dtype=np.float64)
        y_values = np.asarray(y, dtype=np.float64)
        if len(t_values) != len(y_values):
            raise ValueError(
                f"t and y differ in length: {len(t_values)} vs {len(y_values)}"
            )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(GGP_HEADER + "\n")
        for x_value, height in zip(t_values, y_values):
            handle.write(f"{fmt % x_value} {fmt % height}\n")
    return path
