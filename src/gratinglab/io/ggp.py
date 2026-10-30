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

Default form, points only::

    0.000000 0.000000
    0.012346 0.000808
    ...

The shared convention allows a scalar calculation and a PCGrate run to be
driven from identical geometry, which isolates the method difference from any
difference in assumed groove shape.

Five control-string variants exist on disk:

=====================================================  =========================
Control strings                                        Seen in
=====================================================  =========================
none, points only                                      Aditi's, and the default
``3 0 - 20 harmonics - Polygonal type`` + ``Period:``  PCGrate's own exports
``3 0 - Polygonal type`` + ``Period:``                 legacy files written here
``# 3 0 - Polygonal type`` (hashed)                    ``*_fix6/7.ggp``
``# X(normalized...)`` only                            ``*_fix{,2,3,4}.ggp``
=====================================================  =========================

Which of these PCGrate takes is documented, and the vendor's own files agree
with the documentation. The *User's Guide of PCGrate-S(X) Series, version 6.1*
gives the step-by-step instruction for authoring a polygonal ``.ggp`` in
section 10.2.2 (pp. 62-63) and the complete list of format errors the program
rejects in section 10.2.3 (p. 63). Until 2026-10-30 this docstring said the
question had never been tested here.

The first control string is ``3 0 - 20 harmonics - Polygonal type``, verbatim
(10.2.2, step 2). All 35 polygonal files in the installation's own
``BorderProfiles/`` carry exactly that line, with no variation. The second
carries the period, in relative units or as any positive real, and the profile
scaling coefficient: ``Period: 1 PSC: 1``. Until 2026-10-30 this module wrote
and recognised ``3 0 - Polygonal type``, which lacks the ``- 20 harmonics -``
segment and therefore fails the documented check on the first control string;
:data:`GGP_HEADER` is now the documented pair.

A file with no control strings at all, or with control strings in another
format, is explicitly legal input (10.2.3, closing paragraph) -- but through
the Border Profile Editor, which exports it into the program, rather than a
direct load. That is why Aditi's practice works: she loads these through
PCGrate's border-profile import and reports that it does not want a header,
which is why the writer leaves one out unless asked.

The geometry rules are a closed list (10.2.3), and what :func:`write_ggp`
produces already satisfies every one of them:

* exactly two numbers per data string, separated by a space, with any
  fractional part after a decimal point rather than a comma -- at most 16
  digits including the point, and exponential notation is allowed;
* the period within 0.0001 .. 1 000 000, the scaling coefficient within
  0.0001 .. 100;
* every x within 0 .. 1 and never less than its predecessor;
* every y within 0 .. 100;
* the first point exactly ``0 0``, the last exactly ``1 0``;
* between 2 and 5000 data strings, the upper bound ``ggp_n_points`` in
  :mod:`gratinglab.metrology.settings` is validated against.

Twenty-four of the 35 separate their columns with a tab rather than the space
10.2.2 asks for, so the program reads more separators than it documents; this
module reads either. All 35 also end in a line break, which is worth knowing
beside the note in :func:`write_ggp` that a trailing one was reported
rejected.

This module reads all five variants. ``np.savetxt(header=...)`` is what
produces the hashed ones, so it must not be used here.

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

import re
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

#: The documented two control strings, written only on request. Never prefix
#: these lines with ``#``: a hashed control string is one "in another format",
#: which the guide accepts only through the Border Profile Editor.
GGP_HEADER = "3 0 - 20 harmonics - Polygonal type\nPeriod: 1 PSC: 1"

#: The first control string the guide specifies (10.2.2, step 2), and the one
#: every polygonal file in PCGrate's own ``BorderProfiles/`` carries.
_CANONICAL_FIRST = "3 0 - 20 harmonics - Polygonal type"

#: The first control string this module wrote until 2026-10-30. It lacks the
#: ``- 20 harmonics -`` segment, so a direct load rejects it.
_SHORT_FIRST = "3 0 - Polygonal type"

#: The period in the second control string. The guide requires that string to
#: be present and non-empty and to carry the necessary information.
_PERIOD = re.compile(r"Period:\s*[-+0-9.eE]+")


@dataclass(frozen=True, slots=True)
class GgpFile:
    """A parsed ``.ggp``, with the control-string provenance kept.

    ``header_variant`` names the form the control strings are in, and
    ``format_valid`` says whether that form is one the guide documents as
    input: the documented pair, or none at all. It is False for control
    strings in some other format -- the two hashed variants, which
    ``np.savetxt(header=...)`` produces by accident, and the short type line
    this module wrote until 2026-10-30 -- which PCGrate reads only through the
    Border Profile Editor. Neither field is a test against PCGrate; both
    record what the guide says about the form the file is in. See the module
    docstring.
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
    """Parse a ``.ggp``, tolerating every control-string variant on disk.

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
    """Name the control-string variant and say whether the guide documents it.

    Exact matches come first, so a hashed copy of either type line falls
    through to ``commented-type-line`` rather than being read as the form it
    is a comment of. A ``canonical`` file whose second control string is
    missing is still canonical in form, but the guide rejects it: an absent or
    empty second string is the second error 10.2.3 lists.
    """
    if not header_lines:
        return "none", True
    if header_lines[0] == _CANONICAL_FIRST:
        return "canonical", _carries_a_period(header_lines)
    if header_lines[0] == _SHORT_FIRST:
        return "short-type-line", False
    if header_lines[0].lstrip("# ").startswith("3 0"):
        return "commented-type-line", False
    return "commented-columns", False


def _carries_a_period(header_lines: list[str]) -> bool:
    """Whether the second control string holds a period, as 10.2.2 requires."""
    return len(header_lines) > 1 and _PERIOD.search(header_lines[1]) is not None


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
    header: bool = False,
) -> Path:
    """Write a profile as ``x y`` points for PCGrate.

    Pass either a ``profile`` or explicit ``t``/``y`` arrays. A profile is
    reflected onto the file's ``+d-hat`` axis, the inverse of
    :func:`to_profile`. Explicit arrays are written as given, so they are
    file coordinates already, as the metrology pipeline's ``x_norm`` is.

    Values are written to six decimals, single-space separated, one point per
    line, with no line break after the last point. PCGrate rejects a file
    that ends in one, as ``np.savetxt`` output does. The format is pinned by
    ``tests/test_ggp.py``; PCGrate rejects variations on it.

    The points are written alone by default, with no control strings: Aditi
    reports that PCGrate's border-profile import does not want one
    (2026-10-29), and the guide confirms a file without control strings is
    legal input through the Border Profile Editor (10.2.3). ``header=True``
    starts the file with the two documented :data:`GGP_HEADER` lines instead,
    the first of which is byte-for-byte the one PCGrate's own exports carry.
    Reading a headerless file back reports ``header_variant == "none"``.
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
    points = (f"{fmt % x_value} {fmt % height}"
              for x_value, height in zip(t_values, y_values))
    with path.open("w", encoding="utf-8") as handle:
        # No line break after the last point: PCGrate reads a trailing one as
        # an extra, empty point and rejects the file (found by Aditi,
        # 2026-10-27).
        handle.write("\n".join([GGP_HEADER, *points] if header else points))
    return path
