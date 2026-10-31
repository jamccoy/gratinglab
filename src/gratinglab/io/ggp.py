"""Reader and writer for PCGrate ``.ggp`` boundary profiles.

The reader takes all five of PCGrate's border profile types; the writer writes
the polygonal one, which is the only type a measured groove can be expressed
in. See "The other border profile types" below. Until 2026-10-31 the reader
took the polygonal type alone and refused 18 of the 53 files PCGrate ships.

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

The other border profile types
------------------------------

Section 10 (p. 59) names five major types, and the first number of the first
control string selects one: ``0`` Trapezoidal, ``1`` Sine trapezoidal, ``2``
Trigonometric, ``3`` Polygonal, and Section-dependent, which carries no code
in any file here and is therefore refused by code rather than guessed at. The
second number selects a subtype -- General, Triangular, Lamellar and Sawtooth
for Trapezoidal, General and Sinusoidal for Sine trapezoidal (10.1, p. 59).
:data:`PROFILE_TYPES` is the table.

Every type is returned as ``t``/``y``, so a caller does not branch on the type
to get geometry, with the analytic parameters kept alongside.

The guide documents a *file* layout for two of the types only: polygonal in
10.2.2 and trigonometric in 10.3.2 (p. 66), which gives ``C`` or ``S``, a
harmonic number and a Fourier amplitude per line, so the border is
``y(x) = sum C_n cos(2 pi n x) + sum S_n sin(2 pi n x)``. The simple-type files
in ``BorderProfiles/`` are Border Profile Editor exports and their layout is
documented nowhere, so it is read off the files: five numbers, one a line,
``left angle, right angle, ridge, land, depth`` -- the names 10.1 and the
glossary give, with ridge, land and depth relative to the period (10, p. 61,
"values of Gk/d, Bk/d, and hk/d are constant").

Inference is what the 2026-10-30 header correction was about, so that order is
not taken on faith. It *closes*: the ridge, the land and the two slopes'
horizontal extents span exactly one period, and they do for all nine
trapezoidal files to within their six-figure rounding. The sine family is the
same argument with "frequencies F1k and F2k" (10, p. 61) in cycles per period,
so a flank is the half cycle spanning ``1 / (2F)``; that closes for all seven.
Sixteen of sixteen, and :func:`_trapezoid` and :func:`_sine_trapezoid` refuse a
file that does not close, so a wrong order cannot pass unnoticed.

The trigonometric convention is checked against PCGrate rather than assumed:
the installation ships ``Lamellar_converted-to-trigonometric.ggp``, which is
the program's own conversion of ``Lamellar.ggp``. Its odd sine harmonics are
``0.5 / (n pi)`` to one part in 1e9 with ``C 0`` at 0.125 -- the series of a
square wave of amplitude 0.125 about 0.125, which is ``Lamellar.ggp``'s depth
of 0.25 at its 50% duty. The depth read back is 8.951% per side larger, the
Gibbs overshoot of a series truncated at ``n = 100`` against the
Wilbraham-Gibbs 8.9490%; that is a property of the vendor's file and not of
this reader.

The period, and what the file does not carry
--------------------------------------------

The second control string carries a period, and all 53 of PCGrate's files have
one: 34 give 1, and 19 give a real number -- 163, 165.755, 160, 161.3, 176.367,
1000, 11494.3, 2. :attr:`GgpFile.period` is now that number, which this module
discarded until 2026-10-31.

What the file does not carry is the *unit*: a period of 163 is 163 of whatever
the task was set in. So the sidecar is still needed, and a headerless file --
which is this writer's own default, and what Aditi's import wants -- carries no
period at all. Anything downstream that needs a period in nm must still obtain
it from another source, and ``benchmarks/corpus.toml`` is that source, recovered
by hand. Code holding a measured profile should reach for
:meth:`gratinglab.metrology.boundary.BoundaryProfile.to_problem` instead.

The profile scaling coefficient resizes the profile along ``y`` (10.3.1), and
:func:`read_ggp` applies it. No file PCGrate ships exercises it: 21 omit it and
32 give 1, so the path is documented but unvalidated against the program, and
it is a no-op on every real file.
See ``docs/roadmap.md`` on the native boundary format for the file-based
version of the same fix.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from numpy.typing import ArrayLike

from ..profiles import FromProfileData

__all__ = [
    "GgpFile",
    "GGP_HEADER",
    "MAX_HARMONIC",
    "MAX_HARMONIC_AMPLITUDE",
    "PROFILE_TYPES",
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

#: The period's value, and the profile scaling coefficient's. Every one of
#: PCGrate's 53 reference profiles carries a period; 19 carry one other than 1
#: (163, 165.755, 160, 161.3, 176.367, 1000, 11494.3, 2), which is a number
#: this module did not read until 2026-10-31. The unit is still absent, so a
#: period of 163 is 163 of whatever the task is in -- see the module docstring.
_PERIOD_VALUE = re.compile(r"Period:\s*([-+0-9.eE]+)")
_PSC_VALUE = re.compile(r"PSC:\s*([-+0-9.eE]+)")

#: A first control string's leading type code, with an optional ``#`` so a
#: hashed copy is still recognised as the type it is a comment of.
_TYPE_CODE = re.compile(r"^#?\s*(\d+)\s+(\d+)\s*-")

#: The border profile types, by the code the first control string carries.
#: Section 10 (p. 59) names five major types -- Trapezoidal, Sine trapezoidal,
#: Polygonal, Trigonometric and Section-dependent -- and says the first two
#: have subtypes, which the second number selects: General, Triangular,
#: Lamellar and Sawtooth for Trapezoidal, General and Sinusoidal for Sine
#: trapezoidal (10.1, p. 59). Section-dependent carries no code of its own in
#: any file here, so it is deliberately absent rather than guessed at.
PROFILE_TYPES: dict[tuple[int, int], str] = {
    (0, 0): "trapezoidal",
    (0, 1): "lamellar",
    (0, 2): "triangular",
    (0, 3): "sawtooth",
    (1, 0): "sine-trapezoidal",
    (1, 1): "sinusoidal",
    (2, 0): "trigonometric",
    (3, 0): "polygonal",
}

#: Codes whose five numbers are ``left angle, right angle, ridge, land,
#: depth``, and codes whose first two are sine frequencies instead.
_TRAPEZOIDAL_CODES = frozenset({(0, 0), (0, 1), (0, 2), (0, 3)})
_SINE_CODES = frozenset({(1, 0), (1, 1)})

#: 10.3.3's bounds on a trigonometric file: the harmonic number is an integer
#: in 0 .. 5000 and unique within its type, the amplitude a real in -50 .. 50,
#: and the file holds 2 .. 10 001 data strings.
MAX_HARMONIC = 5000
MAX_HARMONIC_AMPLITUDE = 50.0
MAX_HARMONIC_STRINGS = 10_001

#: How many points the analytic types are sampled onto by default. 2000 is
#: `ggp_n_points`'s default in :mod:`gratinglab.metrology.settings` and sits
#: inside 10.2.3's 5000-data-string cap, so a profile read from an analytic
#: file and written back out as polygonal is still a file PCGrate will load.
DEFAULT_SAMPLES = 2000


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
    #: The type the first control string names, from :data:`PROFILE_TYPES`.
    #: ``"polygonal"`` for a file with no recognisable type code, which is
    #: what a headerless points-only file is.
    profile_type: str = "polygonal"
    #: The two numbers the first control string leads with, or ``None`` when
    #: there is no control string to read them from.
    type_code: tuple[int, int] | None = None
    #: The named parameters of an analytic type, empty for the polygonal and
    #: trigonometric types, which carry geometry rather than parameters.
    parameters: Mapping[str, float] = field(default_factory=dict)
    #: A trigonometric file's Fourier coefficients, keyed ``"C"`` (cosine,
    #: the even harmonics) and ``"S"`` (sine, the odd), each mapping harmonic
    #: number to amplitude. Empty for every other type.
    harmonics: Mapping[str, Mapping[int, float]] = field(default_factory=dict)
    #: The period the second control string carries, or ``None`` when there is
    #: no second control string. A bare number: the file does not record its
    #: unit. See the module docstring.
    period: float | None = None
    #: The profile scaling coefficient, 1 when the control string omits it.
    psc: float = 1.0

    @property
    def depth(self) -> float:
        """Peak-to-valley height, in units of the period."""
        return float(np.ptp(self.y))


def read_ggp(path: str | Path, *, samples: int = DEFAULT_SAMPLES) -> GgpFile:
    """Parse a ``.ggp`` of any of PCGrate's border profile types.

    The first control string's type code selects the parse
    (:data:`PROFILE_TYPES`); a file with no recognisable code is read as
    points, which is what a headerless polygonal file is. Every type arrives
    as ``t``/``y`` so a caller need not branch on the type to get geometry,
    with the analytic types' named parameters kept alongside in
    :attr:`GgpFile.parameters` and a trigonometric file's Fourier
    coefficients in :attr:`GgpFile.harmonics`.

    The trapezoidal family is returned as its exact vertices rather than as
    ``samples`` points, because it *is* a polygon and sampling it would only
    add error. ``samples`` therefore applies to the sine trapezoidal and
    trigonometric types, whose flanks are curved.

    Raises
    ------
    ValueError
        On a malformed data line, or on a type code outside
        :data:`PROFILE_TYPES`. The error is deliberate, because a truncated
        point that was silently dropped would shift the profile and change
        every efficiency computed from it.
    """
    path = Path(path)
    lines = _read_text(path).splitlines()
    first = next((line.strip() for line in lines if line.strip()), "")

    code = _type_code(first)
    if code is not None and code not in PROFILE_TYPES:
        raise ValueError(
            f"{path}: border profile type code {code[0]} {code[1]} is not one "
            f"this module knows. Known codes: "
            + ", ".join(f"{a} {b} ({name})" for (a, b), name in PROFILE_TYPES.items())
            + f". The first control string was {first!r}."
        )

    name = PROFILE_TYPES.get(code, "polygonal") if code is not None else "polygonal"
    if code in _TRAPEZOIDAL_CODES or code in _SINE_CODES:
        return _read_simple(path, lines, code, name)
    if code == (2, 0):
        return _read_trigonometric(path, lines, code, name, samples)
    return _read_polygonal(path, lines, code, name)


def _type_code(first_line: str) -> tuple[int, int] | None:
    """The type code a first control string leads with, or ``None``.

    ``None`` for a line that is not a control string at all, which is how a
    headerless points-only file reaches the polygonal parse. A data line
    cannot be mistaken for a code: the pattern needs the ``-`` that separates
    the code from the type's name, and a point's two numbers are not followed
    by one.
    """
    match = _TYPE_CODE.match(first_line)
    return (int(match.group(1)), int(match.group(2))) if match else None


def _control_strings(lines: list[str], *, keep: int) -> list[str]:
    """The first ``keep`` non-blank lines, stripped."""
    found = [line.strip() for line in lines if line.strip()]
    return found[:keep]


def _period_and_psc(control: list[str]) -> tuple[float | None, float]:
    """The period and scaling coefficient the second control string carries.

    PCGrate's own files put both on one line, ``Period: 1 PSC: 1``, and 21 of
    the 53 omit the PSC entirely, which 10.3.1 says means 1. Searching every
    control string rather than only the second keeps a file that splits them
    across lines readable.
    """
    period: float | None = None
    psc = 1.0
    for line in control:
        if period is None and (match := _PERIOD_VALUE.search(line)):
            period = float(match.group(1))
        if match := _PSC_VALUE.search(line):
            psc = float(match.group(1))
    return period, psc


def _read_polygonal(
    path: Path, lines: list[str], code: tuple[int, int] | None, name: str
) -> GgpFile:
    """Parse the polygonal type: two numbers a line, one point each."""
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
    period, psc = _period_and_psc(header_lines)

    return GgpFile(
        t=array[:, 0],
        y=array[:, 1] * psc,
        header_variant=variant,
        format_valid=compatible,
        source=path,
        profile_type=name,
        type_code=code,
        period=period,
        psc=psc,
    )


def _simple_values(path: Path, lines: list[str]) -> list[float]:
    """The five numbers an analytic simple type carries, one a line.

    The guide documents a file layout for the polygonal (10.2.2) and
    trigonometric (10.3.2) types only; the trapezoidal and sine trapezoidal
    files in PCGrate's ``BorderProfiles/`` are Border Profile Editor exports,
    so this layout is read off those files rather than out of the manual. What
    makes it safe to rely on is that it closes: see the module docstring.
    """
    values: list[float] = []
    for number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or _TYPE_CODE.match(stripped) or _PERIOD_VALUE.search(stripped):
            continue
        fields = stripped.split()
        if len(fields) != 1:
            raise ValueError(
                f"{path}:{number}: expected a single number on a parameter "
                f"line, got {stripped!r}"
            )
        try:
            values.append(float(fields[0]))
        except ValueError:
            raise ValueError(
                f"{path}:{number}: {fields[0]!r} is not a number"
            ) from None

    if len(values) != 5:
        raise ValueError(
            f"{path}: expected 5 parameter lines (left angle or frequency, "
            f"right angle or frequency, ridge, land, depth), got {len(values)}"
        )
    return values


def _read_simple(
    path: Path, lines: list[str], code: tuple[int, int], name: str
) -> GgpFile:
    """Parse a trapezoidal- or sine-trapezoidal-family file into vertices."""
    left, right, ridge, land, depth = _simple_values(path, lines)
    period, psc = _period_and_psc(_control_strings(lines, keep=2))

    if depth <= 0:
        raise ValueError(f"{path}: depth must be positive, got {depth}")
    if ridge < 0 or land < 0:
        raise ValueError(
            f"{path}: ridge and land cannot be negative, got {ridge} and {land}"
        )

    if code in _TRAPEZOIDAL_CODES:
        keys = ("left_angle_deg", "right_angle_deg", "ridge", "land", "depth")
        t_values, y_values = _trapezoid(path, left, right, ridge, land, depth)
    else:
        keys = ("left_frequency", "right_frequency", "ridge", "land", "depth")
        t_values, y_values = _sine_trapezoid(path, left, right, ridge, land, depth)

    return GgpFile(
        t=t_values,
        y=y_values * psc,
        header_variant="canonical",
        format_valid=True,
        source=path,
        profile_type=name,
        type_code=code,
        parameters=dict(zip(keys, (left, right, ridge, land, depth))),
        period=period,
        psc=psc,
    )


def _runs(path: Path, left: float, right: float, depth: float) -> tuple[float, float]:
    """The horizontal extent of the two slopes, from their angles.

    The angles are measured from the plane of the grating, so a slope of angle
    ``A`` and height ``depth`` spans ``depth / tan A`` horizontally. Both
    bounds are open: a 0 degree slope has no height and a 90 degree one has no
    width, and the lamellar files approach the latter with 89.9999 rather than
    90, which is itself evidence the program takes the tangent too.
    """
    for label, angle in (("left", left), ("right", right)):
        if not 0 < angle < 90:
            raise ValueError(
                f"{path}: {label} angle must be in (0, 90) degrees, got {angle}"
            )
    return (
        depth / math.tan(math.radians(left)),
        depth / math.tan(math.radians(right)),
    )


def _close_to_a_period(span: float) -> bool:
    """Whether the four horizontal extents span one period, within rounding.

    PCGrate's files carry the five parameters to six or so figures, so the sum
    lands within about 1e-5 of 1 rather than on it.
    """
    return abs(span - 1.0) <= 1e-4


def _trapezoid(
    path: Path, left: float, right: float, ridge: float, land: float, depth: float
) -> tuple[np.ndarray, np.ndarray]:
    """One period of a trapezoid, as its vertices.

    Laid out land first so the profile starts and ends at the valley, which
    makes the first point ``(0, 0)`` and the last ``(1, 0)`` -- the polygonal
    type's own convention, so the result can be written straight back out with
    :func:`write_ggp`. Consecutive duplicates are dropped, which is what turns
    a triangular or sawtooth file (ridge and land both zero) into the three
    vertices it actually has.
    """
    run_left, run_right = _runs(path, left, right, depth)
    span = land + run_left + ridge + run_right
    if not _close_to_a_period(span):
        raise ValueError(
            f"{path}: the profile spans {span:.6f} of a period, not 1. "
            f"ridge {ridge} + land {land} + {depth}/tan({left} deg) + "
            f"{depth}/tan({right} deg) should close, and that it does for "
            "every file PCGrate ships is what fixes the order of the five "
            "numbers. A span other than 1 means they are in another order."
        )
    edges = [0.0, land, land + run_left, land + run_left + ridge, 1.0]
    heights = [0.0, 0.0, depth, depth, 0.0]
    return _drop_repeats(np.asarray(edges), np.asarray(heights))


def _sine_trapezoid(
    path: Path, left: float, right: float, ridge: float, land: float, depth: float
) -> tuple[np.ndarray, np.ndarray]:
    """One period of a sine trapezoid: flat land, sine flank, flat ridge, sine flank.

    ``left`` and ``right`` are the flanks' sine *frequencies* in cycles per
    period (10, p. 61, "frequencies F1k and F2k"), so a flank is the half
    cycle that climbs from valley to ridge and spans ``1 / (2 F)``. That
    reading is what closes: see the module docstring.
    """
    for label, frequency in (("left", left), ("right", right)):
        if frequency <= 0:
            raise ValueError(
                f"{path}: {label} sine frequency must be positive, got {frequency}"
            )
    width_left, width_right = 1.0 / (2.0 * left), 1.0 / (2.0 * right)
    span = land + width_left + ridge + width_right
    if not _close_to_a_period(span):
        raise ValueError(
            f"{path}: the profile spans {span:.6f} of a period, not 1. "
            f"ridge {ridge} + land {land} + 1/(2*{left}) + 1/(2*{right}) "
            "should close, as it does for every sine trapezoidal file "
            "PCGrate ships."
        )

    # Split the sample budget between the two curved flanks in proportion to
    # their width, so a narrow flank is not resolved more finely than a wide
    # one. The flats need two points each and no more.
    per_flank = max(2, DEFAULT_SAMPLES // 2)
    rise = np.linspace(0.0, width_left, per_flank)
    fall = np.linspace(0.0, width_right, per_flank)
    t_values = np.concatenate([
        [0.0],
        land + rise,
        [land + width_left + ridge],
        land + width_left + ridge + fall,
    ])
    heights = np.concatenate([
        [0.0],
        depth / 2.0 * (1.0 - np.cos(2.0 * np.pi * left * rise)),
        [depth],
        depth / 2.0 * (1.0 + np.cos(2.0 * np.pi * right * fall)),
    ])
    t_values[-1] = 1.0
    return _drop_repeats(t_values, heights)


def _read_trigonometric(
    path: Path, lines: list[str], code: tuple[int, int], name: str, samples: int
) -> GgpFile:
    """Parse the trigonometric type and sample its Fourier series.

    The file holds ``C`` or ``S``, a harmonic number and an amplitude per line
    (10.3.2), so the border is
    ``y(x) = sum C_n cos(2 pi n x) + sum S_n sin(2 pi n x)`` with ``x``
    normalised to the period. ``C 0`` is a constant, which is exactly the
    guide's note that "the zeroth even (cosine) harmonic is equal to the shift
    of the zeroth level of the border function" (10.3.3).
    """
    if samples < 2:
        raise ValueError(f"samples must be at least 2, got {samples}")

    control = _control_strings(lines, keep=2)
    period, psc = _period_and_psc(control)
    coefficients: dict[str, dict[int, float]] = {"C": {}, "S": {}}
    strings = 0

    for number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped in control:
            continue
        fields = stripped.split()
        if len(fields) != 3 or fields[0] not in ("C", "S"):
            raise ValueError(
                f"{path}:{number}: expected 'C' or 'S', a harmonic number and "
                f"an amplitude, got {stripped!r}"
            )
        kind = fields[0]
        try:
            order, amplitude = int(fields[1]), float(fields[2])
        except ValueError:
            raise ValueError(
                f"{path}:{number}: {stripped!r} does not carry an integer "
                "harmonic number and a real amplitude"
            ) from None
        if not 0 <= order <= MAX_HARMONIC:
            raise ValueError(
                f"{path}:{number}: harmonic number {order} is outside "
                f"0 .. {MAX_HARMONIC}"
            )
        if abs(amplitude) > MAX_HARMONIC_AMPLITUDE:
            raise ValueError(
                f"{path}:{number}: amplitude {amplitude} is outside "
                f"-{MAX_HARMONIC_AMPLITUDE:g} .. {MAX_HARMONIC_AMPLITUDE:g}"
            )
        if order in coefficients[kind]:
            raise ValueError(
                f"{path}:{number}: harmonic {kind} {order} appears twice, and "
                "10.3.3 requires it to be unique within its type"
            )
        coefficients[kind][order] = amplitude
        strings += 1

    if not 2 <= strings <= MAX_HARMONIC_STRINGS:
        raise ValueError(
            f"{path}: {strings} data string(s); 10.3.3 allows "
            f"2 .. {MAX_HARMONIC_STRINGS}"
        )

    t_values = np.linspace(0.0, 1.0, samples)
    heights = np.zeros_like(t_values)
    for order, amplitude in coefficients["C"].items():
        heights += amplitude * np.cos(2.0 * np.pi * order * t_values)
    for order, amplitude in coefficients["S"].items():
        heights += amplitude * np.sin(2.0 * np.pi * order * t_values)

    return GgpFile(
        t=t_values,
        y=heights * psc,
        header_variant="canonical",
        format_valid=_carries_a_period(control),
        source=path,
        profile_type=name,
        type_code=code,
        harmonics={kind: dict(sorted(v.items())) for kind, v in coefficients.items()},
        period=period,
        psc=psc,
    )


def _drop_repeats(
    t_values: np.ndarray, heights: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Drop points identical to their predecessor, keeping the first."""
    keep = np.ones(len(t_values), dtype=bool)
    keep[1:] = (np.diff(t_values) != 0) | (np.diff(heights) != 0)
    return t_values[keep], heights[keep]


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
