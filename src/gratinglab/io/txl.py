"""Writer for Raith TEXTLIB / LETXTLIB (``.txl``) grating layouts.

A ``.txl`` is a text layout: named structures holding polygons in microns,
plus array references that repeat a structure on a lattice. It is the format
BEAMER imports (``import_txl``), and the format ``ebl_fracture`` reads. This
module writes one for a grating, and writes nothing else: no fields, no
subfields, no dose, no instrument. That line is ``docs/roadmap.md``'s, under
"Why not absorb the fracture tool", and a text layout sits on this side of it
because it carries geometry only. What becomes of that geometry -- proximity
correction, fracture, dose assignment, the machine file -- belongs to the
layout tool.

The grammar written here, which is the intersection of the two dialects seen
on disk::

    !<title>                 optional, '!' starts a comment
    TEXTLIB 9.0.0            or LETXTLIB 1.0.0
    UNIT MICRON
    RESOLVE 0.0001
    BEGLIB
    STRUCT <name>
    LAYER <n>
    B x,y x,y x,y x,y ENDB
    ENDSTRUCT
    ENDLIB

Two shapes of output, because the reference files on disk come both ways:

*explicit* writes one boundary per groove, which is the only option when the
period varies across the pattern, and is what BEAMER's own export looks like
(thousands of boundaries, no references).

*compact* writes one groove and an ``AREF`` that repeats it, which is how a
uniform grating is written by hand -- two structures and a dozen lines instead
of a megabyte. Only available at constant period.

**Units.** Lengths in and out of this module's public functions are in
nanometres (``docs/conventions.md`` §1). The file is in microns, because
``UNIT MICRON`` is what both dialects on disk declare and what BEAMER expects;
the conversion happens here and nowhere else.

**Structure names must not end in a digit**, and :func:`dumps` refuses one
that does. ``ebl_fracture``'s reader recovers an ``AREF``'s repeat counts by
picking the bare integers out of the line, and a trailing digit in the name is
one of them: ``AREF step1 0,0 1000 1,0 1 0,0`` parses with 1000 copies becoming
1, because ``step1`` contributed the first count. It is the reader's bug, not
the format's, and a file this module writes should not be the one that trips
it. ``tests/test_txl.py`` pins the refusal and quotes the parse.
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Callable, Iterable, Sequence

__all__ = ["DIALECTS", "TxlError", "groove_edges", "dumps", "write"]

#: ``(header keyword, default version)`` for each dialect seen on disk.
DIALECTS = {"TEXTLIB": "9.0.0", "LETXTLIB": "1.0.0"}

_NM_PER_UM = 1000.0

#: A name ending in a digit breaks ebl_fracture's AREF count recovery; see the
#: module docstring.
_TRAILING_DIGIT = re.compile(r"\d$")


class TxlError(ValueError):
    """Raised when a layout cannot be written as valid TEXTLIB."""


def groove_edges(
    period_nm: float | Callable[[float], float],
    across_nm: float,
    duty_cycle: float = 0.5,
) -> list[tuple[float, float]]:
    """Where each groove starts and stops, across the dispersion direction.

    Returns ``(start_nm, stop_nm)`` pairs measured from 0, ascending, for the
    grooves that fit inside ``across_nm``. The exposed width of each is
    ``duty_cycle`` of its *local* period, so a varying period varies the
    linewidth with it -- which is what the 5 mm reference layout does, at
    exactly half the local period from 0.954 µm at one edge to 1.050 µm at the
    other.

    ``period_nm`` is either a constant or a callable taking the position of a
    groove's leading edge, in nm from 0, and returning the local period there.
    A callable is stepped rather than inverted: the next groove starts one
    local period beyond this one.
    """
    if not 0.0 < duty_cycle < 1.0:
        raise TxlError(f"duty_cycle must lie in (0, 1), got {duty_cycle}")
    if not (across_nm > 0.0 and math.isfinite(across_nm)):
        raise TxlError(f"across_nm must be positive and finite, got {across_nm}")

    local = period_nm if callable(period_nm) else (lambda _pos, p=float(period_nm): p)

    edges: list[tuple[float, float]] = []
    position = 0.0
    while True:
        period = float(local(position))
        if not (period > 0.0 and math.isfinite(period)):
            raise TxlError(
                f"period must be positive and finite, got {period} at "
                f"{position} nm across the pattern"
            )
        stop = position + period * duty_cycle
        if stop > across_nm:
            break
        edges.append((position, stop))
        position += period
        if position >= across_nm:
            break
    if not edges:
        raise TxlError(
            f"no groove fits: {across_nm} nm across holds none of a "
            f"{float(local(0.0))} nm period at duty cycle {duty_cycle}"
        )
    return edges


def _check_name(name: str) -> str:
    if not name or name.split() != [name]:
        raise TxlError(f"structure name must be a single word, got {name!r}")
    if _TRAILING_DIGIT.search(name):
        raise TxlError(
            f"structure name {name!r} ends in a digit. ebl_fracture's reader "
            "recovers an AREF's repeat counts as the bare integers on the "
            "line, so a trailing digit is read as a count: 'AREF "
            f"{name} 0,0 1000 1,0 1 0,0' parses with 1000 copies becoming "
            f"{name[-1]}. Use a name that does not end in a digit."
        )
    return name


def _rect(
    start: float, stop: float, length_nm: float, grooves_along: str,
    offset: tuple[float, float],
) -> list[tuple[float, float]]:
    """One groove as four corners in nm, counter-clockwise from its origin."""
    dx, dy = offset
    if grooves_along == "x":
        return [(dx, start + dy), (dx + length_nm, start + dy),
                (dx + length_nm, stop + dy), (dx, stop + dy)]
    return [(start + dx, dy), (stop + dx, dy),
            (stop + dx, dy + length_nm), (start + dx, dy + length_nm)]


def _format_points(points: Iterable[Sequence[float]], decimals: int) -> str:
    return " ".join(
        f"{x / _NM_PER_UM:.{decimals}f},{y / _NM_PER_UM:.{decimals}f}"
        for x, y in points
    )


def dumps(
    *,
    period_nm: float | Callable[[float], float],
    length_nm: float,
    across_nm: float,
    duty_cycle: float = 0.5,
    grooves_along: str = "x",
    origin: str = "centre",
    layer: int = 0,
    structure: str = "grating",
    groove_structure: str = "groove",
    title: str | None = None,
    compact: bool = False,
    resolve_um: float = 1e-5,
    dialect: str = "TEXTLIB",
) -> str:
    """A grating layout as TEXTLIB source.

    ``length_nm`` runs along the grooves and ``across_nm`` across them, so
    ``across_nm`` is the dispersion direction and the one the period divides.
    ``grooves_along`` says which axis the grooves themselves run along,
    ``"x"`` or ``"y"``; both appear in the reference files.

    ``origin`` is ``"centre"``, which puts the pattern's centre at ``(0, 0)``
    as BEAMER's export does, or ``"corner"``, which puts its low corner there
    as a hand-written file does.

    ``compact`` writes one groove plus an ``AREF`` instead of a boundary per
    groove. It needs a constant period, since an array repeats one shape on a
    regular lattice and a varying period is not one.

    Coordinates are written to the number of decimals ``resolve_um`` implies,
    since a coordinate finer than the declared grid is not meaningful.
    """
    if dialect not in DIALECTS:
        raise TxlError(f"dialect must be one of {sorted(DIALECTS)}, got {dialect!r}")
    if grooves_along not in ("x", "y"):
        raise TxlError(f"grooves_along must be 'x' or 'y', got {grooves_along!r}")
    if origin not in ("centre", "corner"):
        raise TxlError(f"origin must be 'centre' or 'corner', got {origin!r}")
    if not (length_nm > 0.0 and math.isfinite(length_nm)):
        raise TxlError(f"length_nm must be positive and finite, got {length_nm}")
    if not (resolve_um > 0.0 and math.isfinite(resolve_um)):
        raise TxlError(f"resolve_um must be positive and finite, got {resolve_um}")
    if compact and callable(period_nm):
        raise TxlError(
            "compact output repeats one groove on a regular lattice, so it "
            "needs a constant period; pass compact=False to write a boundary "
            "per groove instead"
        )
    _check_name(structure)
    if compact:
        _check_name(groove_structure)
        if groove_structure == structure:
            raise TxlError(
                "compact output needs two structures, the groove and the array "
                f"that repeats it, but both are named {structure!r}"
            )

    decimals = max(0, -int(math.floor(math.log10(resolve_um))))
    edges = groove_edges(period_nm, across_nm, duty_cycle)

    # The pattern spans the grooves actually written, not the requested extent:
    # a period that does not divide across_nm leaves a remainder unexposed, and
    # centring on the request would put the exposed area off centre. The last
    # groove's own period is recovered from its width, which is the duty cycle
    # of it, so this holds for a varying period too.
    last_start, last_stop = edges[-1]
    span = last_start + (last_stop - last_start) / duty_cycle
    if origin == "centre":
        half_across, half_length = span / 2.0, length_nm / 2.0
        offset = (-half_length, -half_across) if grooves_along == "x" else (-half_across, -half_length)
    else:
        offset = (0.0, 0.0)

    lines: list[str] = []
    if title is not None:
        lines.append(f"!{title}")
    lines += [
        f"{dialect} {DIALECTS[dialect]}",
        "UNIT MICRON",
        f"RESOLVE {resolve_um:.{decimals}f}",
        "BEGLIB",
        "",
    ]

    if compact:
        pitch = float(period_nm)
        lines += [f"STRUCT {groove_structure}", f"LAYER {int(layer)}"]
        corners = _rect(edges[0][0], edges[0][1], length_nm, grooves_along, offset)
        lines += [f"B {_format_points(corners, decimals)} ENDB", "ENDSTRUCT", ""]
        step = (0.0, pitch) if grooves_along == "x" else (pitch, 0.0)
        lines += [
            f"STRUCT {structure}",
            f"AREF {groove_structure} 0,0 {len(edges)} "
            f"{_format_points([step], decimals)} 1 0,0",
            "ENDSTRUCT",
            "",
        ]
    else:
        lines += [f"STRUCT {structure}", f"LAYER {int(layer)}"]
        for start, stop in edges:
            corners = _rect(start, stop, length_nm, grooves_along, offset)
            lines.append(f"B {_format_points(corners, decimals)} ENDB")
        lines += ["ENDSTRUCT", ""]

    lines.append("ENDLIB")
    return "\n".join(lines) + "\n"


def write(path: str | Path, **kwargs) -> Path:
    """Write :func:`dumps` to ``path``, creating parent directories."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dumps(**kwargs), encoding="utf-8")
    return path
