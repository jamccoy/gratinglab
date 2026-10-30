"""The ``.txl`` writer, pinned against the reference layouts on disk.

Nothing here imports ``ebl_fracture``: this package does not depend on it, and
a writer validated only by the one reader it was written beside is not
validated. So the parser below is deliberately independent and deliberately
dumb -- it knows the grammar from ``docs`` and from the two files on disk, and
nothing about how :mod:`gratinglab.io.txl` chose to emit it.

The external pin is the 5 mm variable-line-spacing layout BEAMER exported
(``220412_ThorLabsCyl_5mmX5mm_GratingPattern.txl``). Its boundaries are not
reproduced byte for byte -- it runs its grooves in descending y and this
writer ascends -- but the invariant that makes it a grating is: every groove's
width is exactly half its own local period, while that period sweeps 0.954 to
1.050 µm across the pattern. Those numbers are quoted from the file.
"""

from __future__ import annotations

import math
import re

import pytest

from gratinglab.io.txl import DIALECTS, TxlError, dumps, groove_edges, write

# -- an independent reader ---------------------------------------------------

_PAIR = re.compile(r"(-?\d+(?:\.\d+)?)\s*,\s*(-?\d+(?:\.\d+)?)")


def parse(text):
    """``{structure: {'layer': int, 'boundaries': [[(x, y), ...]], 'arefs': [...]}}``.

    Coordinates come back in microns, as the file states them.
    """
    out, current, layer = {}, None, 0
    header = None
    for raw in text.splitlines():
        line = raw.split("!", 1)[0].strip()
        if not line:
            continue
        head = line.split()[0]
        if head in DIALECTS:
            header = (head, line.split()[1])
        elif head == "STRUCT":
            current = line.split()[1]
            out[current] = {"layer": layer, "boundaries": [], "arefs": []}
        elif head == "ENDSTRUCT":
            current = None
        elif head == "LAYER":
            layer = int(line.split()[1])
            if current is not None:
                out[current]["layer"] = layer
        elif head == "B":
            pts = [(float(a), float(b)) for a, b in _PAIR.findall(line)]
            assert line.endswith("ENDB"), f"boundary not closed: {line!r}"
            out[current]["boundaries"].append(pts)
        elif head == "AREF":
            counts = [int(t) for t in re.findall(r"(?<![,.\d-])(\d+)(?![,.\d])", line)]
            pts = [(float(a), float(b)) for a, b in _PAIR.findall(line)]
            out[current]["arefs"].append(
                {"name": line.split()[1], "origin": pts[0],
                 "n1": counts[0], "v1": pts[1], "n2": counts[1], "v2": pts[2]}
            )
    out["__header__"] = header
    return out


def rect(points):
    """``(x0, x1, y0, y1)`` of an axis-aligned four-point boundary."""
    xs, ys = [p[0] for p in points], [p[1] for p in points]
    assert len(points) == 4, points
    return min(xs), max(xs), min(ys), max(ys)


def across_spans(struct, grooves_along):
    """Each groove's extent across the dispersion direction, ascending."""
    boxes = [rect(b) for b in struct["boundaries"]]
    picked = [(y0, y1) for _, _, y0, y1 in boxes] if grooves_along == "x" else \
             [(x0, x1) for x0, x1, _, _ in boxes]
    return sorted(picked)


UNIFORM = dict(period_nm=1000.0, length_nm=1_000_000.0, across_nm=1_000_000.0)


# -- the envelope ------------------------------------------------------------


def test_the_envelope_is_a_complete_library():
    text = dumps(**UNIFORM, title="a 1 um grating")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    assert lines[0] == "!a 1 um grating"
    assert lines[1] == "TEXTLIB 9.0.0"
    assert lines[2] == "UNIT MICRON"
    assert lines[3] == "RESOLVE 0.00001"
    assert lines[4] == "BEGLIB"
    assert lines[-1] == "ENDLIB"
    assert text.endswith("\n")


@pytest.mark.parametrize("dialect", sorted(DIALECTS))
def test_both_dialects_declare_their_own_version(dialect):
    doc = parse(dumps(**UNIFORM, dialect=dialect))
    assert doc["__header__"] == (dialect, DIALECTS[dialect])


def test_an_unknown_dialect_is_refused():
    with pytest.raises(TxlError, match="dialect must be one of"):
        dumps(**UNIFORM, dialect="GDSII")


# -- geometry ----------------------------------------------------------------


def test_nanometres_in_microns_out():
    """§1 says the API is nanometres; UNIT MICRON says the file is not."""
    doc = parse(dumps(**UNIFORM, origin="corner"))
    spans = across_spans(doc["grating"], "x")
    assert spans[0] == pytest.approx((0.0, 0.5))          # 1000 nm -> 1.0 um
    assert spans[1][0] - spans[0][0] == pytest.approx(1.0)


def test_a_groove_is_the_duty_cycle_of_its_period():
    for duty in (0.2, 0.5, 0.8):
        doc = parse(dumps(**UNIFORM, duty_cycle=duty, origin="corner"))
        spans = across_spans(doc["grating"], "x")
        widths = [b - a for a, b in spans]
        assert widths == pytest.approx([duty] * len(widths), abs=1e-9)


def test_the_groove_count_is_what_fits():
    doc = parse(dumps(period_nm=1000.0, length_nm=5000.0, across_nm=10_400.0))
    assert len(doc["grating"]["boundaries"]) == 10


def test_centre_origin_centres_the_exposed_area():
    doc = parse(dumps(**UNIFORM))
    spans = across_spans(doc["grating"], "x")
    low, high = spans[0][0], spans[-1][1]
    # The pattern cell, not the exposed strip: the last groove's land is part
    # of the grating and is what makes the span symmetric.
    assert low == pytest.approx(-500.0)
    assert high + 0.5 == pytest.approx(500.0)


def test_corner_origin_starts_at_the_origin():
    doc = parse(dumps(**UNIFORM, origin="corner"))
    x0, x1, y0, _ = rect(doc["grating"]["boundaries"][0])
    assert (x0, y0) == pytest.approx((0.0, 0.0))
    assert x1 == pytest.approx(1000.0)


def test_grooves_along_y_swaps_the_axes():
    along_x = parse(dumps(**UNIFORM, grooves_along="x", origin="corner"))
    along_y = parse(dumps(**UNIFORM, grooves_along="y", origin="corner"))
    bx = rect(along_x["grating"]["boundaries"][0])
    by = rect(along_y["grating"]["boundaries"][0])
    assert bx == pytest.approx((0.0, 1000.0, 0.0, 0.5))
    assert by == pytest.approx((0.0, 0.5, 0.0, 1000.0))


def test_the_exposed_area_is_conserved():
    length_nm, across_nm, duty = 400_000.0, 200_000.0, 0.35
    doc = parse(dumps(period_nm=800.0, length_nm=length_nm,
                      across_nm=across_nm, duty_cycle=duty))
    boxes = [rect(b) for b in doc["grating"]["boundaries"]]
    area = sum((x1 - x0) * (y1 - y0) for x0, x1, y0, y1 in boxes)
    n = len(boxes)
    assert n == 250
    assert area == pytest.approx(n * (length_nm / 1000.0) * (800.0 * duty / 1000.0))


def test_the_layer_is_written_and_read_back():
    doc = parse(dumps(**UNIFORM, layer=6))
    assert doc["grating"]["layer"] == 6


# -- variable line spacing ---------------------------------------------------


def test_a_callable_period_varies_the_spacing_and_the_linewidth():
    """The reference layout's invariant: width is half the *local* period.

    ``220412_ThorLabsCyl_5mmX5mm_GratingPattern.txl`` sweeps 0.95403 µm at one
    edge to 1.05035 µm at the other across 5 mm, at exactly 50% duty
    throughout. This reproduces that sweep and checks the same invariant.
    """
    across_nm = 5_000_000.0
    lo, hi = 954.03, 1050.35

    def period(position_nm):
        return lo + (hi - lo) * position_nm / across_nm

    doc = parse(dumps(period_nm=period, length_nm=across_nm,
                      across_nm=across_nm, origin="corner"))
    spans = across_spans(doc["grating"], "x")
    assert len(spans) > 4900

    pitches = [b[0] - a[0] for a, b in zip(spans, spans[1:])]
    assert pitches[0] == pytest.approx(lo / 1000.0, abs=1e-4)
    assert pitches[-1] == pytest.approx(hi / 1000.0, abs=1e-3)
    # Non-decreasing to within one coordinate quantum: the pitches are
    # differences of coordinates already rounded to RESOLVE, so two grooves
    # whose true periods differ by less than that can come back level or a
    # dust-sized step apart in either direction.
    quantum = 2e-5
    assert all(b >= a - quantum for a, b in zip(pitches, pitches[1:]))
    assert pitches[-1] > pitches[0] + 0.09, "the sweep should be unmistakable"

    # Every groove is half its own local period, which is the invariant that
    # distinguishes a VLS grating from a uniform one with jitter.
    for (start, stop), pitch in zip(spans, pitches):
        assert (stop - start) == pytest.approx(pitch / 2.0, abs=2e-5)


def test_a_period_that_goes_non_positive_is_refused():
    with pytest.raises(TxlError, match="period must be positive"):
        dumps(period_nm=lambda pos: 1000.0 - pos, length_nm=1000.0,
              across_nm=2_000_000.0)


# -- the compact form --------------------------------------------------------


def test_compact_writes_one_groove_and_an_array():
    doc = parse(dumps(**UNIFORM, compact=True, origin="corner"))
    assert len(doc["groove"]["boundaries"]) == 1
    assert doc["grating"]["boundaries"] == []
    (aref,) = doc["grating"]["arefs"]
    assert aref["name"] == "groove"
    assert aref["n1"] == 1000
    assert aref["v1"] == pytest.approx((0.0, 1.0))
    assert aref["n2"] == 1


def test_compact_covers_the_same_ground_as_explicit():
    kw = dict(**UNIFORM, origin="corner")
    explicit = parse(dumps(**kw))
    compact = parse(dumps(**kw, compact=True))
    (aref,) = compact["grating"]["arefs"]
    assert aref["n1"] == len(explicit["grating"]["boundaries"])
    assert rect(compact["groove"]["boundaries"][0]) == pytest.approx(
        rect(explicit["grating"]["boundaries"][0]))


def test_compact_refuses_a_varying_period():
    with pytest.raises(TxlError, match="needs a constant period"):
        dumps(period_nm=lambda pos: 1000.0, length_nm=1000.0,
              across_nm=10_000.0, compact=True)


def test_compact_refuses_one_name_for_both_structures():
    with pytest.raises(TxlError, match="needs two structures"):
        dumps(**UNIFORM, compact=True, structure="g", groove_structure="g")


# -- the trailing-digit trap -------------------------------------------------


def test_a_structure_name_ending_in_a_digit_is_refused():
    """Verified against ebl_fracture's reader, 2026-10-30.

    Its ``_parse_ref`` recovers AREF counts as the bare integers on the line,
    and ``step1``'s trailing digit is one of them, so ``AREF step1 0,0 1000
    1,0 1 0,0`` yields counts ``['1', '1000', '1']`` and 1000 copies become 1.
    Every AREF in ``greyscale_1um_7step.txl`` parses with a different count for
    this reason. A file written here should not be the one that trips it.
    """
    with pytest.raises(TxlError, match="ends in a digit"):
        dumps(**UNIFORM, structure="grating1")
    with pytest.raises(TxlError, match="ends in a digit"):
        dumps(**UNIFORM, compact=True, groove_structure="step1")


def test_a_name_with_a_digit_inside_is_fine():
    doc = parse(dumps(**UNIFORM, structure="grating2b"))
    assert "grating2b" in doc


@pytest.mark.parametrize("name", ["", "two words", "   "])
def test_a_name_that_is_not_one_word_is_refused(name):
    with pytest.raises(TxlError, match="single word"):
        dumps(**UNIFORM, structure=name)


# -- resolution and the trivial refusals -------------------------------------


def test_coordinates_are_written_to_the_declared_resolution():
    text = dumps(**UNIFORM, resolve_um=0.001, origin="corner")
    assert "RESOLVE 0.001" in text
    first = [ln for ln in text.splitlines() if ln.startswith("B ")][0]
    for value in re.findall(r"-?\d+\.\d+", first):
        assert len(value.split(".")[1]) == 3


def test_nothing_fits_is_refused_rather_than_written_empty():
    with pytest.raises(TxlError, match="no groove fits"):
        dumps(period_nm=1000.0, length_nm=1000.0, across_nm=100.0)


@pytest.mark.parametrize("duty", [0.0, 1.0, -0.1, 1.5])
def test_an_impossible_duty_cycle_is_refused(duty):
    with pytest.raises(TxlError, match="duty_cycle must lie"):
        dumps(**UNIFORM, duty_cycle=duty)


@pytest.mark.parametrize("bad", [0.0, -1.0, math.inf, math.nan])
def test_an_impossible_length_is_refused(bad):
    with pytest.raises(TxlError, match="length_nm must be positive"):
        dumps(period_nm=1000.0, length_nm=bad, across_nm=10_000.0)


def test_groove_edges_is_usable_on_its_own():
    edges = groove_edges(1000.0, 3500.0, 0.4)
    assert edges == pytest.approx([(0.0, 400.0), (1000.0, 1400.0),
                                   (2000.0, 2400.0), (3000.0, 3400.0)])


def test_write_creates_parent_directories(tmp_path):
    out = write(tmp_path / "nested" / "grating.txl", **UNIFORM)
    assert out.is_file()
    assert parse(out.read_text(encoding="utf-8"))["grating"]["boundaries"]
