"""
Reading Bruker/Veeco Nanoscope files (.spm) directly.

A Nanoscope file is a latin-1 ASCII header of ``\\Key: value`` lines followed by
the binary image planes. One file holds several planes (typically Height Sensor
and Peak Force Error, each recorded in both scan directions), and each plane's
header section says where its bytes are and how to turn them into nanometres.

The reader was written here instead of taken from a library because the format
is small, this project depends on nothing beyond numpy/scipy/matplotlib, and
the conversion has two pitfalls (described below) that are best handled here.

The first is the height scale, which comes from a line like::

    \\@2:Z scale: V [Sens. ZsensSens] (0.0000000001872337 V/LSB) 0.804163 V

The name in brackets must be used verbatim. A file typically defines both::

    \\@Sens. Zsens:     V  32.46000 nm/V
    \\@Sens. ZsensSens: V 166.6319  nm/V

Using ``Zsens`` when the bracket says ``ZsensSens``, an easy mistake since one
looks like the other with a redundant suffix, yields heights 5.13x too small.
The result shows no visible sign of the error: the image is the right shape,
the grooves are in the right places, and every blaze angle is silently
incorrect. The reader was verified against a known export: with the right
parameter, the peak-to-peak height matches to 1.0000.

The second is the direction of the slow axis. The rows are returned reversed
against the order the plane stores them, because that is what Gwyddion does
and this project's measurements all came through Gwyddion exports. Gwyddion
reverses every Nanoscope import unconditionally -- ``modules/file/nanoscope.c``
calls ``gwy_data_field_invert(dfield, TRUE, FALSE, FALSE)`` straight after
scaling the plane, and that first argument is ``yflipped``, documented as "the
image will be flipped upside down". It does not consult ``\\Frame direction``,
which reads ``Down`` on the files here.

This showed up as a sign: the same grating measured -0.3722 deg of groove
rotation read as a .spm and +0.3721 deg read as the Gwyddion .txt export of
that scan. A flip does not move the blaze angle, so it went unnoticed until
``groove_orientation`` began correcting for the rotation it reported.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["SpmChannel", "read_spm", "list_channels", "is_nanoscope_file",
           "scan_y_size_um", "DEFAULT_CHANNEL", "DEFAULT_DIRECTION"]

DEFAULT_CHANNEL = "Height Sensor"

#: Retrace by default because it is the plane the project's existing Gwyddion
#: exports were taken from - so reading a .spm and reading its .txt export agree
#: by construction rather than by luck.
DEFAULT_DIRECTION = "Retrace"

#: The header is ASCII but the file is mostly binary; decode only far enough to
#: cover it. Real files put the first plane at ~80 kB, so this is generous.
_HEADER_BYTES = 262144

_MAGIC = b"\\*File list"

_SECTION_SPLIT = "\\*Ciao image list"

_BYTES_TO_DTYPE = {1: "<i1", 2: "<i2", 4: "<i4", 8: "<i8"}


@dataclass(frozen=True)
class SpmChannel:
    """One image plane in the file."""

    channel: str
    direction: str
    rows: int
    cols: int
    scan_x_size_um: float

    def describe(self) -> str:
        return (f"{self.channel} / {self.direction} "
                f"({self.rows}x{self.cols}, {self.scan_x_size_um:g} um)")


def is_nanoscope_file(path) -> bool:
    """
    True when this looks like a Nanoscope file.

    The check reads the content and ignores the extension. Nanoscope writes
    companion files with no extension at all (``sample_flatten.0_00003``), and
    those are the same format, so refusing them for want of a suffix would be
    arbitrary.
    """
    try:
        with open(path, "rb") as handle:
            return handle.read(len(_MAGIC)) == _MAGIC
    except OSError:
        return False


def _header_text(raw: bytes) -> str:
    return raw[:_HEADER_BYTES].decode("latin-1", errors="replace")


def _sensitivities(header: str) -> dict[str, float]:
    """Every ``\\@Sens. <name>: V <value> <unit>`` in the file."""
    return {name: float(value) for name, value in
            re.findall(r"\\@Sens\. (\w+): V\s+([\d.eE+-]+)", header)}


def _scan_size_um(section: str, header: str) -> float:
    """
    Scan width in microns.

    The image section carries ``Scan Size: 2.000 2.000 ~m`` (the ``~m`` is a
    mangled µm), and the scan list carries ``Scan Size: 2000 nm``. The section
    value is preferred, with the file-level value as the fallback.
    """
    for text in (section, header):
        match = re.search(r"\\Scan Size:\s*([\d.]+)(?:\s+[\d.]+)?\s*(~m|um|µm|nm)",
                          text)
        if match:
            value, unit = float(match.group(1)), match.group(2)
            return value / 1000.0 if unit == "nm" else value
    raise ValueError("no Scan Size found in the header")


#: The scan shape, ``\\Aspect Ratio: 4:1``. The colon is required: see
#: :func:`_aspect_ratio` for the other field of the same name.
_ASPECT_RE = re.compile(r"\\Aspect Ratio:\s*([\d.]+):([\d.]+)")


def _aspect_ratio(header: str) -> float:
    """
    Fast-axis extent divided by slow-axis extent. 1.0 when not stated.

    Nanoscope writes the *same* number twice in ``Scan Size`` even for a
    rectangular scan, and records the rectangle here instead: a 1.5 µm scan at
    ``4:1`` covers 1.5 µm along the fast axis and 0.375 µm along the slow one.
    Gwyddion reads it that way, and its export of such a file states
    ``# Width: 1.500 µm`` against ``# Height: 0.375 µm``.

    **Two different fields carry this name, and only one is the scan shape.**
    The scan list states it as a ratio with a colon, ``1:1`` or ``4:1``. Each
    plane section also carries an ``\\Aspect Ratio`` that is a bare number --
    0.352542 on the suite's own TASTE file -- which is not a scan shape and
    must not be read as one. Taking the section's value in preference, which
    is the obvious way to write this function given every other field here
    works that way, turns that file's 1.57 µm slow axis into 4.44 µm. So the
    colon is part of the pattern, and this reads the header rather than a
    section.

    The scan shape is independent of the sampling density. ``Samps/line`` and
    ``Lines`` are in the same ratio on most files, which is what makes their
    pixels square, but 25 scans in this project's archive disagree, and
    non-square pixels are an ordinary thing to record. The two are not
    interchangeable, and it is the aspect ratio that describes the area.

    A missing or unusable field gives 1.0, which is what this reader assumed
    unconditionally before.
    """
    match = _ASPECT_RE.search(header)
    if match:
        fast, slow = float(match.group(1)), float(match.group(2))
        if fast > 0 and slow > 0:
            return fast / slow
    return 1.0


def _slow_size_um(section: str, header: str, rows: int) -> float | None:
    """
    Physical length of the recorded lines along the slow axis, in microns.

    Two corrections to the stated ``Scan Size``, and both matter.

    The aspect ratio first. ``Scan Size`` carries two values, fast then slow
    (``2.000 2.000 ~m``), and on a rectangular scan they are *equal*, with the
    shape held in ``\\Aspect Ratio`` -- see :func:`_aspect_ratio`. Taking the
    second value at face value made the slow axis too long by that ratio on
    1013 of the 1245 scans in this project's archive, 911 of them at 4:1. Since
    a groove rotation is ``atan`` of a ratio of the two axes, every such scan
    reported a rotation too small by the same factor, and some declined a
    correction they needed because the understated angle fell below the
    threshold. The suite's own TASTE file is 1:1, which is why nothing caught
    it; `tests/metrology/test_spm_reader.py` now checks the arithmetic
    directly and, where the archive is reachable, against Gwyddion's own
    numbers.

    Then the lines actually recorded. The slow extent is what the scan was
    *set* to cover in ``\\Lines:`` lines, and a scan stopped early records
    fewer: the TASTE file holds 401 lines of a 512-line scan, so the recorded
    height is ``extent * rows / lines``.

    Returns None, rather than guessing, when the header carries only one
    ``Scan Size`` value or no ``\\Lines:`` count.
    """
    size = None
    for text in (section, header):
        match = re.search(r"\\Scan Size:\s*[\d.]+\s+([\d.]+)\s*(~m|um|µm|nm)", text)
        if match:
            value, unit = float(match.group(1)), match.group(2)
            size = value / 1000.0 if unit == "nm" else value
            break
    lines = re.search(r"\\Lines:\s*(\d+)", header)
    if size is None or lines is None or int(lines.group(1)) <= 0:
        return None
    extent = size / _aspect_ratio(header)
    return extent * rows / int(lines.group(1))


def scan_y_size_um(path, channel: str = DEFAULT_CHANNEL,
                   direction: str = DEFAULT_DIRECTION) -> float | None:
    """The slow-axis length of one plane's recorded lines, or None if unknown."""
    header = _header_text(Path(path).read_bytes())
    for section in _plane_sections(header):
        described = _describe_section(section, header)
        if (described is not None and described.channel == channel
                and described.direction == direction):
            return _slow_size_um(section, header, described.rows)
    return None


def _plane_sections(header: str) -> list[str]:
    return header.split(_SECTION_SPLIT)[1:]


def _describe_section(section: str, header: str) -> SpmChannel | None:
    channel = re.search(r'\\@2:Image Data: S \[\w+\] "([^"]+)"', section)
    rows = re.search(r"\\Number of lines: (\d+)", section)
    cols = re.search(r"\\Samps/line: (\d+)", section)
    if not (channel and rows and cols):
        return None
    direction = re.search(r"\\Line Direction: (\w+)", section)
    return SpmChannel(
        channel=channel.group(1),
        direction=direction.group(1) if direction else "Unknown",
        rows=int(rows.group(1)),
        cols=int(cols.group(1)),
        scan_x_size_um=_scan_size_um(section, header),
    )


def list_channels(path) -> tuple[SpmChannel, ...]:
    """Every image plane in the file, in the order the header lists them."""
    raw = Path(path).read_bytes()
    header = _header_text(raw)
    found = (_describe_section(s, header) for s in _plane_sections(header))
    return tuple(c for c in found if c is not None)


def read_spm(path, channel: str = DEFAULT_CHANNEL,
             direction: str = DEFAULT_DIRECTION):
    """
    Read one image plane.

    Parameters:
        channel: e.g. "Height Sensor". Only topography channels make sense for
            blaze-angle work; "Peak Force Error" is not a height map.
        direction: "Retrace" or "Trace".

    Returns:
        (data, scan_x_size_um), with data in metres, matching what the text
        loader returns, so nothing downstream needs to know which format it
        came from.
    """
    path = Path(path)
    raw = path.read_bytes()
    header = _header_text(raw)

    if not raw.startswith(_MAGIC):
        raise ValueError(f"{path.name} is not a Nanoscope file "
                         f"(expected it to start with '{_MAGIC.decode()}')")

    sensitivities = _sensitivities(header)

    for section in _plane_sections(header):
        described = _describe_section(section, header)
        if described is None:
            continue
        if described.channel != channel or described.direction != direction:
            continue

        offset = int(re.search(r"\\Data offset: (\d+)", section).group(1))
        bytes_per_pixel = int(re.search(r"\\Bytes/pixel: (\d+)", section).group(1))
        dtype = _BYTES_TO_DTYPE.get(bytes_per_pixel)
        if dtype is None:
            raise ValueError(f"unsupported Bytes/pixel: {bytes_per_pixel}")

        scale = re.search(
            r"\\@2:Z scale: V \[Sens\. (\w+)\] \(([\d.eE+-]+) V/LSB\)", section)
        if scale is None:
            raise ValueError(f"no '@2:Z scale' for {described.describe()}")

        # Verbatim. See the module docstring: 'Zsens' and 'ZsensSens' are both
        # present and differ by 5.13x.
        parameter, volts_per_lsb = scale.group(1), float(scale.group(2))
        if parameter not in sensitivities:
            raise ValueError(
                f"'@2:Z scale' refers to '\\@Sens. {parameter}', which this file "
                f"does not define. Known: {', '.join(sorted(sensitivities))}")
        nm_per_volt = sensitivities[parameter]

        count = described.rows * described.cols
        expected = count * bytes_per_pixel
        if offset + expected > len(raw):
            raise ValueError(
                f"{path.name} is truncated: {described.describe()} needs "
                f"{expected} bytes at offset {offset}, file is {len(raw)}")

        plane = np.frombuffer(raw, dtype=dtype, count=count, offset=offset)
        heights_nm = plane.reshape(described.rows, described.cols).astype(float) \
            * volts_per_lsb * nm_per_volt

        # Reversed: see the module docstring. Gwyddion's importer flips every
        # Nanoscope plane upside down, so its text export of a scan is
        # row-reversed against the plane as stored, and this reader exists to
        # agree with those exports -- the same reason DEFAULT_DIRECTION is
        # Retrace.
        return heights_nm[::-1] * 1e-9, described.scan_x_size_um

    available = list_channels(path)
    raise ValueError(
        f"{path.name} has no '{channel}' / '{direction}' plane. It contains: "
        + "; ".join(c.describe() for c in available))
