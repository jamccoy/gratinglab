"""
Reading Bruker/Veeco Nanoscope files (.spm) directly.

A Nanoscope file is a latin-1 ASCII header of ``\\Key: value`` lines followed by
the binary image planes. One file holds several planes - typically Height Sensor
and Peak Force Error, each recorded in both scan directions - and each plane's
header section says where its bytes are and how to turn them into nanometres.

Written rather than taken from a library because the format is small, this
project deliberately depends on nothing beyond numpy/scipy/matplotlib, and the
conversion has one trap (below) that is worth owning explicitly.

**The trap.** The height scale comes from a line like::

    \\@2:Z scale: V [Sens. ZsensSens] (0.0000000001872337 V/LSB) 0.804163 V

The name in brackets must be used **verbatim**. A file typically defines both::

    \\@Sens. Zsens:     V  32.46000 nm/V
    \\@Sens. ZsensSens: V 166.6319  nm/V

Reaching for ``Zsens`` when the bracket says ``ZsensSens`` - which is easy, since
one looks like the other with a redundant suffix - yields heights 5.13x too small.
Nothing about the result looks wrong: the image is the right shape, the grooves
are in the right places, and every blaze angle is quietly incorrect. Verified
against a known export: with the right parameter, peak-to-peak matches to 1.0000.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["SpmChannel", "read_spm", "list_channels", "is_nanoscope_file",
           "DEFAULT_CHANNEL", "DEFAULT_DIRECTION"]

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

    Checked by content, not by extension. Nanoscope writes companion files with
    no extension at all (``sample_flatten.0_00003``), and those are the same
    format - refusing them for want of a suffix would be arbitrary.
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

    The image section carries ``Scan Size: 2.000 2.000 ~m`` - the ``~m`` is a
    mangled µm - and the scan list carries ``Scan Size: 2000 nm``. Prefer the
    section, fall back to the file-level value.
    """
    for text in (section, header):
        match = re.search(r"\\Scan Size:\s*([\d.]+)(?:\s+[\d.]+)?\s*(~m|um|µm|nm)",
                          text)
        if match:
            value, unit = float(match.group(1)), match.group(2)
            return value / 1000.0 if unit == "nm" else value
    raise ValueError("no Scan Size found in the header")


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
        (data, scan_x_size_um) - data in **metres**, matching what the text
        loader returns, so nothing downstream needs to know which format it came
        from.
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

        return heights_nm * 1e-9, described.scan_x_size_um

    available = list_channels(path)
    raise ValueError(
        f"{path.name} has no '{channel}' / '{direction}' plane. It contains: "
        + "; ".join(c.describe() for c in available))
