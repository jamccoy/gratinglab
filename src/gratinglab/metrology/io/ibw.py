"""
Reading Asylum Research files (.ibw) directly.

Asylum's MFP-3D and Cypher software is written in Igor Pro and saves each scan
as an Igor binary wave: a 64-byte bin header, a 320-byte wave header, the data,
and then a run of optional sections, among them a text note and the dimension
labels. An image is one 3-D wave whose layers are the channels
(``HeightRetrace``, ``ZSensorTrace``, ``AmplitudeRetrace``, ...), named in the
labels of the third dimension, with the scan parameters in the note as
``Key: value`` lines.

The reader was written here instead of taken from a library, for the reason
``spm.py`` gives: the format is small and this project depends on nothing
beyond numpy/scipy/matplotlib. It follows WaveMetrics Technical Note PTN003
for version 5 of the format, the version Asylum writes; versions 1 to 3 are
recognised and refused by name.

**It has not been checked against a file from an instrument.** None was
available when it was written (2026-10-24), so the tests build files to the
specification, which tests the reader against its author's reading of that
specification and nothing more. Two details of the specification were not
certain at the time, and the reader is written so that neither can misplace
data: the bytes after the wave are located from the header's ``wfmSize``
rather than from an assumed padding, and the header checksum is not used to
accept or refuse a file. The first real file is what turns this paragraph
into a verification, as the Gwyddion export did for ``.spm``.

The pitfall that ``spm.py`` guards against, a height scale off by a constant
factor with nothing visibly wrong, has an equivalent here: Asylum stores heights
in metres, so the data are used as they are. A topography channel whose values
reach a millimetre is refused rather than read, because that is a channel stored
in other units, and reading it as metres would scale every height silently.

Orientation: Igor stores a matrix with its first dimension varying fastest. For
an Asylum image that is X, the fast scan axis, so layer ``l`` reshaped C-order
as ``(lines, points)`` gives rows that are scan lines, which is what the text
and ``.spm`` loaders return.
"""
from __future__ import annotations

import re
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = ["IbwChannel", "read_ibw", "list_channels", "is_igor_wave_file",
           "TOPOGRAPHY_PREFERENCE", "DEFAULT_DIRECTION"]

#: Topography channels in order of preference. ZSensor is the closed-loop
#: sensor reading, the counterpart of Nanoscope's "Height Sensor"; Height is
#: derived from the piezo drive and is the fallback when no sensor channel was
#: recorded.
TOPOGRAPHY_PREFERENCE = ("ZSensor", "Height")

DEFAULT_DIRECTION = "Retrace"

_BIN_HEADER = 64
_WAVE_HEADER = 320
_DATA_START = _BIN_HEADER + _WAVE_HEADER

#: Igor number types. The complex bit (1) is refused rather than mapped.
_TYPES = {2: "f4", 4: "f8", 8: "i1", 0x10: "i2", 0x20: "i4",
          0x48: "u1", 0x50: "u2", 0x60: "u4"}

#: Each dimension label is a fixed 32-byte field (MAX_DIM_NAME + 1).
_LABEL_BYTES = 32

#: A height in metres beyond this is not a height in metres.
_MAX_HEIGHT_M = 1e-3


@dataclass(frozen=True)
class IbwChannel:
    """One layer of the image wave."""

    channel: str
    direction: str
    layer: int
    rows: int
    cols: int
    scan_x_size_um: float

    @property
    def label(self) -> str:
        return f"{self.channel}{self.direction}"

    def describe(self) -> str:
        return (f"{self.channel} / {self.direction} "
                f"({self.rows}x{self.cols}, {self.scan_x_size_um:g} um)")


@dataclass(frozen=True)
class _Wave:
    order: str          # "<" or ">"
    dims: tuple[int, int, int, int]
    dtype: np.dtype
    sf_a: tuple[float, float, float, float]
    data_units: str
    wfm_size: int
    sizes: dict


def _parse(raw: bytes) -> _Wave:
    """The two headers, in whichever byte order makes them consistent."""
    if len(raw) < _DATA_START:
        raise ValueError("too short to be an Igor binary wave")
    for order in ("<", ">"):
        version = struct.unpack_from(order + "h", raw, 0)[0]
        if version in (1, 2, 3):
            raise ValueError(
                f"Igor binary wave version {version} is not supported; Asylum "
                "software writes version 5")
        if version != 5:
            continue
        (_, _, wfm_size, formula, note, data_eunits, *rest) = struct.unpack_from(
            order + "hhiiii4i4iiii", raw, 0)
        dim_eunits, dim_labels = rest[0:4], rest[4:8]
        indices, options1, options2 = rest[8:11]
        npnts = struct.unpack_from(order + "i", raw, _BIN_HEADER + 12)[0]
        kind = struct.unpack_from(order + "h", raw, _BIN_HEADER + 16)[0]
        dims = struct.unpack_from(order + "4i", raw, _BIN_HEADER + 68)
        sf_a = struct.unpack_from(order + "4d", raw, _BIN_HEADER + 84)
        units = raw[_BIN_HEADER + 148:_BIN_HEADER + 152].split(b"\0")[0]

        if kind & 1:
            raise ValueError("complex waves are not AFM images")
        if kind not in _TYPES:
            continue
        dtype = np.dtype(order + _TYPES[kind])
        used = [d for d in dims if d > 0]
        if not used or int(np.prod(used)) != npnts:
            continue
        if wfm_size < _WAVE_HEADER + npnts * dtype.itemsize:
            continue
        return _Wave(
            order=order, dims=tuple(dims), dtype=dtype, sf_a=tuple(sf_a),
            data_units=units.decode("latin-1"), wfm_size=wfm_size,
            sizes={"formula": formula, "note": note, "data_eunits": data_eunits,
                   "dim_eunits": tuple(dim_eunits),
                   "dim_labels": tuple(dim_labels), "indices": indices,
                   "options": (options1, options2)},
        )
    raise ValueError("not an Igor binary wave (no consistent version 5 header)")


def is_igor_wave_file(path) -> bool:
    """True when the content is an Igor binary wave this module can describe.

    Like ``spm.is_nanoscope_file`` it reads the content and ignores the
    extension. A version 1 to 3 wave also answers True, so that loading one
    reports the unsupported version instead of falling through to the text
    reader.
    """
    try:
        with open(path, "rb") as handle:
            head = handle.read(_DATA_START)
    except OSError:
        return False
    try:
        _parse(head)
    except ValueError as refused:
        return "not supported" in str(refused) or "complex" in str(refused)
    return True


def _sections(raw: bytes, wave: _Wave) -> tuple[str, list[list[str]]]:
    """The note, and the labels of each dimension, from after the data."""
    at = _BIN_HEADER + wave.wfm_size
    at += wave.sizes["formula"]
    note = raw[at:at + wave.sizes["note"]].decode("latin-1", errors="replace")
    # Igor separates note lines with a bare carriage return.
    note = note.replace("\r\n", "\n").replace("\r", "\n")
    at += wave.sizes["note"] + wave.sizes["data_eunits"]
    at += sum(wave.sizes["dim_eunits"])
    labels: list[list[str]] = []
    for size in wave.sizes["dim_labels"]:
        block = raw[at:at + size]
        at += size
        labels.append([
            block[i:i + _LABEL_BYTES].split(b"\0")[0].decode("latin-1")
            for i in range(0, len(block) - _LABEL_BYTES + 1, _LABEL_BYTES)
        ])
    return note, labels


def _note_value(note: str, key: str) -> str | None:
    match = re.search(rf"^{re.escape(key)}:\s*(.+?)\s*$", note, re.M)
    return match.group(1) if match else None


def _scan_x_size_um(wave: _Wave, note: str) -> float:
    """Scan width from the X scaling, or from the note's ``ScanSize``.

    Asylum sets the wave's X scaling to metres per point. A wave left with
    Igor's default scaling of 1 per point has no physical width, and then the
    note's ``ScanSize`` (metres) is the only source.
    """
    delta = wave.sf_a[0]
    if 0.0 < delta < 1e-3:
        return delta * wave.dims[0] * 1e6
    scan = _note_value(note, "ScanSize")
    if scan is not None:
        return float(scan) * 1e6
    raise ValueError("no physical X scaling and no ScanSize in the note")


def _split_label(label: str) -> tuple[str, str]:
    for direction in ("Retrace", "Trace"):
        if label.endswith(direction) and len(label) > len(direction):
            return label[: -len(direction)], direction
    return label, "Unknown"


def _describe(raw: bytes) -> tuple[_Wave, str, list[IbwChannel]]:
    wave = _parse(raw)
    note, labels = _sections(raw, wave)
    cols, rows = wave.dims[0], wave.dims[1]
    if rows <= 0:
        raise ValueError("a one-dimensional wave is not an image")
    layers = max(wave.dims[2], 1)
    names = labels[2][1:] if len(labels) > 2 and len(labels[2]) > 1 else []
    width = _scan_x_size_um(wave, note)
    channels = []
    for layer in range(layers):
        label = names[layer] if layer < len(names) and names[layer] else f"layer{layer}"
        channel, direction = _split_label(label)
        channels.append(IbwChannel(channel, direction, layer, rows, cols, width))
    return wave, note, channels


def list_channels(path) -> tuple[IbwChannel, ...]:
    """Every layer of the image, in file order."""
    return tuple(_describe(Path(path).read_bytes())[2])


def _choose(channels: list[IbwChannel], channel: str | None,
            direction: str) -> IbwChannel:
    if channel is not None:
        wanted = [c for c in channels if c.channel == channel
                  and c.direction == direction]
        if not wanted:
            known = ", ".join(c.label for c in channels)
            raise ValueError(f"no {channel} / {direction} channel; this file has "
                             f"{known}")
        return wanted[0]
    for preferred in TOPOGRAPHY_PREFERENCE:
        for c in channels:
            if c.channel == preferred and c.direction == direction:
                return c
    known = ", ".join(c.label for c in channels)
    raise ValueError(
        f"no topography channel ({' or '.join(TOPOGRAPHY_PREFERENCE)}) in the "
        f"{direction} direction; this file has {known}")


def read_ibw(path, channel: str | None = None,
             direction: str = DEFAULT_DIRECTION):
    """
    Read one channel of an Asylum image.

    Parameters:
        channel: e.g. "ZSensor" or "Height", without the direction suffix.
            ``None`` takes the first of :data:`TOPOGRAPHY_PREFERENCE` present.
        direction: "Retrace" or "Trace".

    Returns:
        (data, scan_x_size_um), with data in metres, shaped (lines, points),
        matching what the text and ``.spm`` loaders return.
    """
    path = Path(path)
    raw = path.read_bytes()
    wave, _note, channels = _describe(raw)
    chosen = _choose(channels, channel, direction)

    count = int(np.prod([d for d in wave.dims if d > 0]))
    end = _DATA_START + count * wave.dtype.itemsize
    if end > len(raw):
        raise ValueError(f"{path.name} is truncated: the wave needs {end} bytes, "
                         f"the file has {len(raw)}")
    values = np.frombuffer(raw, dtype=wave.dtype, count=count, offset=_DATA_START)
    layers = max(wave.dims[2], 1)
    cube = values.reshape(layers, chosen.rows, chosen.cols)
    data = np.array(cube[chosen.layer], dtype=np.float64)

    if wave.data_units and wave.data_units not in ("m", ""):
        raise ValueError(
            f"{chosen.describe()} is in '{wave.data_units}'; heights are read "
            "as metres and this file does not store them that way")
    finite = data[np.isfinite(data)]
    if finite.size and np.abs(finite).max() > _MAX_HEIGHT_M:
        raise ValueError(
            f"{chosen.describe()} reaches {np.abs(finite).max():.3g}, which is "
            "not a height in metres; Asylum stores topography in metres, so "
            "this channel is in other units and reading it as metres would "
            "scale every height")
    return data, chosen.scan_x_size_um
