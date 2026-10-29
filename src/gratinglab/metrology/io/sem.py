"""
Reading Zeiss SmartSEM electron micrographs (.tif) directly.

A SmartSEM image is a baseline TIFF carrying the whole instrument state in a
private ASCII tag, 34118. The reader exists because an SEM micrograph of a
grating is a *measurement* and not a picture: the tag records the pixel size
in metres, so a cross-section can be turned into nanometres without anyone
reading a scale bar off the image.

What that buys the rest of the package is stated in ``docs/findings.md``
("What a cross-section shows that a scan cannot"): a cleaved cross-section is
the one groove measurement with no tip to deconvolve, it is the only one that
can see an undercut, and a top-down view measures the width of a crest, which
is exactly the feature an AFM tip rounds away.

The reader is written here, rather than taken from Pillow or tifffile, for the
reason given in :mod:`~gratinglab.metrology.io.spm`: the format actually in
front of us is small, and the package depends on nothing beyond
numpy/scipy/matplotlib. Only what SmartSEM writes is supported --
single-strip, uncompressed, 8-bit, one sample per pixel -- and anything else
raises rather than guessing.

Four traps are handled here because each produces a wrong number rather than
an error.

**The tag's full-precision pixel size is the display pixel, not the image
pixel.** They are equal only at a store resolution of 1024 x 768. Scan at
2048 x 1536 and the number is twice too large, so every length measured from
the image doubles while the image still looks entirely reasonable. The field
``AP_IMAGE_PIXEL_SIZE`` is the one about the stored image, and
:func:`_image_pixel_size_nm` uses it to anchor the scale before taking the
block's extra digits.

**The palette is not a greyscale ramp.** The images are TIFF palette images
whose colour table is the identity (entry ``i`` is grey ``i``) at every index
except multiples of nine, which SmartSEM reserves for the coloured annotation
in the data bar. The image itself never uses the reserved indices, so the raw
index *is* the detector level; mapping it through the palette would be
harmless but converting the result to luminance would quietly recolour the
reserved entries into plausible grey values.

**The data bar overwrites the image.** The bar is painted over the bottom of
the scan, not appended to it, so the stated field of view covers the full
stored height and the pixels beneath the bar are gone. Measuring a depth
across it would include the bar's own geometry. :func:`read_sem` detects it
and crops it off, and records how many rows it removed. Some images carry no
bar at all, and those are returned whole.

**A tilted stage foreshortens the image.** With the stage tilted by ``T`` and
SmartSEM's tilt correction off, a length across the tilt axis appears shorter
than it is by ``cos T``, and the recorded pixel size applies to the image, not
to the specimen. With the correction on, the scan has already been stretched
and the recorded pixel size is right as it stands. The reader reports
:attr:`SemImage.foreshortening` and refuses to apply it silently, in the same
spirit as the ``Sens. ZsensSens`` trap in the Nanoscope reader: an image taken
at 45 degrees with the correction off understates a depth by 29%, and looks
perfectly reasonable while doing it.
"""
from __future__ import annotations

import re
import struct
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

__all__ = [
    "SemImage",
    "read_sem",
    "read_sem_tags",
    "is_zeiss_sem_file",
    "parse_zeiss_value",
    "ZEISS_TAG",
]

#: The private TIFF tag holding SmartSEM's instrument state, as latin-1 text.
#: Tag 34119 holds the same content as UTF-16 and is ignored.
ZEISS_TAG = 34118

_TIFF_MAGIC = 42

#: TIFF field types, as (struct code, bytes). Only those SmartSEM writes.
_FIELD_TYPES = {
    1: ("B", 1),   # BYTE
    2: ("s", 1),   # ASCII
    3: ("H", 2),   # SHORT
    4: ("I", 4),   # LONG
    5: ("II", 8),  # RATIONAL
    7: ("B", 1),   # UNDEFINED
}

# TIFF tags this reader uses.
_IMAGE_WIDTH = 256
_IMAGE_LENGTH = 257
_BITS_PER_SAMPLE = 258
_COMPRESSION = 259
_STRIP_OFFSETS = 273
_SAMPLES_PER_PIXEL = 277
_STRIP_BYTE_COUNTS = 279
_PLANAR_CONFIG = 284

#: Fraction of a row that must share one value for it to be the rule SmartSEM
#: draws between the micrograph and the data bar. Detector noise puts a real
#: scan far below this even over a featureless field, so the rule is
#: unmistakable; a few stray pixels of overlapping annotation are tolerated.
_SEPARATOR_ROW_FRACTION = 0.98

#: Fewest rows a data bar can occupy. Guards against a single saturated or
#: blanked scan line at the foot of an image being read as a bar.
_MIN_DATABAR_ROWS = 8

#: How far up from the bottom edge to look for the bar.
_DATABAR_SEARCH_FRACTION = 0.2

#: Length units SmartSEM writes, in nanometres.
_LENGTH_UNITS_NM = {
    "fm": 1e-6, "pm": 1e-3, "nm": 1.0, "um": 1e3, "µm": 1e3,
    "μm": 1e3, "mm": 1e6, "cm": 1e7, "m": 1e9,
}

#: Other unit suffixes, as a multiplier on the bare number.
_SCALAR_UNITS = {
    "kv": 1e3, "v": 1.0, "mv": 1e-3,
    "k x": 1e3, "m x": 1e6, "x": 1.0,
    "°": 1.0, "%": 1.0, "ka": 1e3, "a": 1.0, "pa": 1.0,
    "ns": 1.0, "us": 1e3, "ms": 1e6, "s": 1.0,
}

_NUMBER = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

#: A SmartSEM field key: two letters, an underscore, then upper-case and
#: digits. ``AP_`` is an analogue parameter, ``DP_`` a digital one and
#: ``SV_`` a stored value. No field *value* has this shape -- every one
#: carries a label and a separator -- so a key is unambiguous.
_KEY = re.compile(r"[A-Z]{2}_[A-Z0-9_]+")


class _TiffError(ValueError):
    """The file is not the kind of TIFF this reader supports."""


def _read_ifd(blob: bytes) -> tuple[dict[int, object], str]:
    """Tag values of the first IFD, and the byte-order code."""
    if len(blob) < 8:
        raise _TiffError("file is too short to be a TIFF")
    order = blob[:2]
    if order == b"II":
        bo = "<"
    elif order == b"MM":
        bo = ">"
    else:
        raise _TiffError(f"not a TIFF: byte-order mark is {order!r}")

    magic, ifd_offset = struct.unpack(bo + "HI", blob[2:8])
    if magic != _TIFF_MAGIC:
        raise _TiffError(f"not a TIFF: magic number is {magic}, expected 42")

    (count,) = struct.unpack(bo + "H", blob[ifd_offset : ifd_offset + 2])
    tags: dict[int, object] = {}
    for i in range(count):
        entry = ifd_offset + 2 + 12 * i
        tag, field_type, n = struct.unpack(bo + "HHI", blob[entry : entry + 8])
        payload = blob[entry + 8 : entry + 12]
        if field_type not in _FIELD_TYPES:
            continue  # A tag this reader has no use for, in a type it cannot read.
        code, size = _FIELD_TYPES[field_type]
        total = n * size
        if total > 4:
            (offset,) = struct.unpack(bo + "I", payload)
            raw = blob[offset : offset + total]
        else:
            raw = payload[:total]
        if len(raw) < total:
            raise _TiffError(f"tag {tag} runs past the end of the file")

        if field_type == 2:
            tags[tag] = raw.rstrip(b"\x00").decode("latin-1")
        elif field_type == 5:
            pairs = struct.unpack(bo + code * n, raw)
            tags[tag] = tuple(
                num / den if den else 0.0
                for num, den in zip(pairs[0::2], pairs[1::2])
            )
        else:
            values = struct.unpack(bo + code * n, raw)
            tags[tag] = values[0] if n == 1 else values
    return tags, bo


def _tag(tags: dict[int, object], key: int, default: object = None) -> object:
    value = tags.get(key, default)
    return value[0] if isinstance(value, tuple) and len(value) == 1 else value


def _read_pixels(blob: bytes, tags: dict[int, object], bo: str) -> NDArray[np.uint8]:
    """The single uncompressed 8-bit plane, as stored."""
    width = int(_tag(tags, _IMAGE_WIDTH))
    height = int(_tag(tags, _IMAGE_LENGTH))
    bits = int(_tag(tags, _BITS_PER_SAMPLE, 8))
    compression = int(_tag(tags, _COMPRESSION, 1))
    samples = int(_tag(tags, _SAMPLES_PER_PIXEL, 1))
    planar = int(_tag(tags, _PLANAR_CONFIG, 1))

    if compression != 1:
        raise _TiffError(
            f"compressed TIFF (compression={compression}); this reader handles "
            "only the uncompressed images SmartSEM writes"
        )
    if bits != 8 or samples != 1 or planar != 1:
        raise _TiffError(
            f"unsupported layout: {bits} bits, {samples} sample(s) per pixel, "
            f"planar configuration {planar}. SmartSEM writes 8-bit single-sample "
            "images and this reader supports no other"
        )

    offsets = tags.get(_STRIP_OFFSETS)
    counts = tags.get(_STRIP_BYTE_COUNTS)
    if offsets is None or counts is None:
        raise _TiffError("no strip offsets; tiled TIFFs are not supported")
    if not isinstance(offsets, tuple):
        offsets, counts = (offsets,), (counts,)

    data = b"".join(
        blob[int(o) : int(o) + int(c)] for o, c in zip(offsets, counts)
    )
    if len(data) < width * height:
        raise _TiffError(
            f"image data is {len(data)} bytes, short of the {width * height} "
            f"a {width}x{height} 8-bit image needs"
        )
    return np.frombuffer(data[: width * height], dtype=np.uint8).reshape(height, width)


def _databar_rows(pixels: NDArray[np.uint8]) -> int:
    """How many rows at the foot of the image the data bar occupies.

    The bar is found by the rule SmartSEM draws above it: a full-width row of
    one saturated value, pure black or pure white, in the bottom fifth of the
    frame. The bar's own contents are not a reliable signal, because the two
    instruments draw it differently. A Gemini 500 may fill it with a flat
    background and coloured annotation or leave the micrograph showing
    through; a Merlin always leaves it showing through and writes the
    annotation over it, so most of its bar is indistinguishable from specimen
    row by row. The rule above it is common to all of them, and its absence is
    how an image saved with the bar switched off is recognised.

    Both halves of the test earn their place. Dropping the saturation
    requirement admits the faint full-width scan artefacts that cross some
    micrographs, which would crop away good specimen; dropping the
    full-width requirement admits a charging bright edge, which can put a
    tenth of a row at 255 without being a rule at all.
    """
    height, width = pixels.shape
    first = int(height * (1.0 - _DATABAR_SEARCH_FRACTION))
    for row in range(first, height - _MIN_DATABAR_ROWS):
        counts = np.bincount(pixels[row], minlength=256)
        if counts.max() / width >= _SEPARATOR_ROW_FRACTION and counts.argmax() in (
            0,
            255,
        ):
            return height - row
    return 0


def parse_zeiss_value(text: str) -> float | None:
    """The number in a SmartSEM field, in nm, volts, degrees or bare multiples.

    Fields look like ``'WD =  3.6 mm'`` or ``'Mag =   50.00 K X'``. Lengths
    are converted to nanometres, so the caller never has to notice that the
    pixel size of one image is quoted in nm and another's in µm. A field with
    no number in it returns ``None``.

    >>> parse_zeiss_value('Image Pixel Size = 2.233 nm')
    2.233
    >>> parse_zeiss_value('WD =  3.6 mm')
    3600000.0
    >>> parse_zeiss_value('Mag =   50.00 K X')
    50000.0
    >>> parse_zeiss_value('Sample ID =') is None
    True
    """
    _, _, value = text.partition("=" if "=" in text else ":")
    value = value.strip()
    match = _NUMBER.search(value)
    if match is None:
        return None
    number = float(match.group())
    unit = value[match.end() :].strip()
    if unit in _LENGTH_UNITS_NM:
        return number * _LENGTH_UNITS_NM[unit]
    return number * _SCALAR_UNITS.get(unit.lower(), 1.0)


def _parse_zeiss_tag(text: str) -> tuple[dict[str, str], float | None]:
    """SmartSEM's instrument state: the keyed fields, and the display pixel size.

    The tag opens with three zeros and then four seven-number blocks, one per
    display zone, each separated by its index. The first number of each block
    is a pixel size in metres at full precision. After the blocks comes a
    count, and then that many ``key`` / ``label = value`` line pairs.

    The returned pixel size is the *display* one, which is what the block
    holds; :func:`_image_pixel_size_nm` turns it into the image one.
    """
    lines = text.split("\r\n")
    exact_pixel_m: float | None = None
    try:
        exact_pixel_m = float(lines[3])
    except (IndexError, ValueError):
        exact_pixel_m = None

    # The records start at the first SmartSEM key. Counting in from the top
    # instead would be wrong: the per-zone blocks are separated by their own
    # index digits, so the count is not the first bare integer in the tag, and
    # the number of leading blocks is not fixed across SmartSEM versions.
    fields: dict[str, str] = {}
    start = next(
        (i for i, line in enumerate(lines) if _KEY.fullmatch(line)), None
    )
    if start is not None:
        records = lines[start:]
        for key, value in zip(records[0::2], records[1::2]):
            if _KEY.fullmatch(key):
                fields[key] = value

    if exact_pixel_m is not None and not 0.0 < exact_pixel_m < 1.0:
        exact_pixel_m = None
    return fields, exact_pixel_m


def _image_pixel_size_nm(
    fields: dict[str, str], display_nm: float | None
) -> float | None:
    """The edge of one pixel *of the stored image*, in nm.

    Three sources disagree, and only one of them is the right answer:

    ``AP_IMAGE_PIXEL_SIZE``
        the image pixel, rounded to four significant figures.
    ``AP_PIXEL_SIZE``
        the *display* pixel. Equal to the image pixel only when the store
        resolution matches the display's 1024 x 768.
    the leading numeric block
        the display pixel again, but at full precision.

    Scanning at 2048 x 1536 therefore makes the block exactly twice the real
    image pixel, and every length measured from such an image comes out
    doubled with nothing whatever to show for it -- the image is the right
    shape, the grooves are in the right places, and the period is wrong by a
    factor of two. The quoted image pixel is what settles the scale; the
    block then supplies the precision it was rounded away from, but only
    after the integer ratio between them has been confirmed.
    """
    quoted = fields.get("AP_IMAGE_PIXEL_SIZE")
    image_nm = parse_zeiss_value(quoted) if quoted else None
    if image_nm is None or image_nm <= 0:
        # No image pixel quoted: the display pixel is all there is, and on a
        # store resolution of 1024 x 768 it is also the right answer.
        fallback = fields.get("AP_PIXEL_SIZE")
        return display_nm or (parse_zeiss_value(fallback) if fallback else None)
    if display_nm is None or display_nm <= 0:
        return image_nm

    ratio = display_nm / image_nm
    binning = round(ratio) if ratio >= 1.0 else 1.0 / round(1.0 / ratio)
    if binning and abs(ratio - binning) <= 0.01 * max(binning, 1.0):
        return display_nm / binning
    # The two sources are not related by a whole number of display pixels per
    # image pixel, so the assumption behind the refinement does not hold here.
    # Keep the quoted value, which is coarse but is about the right quantity.
    return image_nm


@dataclass(frozen=True)
class SemImage:
    """One micrograph, with the instrument state that makes it a measurement.

    ``data`` holds detector levels on the instrument's own 0--255 scale, which
    the reader does not rescale: the grey level of an SEM image is a yield,
    not a height, and normalising it would imply otherwise.
    """

    #: Detector levels, rows top to bottom, with the data bar removed.
    data: NDArray[np.float64]

    #: Edge length of a pixel in the image plane, nm.
    pixel_size_nm: float

    #: SmartSEM's fields, keyed as the instrument keys them (``'AP_WD'``), each
    #: holding the whole labelled line (``'WD =  3.6 mm'``).
    fields: dict[str, str]

    #: Rows removed from the foot of the image; 0 when it carried no data bar.
    databar_rows: int

    path: Path

    def field(self, key: str) -> str | None:
        """The labelled line for a SmartSEM key, or ``None`` if absent."""
        return self.fields.get(key)

    def value(self, key: str) -> float | None:
        """The number in a SmartSEM field: nm for lengths, else as written."""
        text = self.fields.get(key)
        return None if text is None else parse_zeiss_value(text)

    def text(self, key: str) -> str | None:
        """The part of a field after its label, unparsed."""
        line = self.fields.get(key)
        if line is None:
            return None
        _, _, value = line.partition("=" if "=" in line else ":")
        return value.strip()

    @property
    def shape(self) -> tuple[int, int]:
        return self.data.shape

    @property
    def width_nm(self) -> float:
        return self.data.shape[1] * self.pixel_size_nm

    @property
    def height_nm(self) -> float:
        """Extent of the returned image, which excludes any cropped data bar."""
        return self.data.shape[0] * self.pixel_size_nm

    @property
    def magnification(self) -> float | None:
        return self.value("AP_MAG")

    @property
    def eht_kv(self) -> float | None:
        """Accelerating voltage in kV. Stored in volts, reported in kV."""
        volts = self.value("AP_ACTUALKV")
        return None if volts is None else volts / 1e3

    @property
    def working_distance_mm(self) -> float | None:
        nanometres = self.value("AP_WD")
        return None if nanometres is None else nanometres / 1e6

    @property
    def stage_tilt_deg(self) -> float | None:
        """Physical tilt of the stage, degrees. Distinct from the correction."""
        return self.value("AP_STAGE_AT_T")

    @property
    def tilt_correction(self) -> bool | None:
        """Whether SmartSEM was stretching the scan to undo the foreshortening."""
        text = self.text("DP_TILT_CORRECTION")
        return None if text is None else text.strip().lower() == "on"

    @property
    def tilt_angle_deg(self) -> float | None:
        """The angle entered for the tilt correction, not the stage's own tilt."""
        return self.value("AP_TILT_ANGLE")

    @property
    def tilt_axis_deg(self) -> float | None:
        """Azimuth of the tilt axis in the image, degrees."""
        return self.value("AP_TILT_AXIS")

    @property
    def detector(self) -> str | None:
        return self.text("DP_DETECTOR_CHANNEL") or self.text("AP_DETECTOR")

    @property
    def instrument(self) -> str | None:
        return self.text("SV_SERIAL_NUMBER")

    @property
    def operator(self) -> str | None:
        return self.text("SV_USER_NAME")

    @property
    def recorded(self) -> datetime | None:
        """When the image was taken, or ``None`` if either field is unreadable."""
        date, time = self.text("AP_DATE"), self.text("AP_TIME")
        if not date or not time:
            return None
        for fmt in ("%d %b %Y %H:%M:%S", "%d-%b-%Y %H:%M:%S", "%d %B %Y %H:%M:%S"):
            try:
                return datetime.strptime(f"{date} {time}", fmt)
            except ValueError:
                continue
        return None

    @property
    def foreshortening(self) -> float:
        """Factor by which the image compresses lengths across the tilt axis.

        ``cos`` of the stage tilt when SmartSEM's tilt correction was off, and
        1.0 when it was on or the stage was flat. A length measured across the
        tilt axis is the image length divided by this.

        It is reported and not applied. Which image direction lies across the
        tilt axis depends on :attr:`tilt_axis_deg` and on how the specimen was
        loaded, and a reader who has not thought about that should not be
        handed a silently stretched image.
        """
        tilt = self.stage_tilt_deg
        if tilt is None or self.tilt_correction:
            return 1.0
        return float(np.cos(np.radians(tilt)))

    @property
    def summary(self) -> str:
        """One-line description, for a status bar."""
        rows, cols = self.data.shape
        bar = f", {self.databar_rows} data-bar rows cropped" if self.databar_rows else ""
        return (
            f"{cols}x{rows} at {self.pixel_size_nm:.4g} nm/pixel "
            f"({self.width_nm / 1e3:.3g} x {self.height_nm / 1e3:.3g} um){bar}"
        )


def _load(path: Path) -> tuple[bytes, dict[int, object], str]:
    blob = Path(path).read_bytes()
    tags, bo = _read_ifd(blob)
    return blob, tags, bo


def is_zeiss_sem_file(path: str | Path) -> bool:
    """True when the file is a TIFF carrying SmartSEM's instrument tag."""
    try:
        _, tags, _ = _load(Path(path))
    except (OSError, _TiffError, struct.error):
        return False
    return ZEISS_TAG in tags


def read_sem_tags(path: str | Path) -> dict[str, str]:
    """SmartSEM's instrument state, without decoding the image.

    Keyed as the instrument keys it, each value the whole labelled line::

        >>> tags = read_sem_tags(path)                      # doctest: +SKIP
        >>> tags['AP_IMAGE_PIXEL_SIZE']                     # doctest: +SKIP
        'Image Pixel Size = 2.233 nm'
    """
    _, tags, _ = _load(Path(path))
    text = tags.get(ZEISS_TAG)
    if not isinstance(text, str):
        raise ValueError(
            f"{path} has no SmartSEM tag {ZEISS_TAG}; it is a TIFF but not a "
            "Zeiss micrograph"
        )
    fields, _ = _parse_zeiss_tag(text)
    return fields


def read_sem(path: str | Path, *, crop_databar: bool = True) -> SemImage:
    """Read a SmartSEM micrograph and the instrument state that scales it.

    Parameters
    ----------
    path
        A Zeiss SmartSEM TIFF.
    crop_databar
        Remove the burnt-in data bar from the foot of the image, which is the
        default because the bar is not specimen. Pass ``False`` to get the
        frame exactly as stored.

    Raises
    ------
    ValueError
        If the file is not a SmartSEM TIFF, or records no pixel size. A
        micrograph without a pixel size is a picture, and this package has no
        use for one.
    """
    path = Path(path)
    blob, tags, bo = _load(path)

    text = tags.get(ZEISS_TAG)
    if not isinstance(text, str):
        raise ValueError(
            f"{path} has no SmartSEM tag {ZEISS_TAG}; it is a TIFF but not a "
            "Zeiss micrograph"
        )
    fields, exact_pixel_m = _parse_zeiss_tag(text)
    display_nm = None if exact_pixel_m is None else exact_pixel_m * 1e9
    pixel_size_nm = _image_pixel_size_nm(fields, display_nm)
    if not pixel_size_nm or pixel_size_nm <= 0:
        raise ValueError(
            f"{path} records no usable pixel size; nothing can be measured "
            "from it"
        )

    pixels = _read_pixels(blob, tags, bo)
    rows = _databar_rows(pixels) if crop_databar else 0
    kept = pixels[: pixels.shape[0] - rows] if rows else pixels

    return SemImage(
        data=kept.astype(np.float64),
        pixel_size_nm=float(pixel_size_nm),
        fields=fields,
        databar_rows=int(rows),
        path=path,
    )
