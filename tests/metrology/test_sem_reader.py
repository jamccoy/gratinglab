"""
Reading Zeiss SmartSEM micrographs.

Unlike the Nanoscope reader next door, this one tests against files it builds
itself. The .spm tests need a real instrument file because the trap there is a
pair of similarly named header keys that only a real file contains. Here the
traps are all in geometry and arithmetic -- which of two pixel sizes is the
image's, where the data bar begins, what a tilted stage does to a length -- and
each can be set up exactly, with the right answer known in advance, in a file
of a few kilobytes. A synthetic fixture also lets the suite assert on cases the
group's own micrographs do not happen to contain.

The builder below writes the same thing SmartSEM writes: a little-endian
baseline TIFF, single uncompressed strip, 8-bit, with the instrument state in
private tag 34118.
"""
import os
import struct
import sys

import numpy as np
import pytest

# src/ is placed on the path by tests/conftest.py; repeated here so the file
# also runs directly as a script, not only under pytest.
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src'))

from gratinglab.metrology.io.sem import (                               # noqa: E402
    SemImage, is_zeiss_sem_file, parse_zeiss_value, read_sem, read_sem_tags)

# ── Building a SmartSEM TIFF ─────────────────────────────────────────────────

#: The fields a micrograph of the group's own gratings carries, as SmartSEM
#: writes them: key, then the whole labelled line.
BASE_FIELDS = {
    "AP_IMAGE_PIXEL_SIZE": "Image Pixel Size = 2.233 nm",
    "AP_PIXEL_SIZE": "Pixel Size = 2.233 nm",
    "AP_MAG": "Mag =   50.00 K X",
    "AP_ACTUALKV": "EHT =  5.00 kV",
    "AP_WD": "WD =  3.6 mm",
    "AP_STAGE_AT_T": "Stage at T =   0.0 °",
    "AP_TILT_ANGLE": "Tilt Angle =   0.0 °",
    "AP_TILT_AXIS": "Tilt Axis =  90.0 °",
    "DP_TILT_CORRECTION": "Tilt Corrn. = Off",
    "DP_DETECTOR_CHANNEL": "Signal A = InLens",
    "SV_SERIAL_NUMBER": "Serial No. = GeminiSEM 500-8203017174",
    "SV_USER_NAME": "User Name = JAM1117",
    "AP_DATE": "Date: 10 Jun 2021",
    "AP_TIME": "Time: 15:37:59",
}


def zeiss_tag(fields, display_pixel_m=2.232998e-09):
    """SmartSEM's tag 34118: leading per-zone blocks, a count, then the pairs."""
    lines = ["0", "0", "0"]
    for zone in range(4):
        if zone:
            lines.append(str(zone))
        lines += [f"{display_pixel_m:e}", "5.000003e+04", "6", "5.000000e+03",
                  "2.310000e+00", "0.000000e+00", "3.566873e-03"]
    lines.append(str(len(fields)))
    for key, value in fields.items():
        lines += [key, value]
    return "\r\n".join(lines) + "\r\n"


def write_tiff(path, pixels, tag_text=None, *, palette=True):
    """A baseline little-endian TIFF of the kind SmartSEM writes."""
    height, width = pixels.shape
    entries = [
        (256, 3, 1, width),            # ImageWidth
        (257, 3, 1, height),           # ImageLength
        (258, 3, 1, 8),                # BitsPerSample
        (259, 3, 1, 1),                # Compression: none
        (262, 3, 1, 3 if palette else 1),   # PhotometricInterpretation
        (277, 3, 1, 1),                # SamplesPerPixel
        (278, 4, 1, 0xFFFFFFFF),       # RowsPerStrip: one strip
        (284, 3, 1, 1),                # PlanarConfiguration
    ]
    blobs = []          # (tag, type, count, bytes) for values that do not fit

    if palette:
        # The identity ramp SmartSEM uses, with every ninth entry reserved for
        # the coloured annotation in the data bar.
        table = []
        for channel in range(3):
            for i in range(256):
                reserved = (i % 9 == 0) and i > 0
                table.append((255 if channel == 0 else 0) * 257 if reserved else i * 257)
        blobs.append((320, 3, len(table), struct.pack("<" + "H" * len(table), *table)))
    if tag_text is not None:
        raw = tag_text.encode("latin-1") + b"\x00"
        blobs.append((34118, 2, len(raw), raw))

    # +2 for StripOffsets and StripByteCounts, whose values need the layout.
    n = len(entries) + len(blobs) + 2
    cursor = 8 + 2 + 12 * n + 4
    placed = []
    for tag, field_type, count, raw in blobs:
        placed.append((tag, field_type, count, cursor))
        cursor += len(raw)
    strip_offset = cursor

    all_entries = sorted(
        entries
        + placed
        + [(273, 4, 1, strip_offset), (279, 4, 1, pixels.size)],
        key=lambda e: e[0],
    )
    out = bytearray()
    out += b"II" + struct.pack("<HI", 42, 8)
    out += struct.pack("<H", len(all_entries))
    for tag, field_type, count, value in all_entries:
        out += struct.pack("<HHI", tag, field_type, count)
        size = {1: 1, 2: 1, 3: 2, 4: 4}[field_type] * count
        if size <= 4:
            out += struct.pack("<H", value)[:2].ljust(4, b"\x00") if field_type == 3 \
                else struct.pack("<I", value)
        else:
            out += struct.pack("<I", value)
    out += struct.pack("<I", 0)                 # no next IFD
    for (_, _, _, raw), _ in zip(blobs, placed):
        out += raw
    out += pixels.astype(np.uint8).tobytes()
    path.write_bytes(bytes(out))
    return path


def micrograph(height=768, width=1024, seed=0):
    """A plausible noisy field with no saturated pixels anywhere."""
    rng = np.random.default_rng(seed)
    return rng.integers(40, 200, size=(height, width), dtype=np.uint8)


def with_databar(pixels, rows=77, *, transparent=False):
    """Paint a SmartSEM data bar over the foot, rule and all."""
    out = pixels.copy()
    top = out.shape[0] - rows
    out[top] = 0                                  # the rule above the bar
    if not transparent:
        out[top + 1:] = 30                        # a Gemini's flat background
    out[top + 20: top + 28, 100:400] = 255        # annotation, either way
    return out


@pytest.fixture
def sem_file(tmp_path):
    """A complete SmartSEM micrograph with a data bar, 2.233 nm/pixel."""
    pixels = with_databar(micrograph())
    return write_tiff(tmp_path / "xsec.tif", pixels, zeiss_tag(BASE_FIELDS))


# ── Recognition ──────────────────────────────────────────────────────────────

class TestRecognition:
    def test_recognises_a_smartsem_file_by_content(self, sem_file):
        assert is_zeiss_sem_file(sem_file)

    def test_a_tiff_without_the_tag_is_not_a_micrograph(self, tmp_path):
        """The group's own folder holds one of these: a micrograph re-saved by
        an image editor, which kept the pixels and dropped the instrument."""
        plain = write_tiff(tmp_path / "copy.tif", micrograph(64, 64), None)
        assert not is_zeiss_sem_file(plain)
        with pytest.raises(ValueError, match="not a Zeiss micrograph"):
            read_sem(plain)

    def test_a_file_that_is_not_a_tiff_is_rejected(self, tmp_path):
        path = tmp_path / "notatiff.tif"
        path.write_bytes(b"this is not a TIFF at all, not even close")
        assert not is_zeiss_sem_file(path)

    def test_a_missing_file_is_not_a_micrograph(self, tmp_path):
        assert not is_zeiss_sem_file(tmp_path / "absent.tif")


# ── The pixel size, which is the whole point ─────────────────────────────────

class TestPixelSize:
    def test_takes_full_precision_from_the_leading_block(self, sem_file):
        """The keyed field is rounded to four figures; the block is not."""
        image = read_sem(sem_file)
        assert image.pixel_size_nm == pytest.approx(2.232998, abs=1e-9)

    def test_the_block_is_the_display_pixel_and_is_scaled_to_the_image(
        self, tmp_path
    ):
        """The trap. At a store resolution of 2048 x 1536 the block holds twice
        the image pixel size, and taking it as written doubles every length
        measured from the micrograph.
        """
        fields = dict(BASE_FIELDS)
        fields["AP_IMAGE_PIXEL_SIZE"] = "Image Pixel Size = 5.582 nm"
        fields["AP_PIXEL_SIZE"] = "Pixel Size = 11.16 nm"
        path = write_tiff(
            tmp_path / "big.tif",
            micrograph(1536, 2048),
            zeiss_tag(fields, display_pixel_m=1.1165e-08),
        )
        image = read_sem(path)
        assert image.pixel_size_nm == pytest.approx(5.5825, abs=1e-6)

    def test_falls_back_to_the_quoted_field_when_the_two_disagree(self, tmp_path):
        """No whole number of display pixels per image pixel means the
        refinement's premise is false, so the coarse value stands rather than
        a confident wrong one."""
        fields = dict(BASE_FIELDS)
        fields["AP_IMAGE_PIXEL_SIZE"] = "Image Pixel Size = 5.000 nm"
        path = write_tiff(
            tmp_path / "odd.tif", micrograph(64, 64),
            zeiss_tag(fields, display_pixel_m=7.3e-09),
        )
        assert read_sem(path).pixel_size_nm == pytest.approx(5.0)

    def test_a_micrograph_without_a_pixel_size_is_refused(self, tmp_path):
        """Without a scale it is a picture, and this package has no use for
        one. Better to refuse than to return an image that cannot be measured."""
        fields = {k: v for k, v in BASE_FIELDS.items() if "PIXEL_SIZE" not in k}
        path = write_tiff(
            tmp_path / "noscale.tif", micrograph(64, 64),
            zeiss_tag(fields, display_pixel_m=0.0),
        )
        with pytest.raises(ValueError, match="no usable pixel size"):
            read_sem(path)

    def test_field_of_view_follows_from_the_pixel_size(self, sem_file):
        image = read_sem(sem_file)
        assert image.width_nm == pytest.approx(1024 * image.pixel_size_nm)
        assert image.height_nm == pytest.approx(
            (768 - image.databar_rows) * image.pixel_size_nm
        )


# ── The data bar ─────────────────────────────────────────────────────────────

class TestDataBar:
    def test_crops_a_flat_data_bar(self, sem_file):
        image = read_sem(sem_file)
        assert image.databar_rows == 77
        assert image.data.shape == (768 - 77, 1024)

    def test_crops_a_transparent_data_bar(self, tmp_path):
        """A Merlin leaves the micrograph showing through the bar and writes
        over it, so the bar is not flat and only the rule above it marks where
        the specimen ends."""
        pixels = with_databar(micrograph(), rows=63, transparent=True)
        path = write_tiff(tmp_path / "merlin.tif", pixels, zeiss_tag(BASE_FIELDS))
        assert read_sem(path).databar_rows == 63

    def test_an_image_saved_without_a_bar_is_returned_whole(self, tmp_path):
        path = write_tiff(tmp_path / "nobar.tif", micrograph(), zeiss_tag(BASE_FIELDS))
        image = read_sem(path)
        assert image.databar_rows == 0
        assert image.data.shape == (768, 1024)

    def test_a_faint_scan_artefact_is_not_mistaken_for_the_rule(self, tmp_path):
        """One of the group's own micrographs carries a full-width grey line
        across it above the bar. Cropping there would throw away good specimen,
        so the rule must be saturated as well as full width."""
        pixels = with_databar(micrograph())
        pixels[768 - 114] = 128
        path = write_tiff(tmp_path / "artefact.tif", pixels, zeiss_tag(BASE_FIELDS))
        assert read_sem(path).databar_rows == 77

    def test_a_charging_bright_edge_is_not_mistaken_for_the_rule(self, tmp_path):
        """Saturation alone is not enough either: a charged edge can put a
        tenth of a row at 255 without being a rule."""
        pixels = micrograph()
        pixels[700, :120] = 255
        path = write_tiff(tmp_path / "charged.tif", pixels, zeiss_tag(BASE_FIELDS))
        assert read_sem(path).databar_rows == 0

    def test_the_bar_can_be_kept(self, sem_file):
        image = read_sem(sem_file, crop_databar=False)
        assert image.databar_rows == 0
        assert image.data.shape == (768, 1024)

    def test_pixels_are_detector_levels_and_are_not_rescaled(self, tmp_path):
        """The grey level of an SEM image is a yield, not a height."""
        pixels = micrograph()
        path = write_tiff(tmp_path / "raw.tif", pixels, zeiss_tag(BASE_FIELDS))
        assert np.array_equal(read_sem(path).data, pixels.astype(float))


# ── The instrument state ─────────────────────────────────────────────────────

class TestInstrumentState:
    def test_reports_the_conditions_the_image_was_taken_under(self, sem_file):
        image = read_sem(sem_file)
        assert image.magnification == pytest.approx(50_000.0)
        assert image.eht_kv == pytest.approx(5.0)
        assert image.working_distance_mm == pytest.approx(3.6)
        assert image.detector == "InLens"
        assert image.instrument == "GeminiSEM 500-8203017174"
        assert image.operator == "JAM1117"
        assert image.recorded.isoformat() == "2021-06-10T15:37:59"

    def test_tags_can_be_read_without_decoding_the_image(self, sem_file):
        fields = read_sem_tags(sem_file)
        assert fields["AP_WD"] == "WD =  3.6 mm"
        assert len(fields) == len(BASE_FIELDS)

    def test_an_absent_field_is_none_rather_than_a_guess(self, tmp_path):
        """SmartSEM versions differ in which keys they write."""
        fields = {k: v for k, v in BASE_FIELDS.items() if k != "AP_WD"}
        path = write_tiff(
            tmp_path / "sparse.tif", micrograph(64, 64), zeiss_tag(fields)
        )
        image = read_sem(path)
        assert image.working_distance_mm is None
        assert image.magnification == pytest.approx(50_000.0)

    def test_summary_states_the_scale_and_what_was_cropped(self, sem_file):
        summary = read_sem(sem_file).summary
        assert "2.233 nm/pixel" in summary
        assert "77 data-bar rows cropped" in summary


# ── The tilt, reported and not applied ───────────────────────────────────────

class TestForeshortening:
    def test_a_flat_stage_foreshortens_nothing(self, sem_file):
        assert read_sem(sem_file).foreshortening == pytest.approx(1.0)

    def test_a_tilted_stage_foreshortens_by_the_cosine(self, tmp_path):
        fields = dict(BASE_FIELDS)
        fields["AP_STAGE_AT_T"] = "Stage at T =  45.0 °"
        path = write_tiff(
            tmp_path / "tilted.tif", micrograph(64, 64), zeiss_tag(fields)
        )
        image = read_sem(path)
        assert image.stage_tilt_deg == pytest.approx(45.0)
        assert image.foreshortening == pytest.approx(np.cos(np.radians(45.0)))

    def test_the_instruments_own_correction_is_honoured(self, tmp_path):
        """With the correction on, SmartSEM already stretched the scan and the
        recorded pixel size is right as it stands."""
        fields = dict(BASE_FIELDS)
        fields["AP_STAGE_AT_T"] = "Stage at T =  45.0 °"
        fields["DP_TILT_CORRECTION"] = "Tilt Corrn. = On"
        path = write_tiff(
            tmp_path / "corrected.tif", micrograph(64, 64), zeiss_tag(fields)
        )
        image = read_sem(path)
        assert image.tilt_correction is True
        assert image.foreshortening == pytest.approx(1.0)

    def test_the_pixel_size_is_not_quietly_stretched(self, tmp_path):
        """The reader reports the foreshortening and leaves it to the caller.
        Which image direction lies across the tilt axis depends on how the
        specimen was loaded, which the file does not say."""
        fields = dict(BASE_FIELDS)
        fields["AP_STAGE_AT_T"] = "Stage at T =  45.0 °"
        path = write_tiff(
            tmp_path / "tilted2.tif", micrograph(64, 64), zeiss_tag(fields)
        )
        image = read_sem(path)
        assert image.pixel_size_nm == pytest.approx(2.232998, abs=1e-9)


# ── Unit parsing ─────────────────────────────────────────────────────────────

class TestValueParsing:
    @pytest.mark.parametrize(
        "text, expected",
        [
            ("Image Pixel Size = 2.233 nm", 2.233),
            ("Width = 2.287 µm", 2287.0),
            ("WD =  3.6 mm", 3.6e6),
            ("FIB Pixel Size = 0.0 fm", 0.0),
            ("Mag =   50.00 K X", 50_000.0),
            ("Min Mag =    157 X", 157.0),
            ("EHT =  5.00 kV", 5000.0),
            ("Stage at T =  -0.8 °", -0.8),
            ("Brightness =  50.0 %", 50.0),
            ("Date: 10 Jun 2021", 10.0),
        ],
    )
    def test_reads_the_number_and_its_unit(self, text, expected):
        assert parse_zeiss_value(text) == pytest.approx(expected)

    def test_an_empty_field_has_no_number(self):
        assert parse_zeiss_value("Sample ID =") is None
        assert parse_zeiss_value("User Text = ") is None


# ── Refusing what SmartSEM does not write ────────────────────────────────────

class TestUnsupportedLayouts:
    def test_a_compressed_tiff_is_refused_rather_than_misread(self, tmp_path):
        path = write_tiff(tmp_path / "c.tif", micrograph(64, 64), zeiss_tag(BASE_FIELDS))
        raw = bytearray(path.read_bytes())
        # Rewrite the Compression tag's value in place: 5 is LZW.
        head = raw.index(struct.pack("<HHI", 259, 3, 1))
        raw[head + 8: head + 10] = struct.pack("<H", 5)
        path.write_bytes(bytes(raw))
        with pytest.raises(ValueError, match="compressed"):
            read_sem(path)

    def test_a_sixteen_bit_image_is_refused(self, tmp_path):
        path = write_tiff(tmp_path / "d.tif", micrograph(64, 64), zeiss_tag(BASE_FIELDS))
        raw = bytearray(path.read_bytes())
        head = raw.index(struct.pack("<HHI", 258, 3, 1))
        raw[head + 8: head + 10] = struct.pack("<H", 16)
        path.write_bytes(bytes(raw))
        with pytest.raises(ValueError, match="unsupported layout"):
            read_sem(path)

    def test_truncated_image_data_is_refused(self, tmp_path):
        path = write_tiff(tmp_path / "t.tif", micrograph(64, 64), zeiss_tag(BASE_FIELDS))
        raw = path.read_bytes()
        path.write_bytes(raw[: len(raw) - 500])
        with pytest.raises(ValueError):
            read_sem(path)


def test_sem_image_is_frozen(sem_file):
    image = read_sem(sem_file)
    assert isinstance(image, SemImage)
    with pytest.raises(Exception):
        image.pixel_size_nm = 1.0
