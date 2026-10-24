"""
Reading Asylum Research .ibw files.

No file from an instrument was available when the reader was written
(2026-10-24), so every file here is built to WaveMetrics Technical Note PTN003
by ``_write_ibw`` below. That tests the reader against the specification as
this file's author read it, and the module docstring of ``io/ibw.py`` says so.
What it can still establish is the part that silently corrupts results when
wrong: orientation (lines against points), the X scale, byte order, channel
choice, and that a height comes through as a height. The end-to-end test does
that physically, by measuring a known blaze angle and period through the full
pipeline. The first real Asylum scan should be added beside these as a fixture,
with an export of the same scan, as the .spm reader's was.
"""
import contextlib
import io
import struct

import numpy as np
import pytest

from gratinglab.metrology.core.processing import load_afm_data
from gratinglab.metrology.io.ibw import (
    is_igor_wave_file, list_channels, read_ibw)
from gratinglab.metrology.io.spm import is_nanoscope_file


def _write_ibw(path, layers, *, dx_m, order="<", units="m", note="",
               version=5, kind=2):
    """An Igor version 5 binary wave holding ``layers`` (name -> 2-D array).

    Each array is (lines, points). Igor stores the first dimension fastest,
    so the cube is written as (layers, lines, points) in C order, which puts
    points first in memory, as an Asylum image does.
    """
    names = list(layers)
    cube = np.stack([np.asarray(layers[n]) for n in names])
    n_layers, rows, cols = cube.shape
    dtype = {2: "f4", 4: "f8", 0x10: "i2"}[kind]
    data = cube.astype(order + dtype).tobytes()
    npnts = cube.size

    labels = [b"", b"", b"", b""]
    labels[2] = b"".join(
        name.encode().ljust(32, b"\0") for name in [""] + names)
    note_bytes = note.encode("latin-1")
    wfm_size = 320 + len(data)

    bin_header = struct.pack(
        order + "hhiiii4i4iiii", version, 0, wfm_size, 0, len(note_bytes), 0,
        0, 0, 0, 0, *(len(b) for b in labels), 0, 0, 0)
    wave = bytearray(320)
    struct.pack_into(order + "i", wave, 12, npnts)
    struct.pack_into(order + "h", wave, 16, kind)
    struct.pack_into(order + "4i", wave, 68, cols, rows, n_layers, 0)
    struct.pack_into(order + "4d", wave, 84, dx_m, dx_m, 1.0, 1.0)
    wave[148:152] = units.encode().ljust(4, b"\0")[:4]
    path.write_bytes(bin_header + bytes(wave) + data + note_bytes
                     + b"".join(labels))
    return path


def _ramp(rows=8, cols=16):
    """Distinct along each axis, so a transposed read cannot pass."""
    y, x = np.mgrid[0:rows, 0:cols]
    return (x * 1e-9 + y * 1e-7).astype(float)


def _sawtooth(period_nm, blaze_deg, antiblaze_deg, scan_um, cols, rows, seed=0):
    """A blazed grating, heights in metres, grooves running along the lines."""
    x = np.linspace(0.0, scan_um * 1e3, cols, endpoint=False)
    t = np.mod(x, period_nm)
    tan_b, tan_a = np.tan(np.radians(blaze_deg)), np.tan(np.radians(antiblaze_deg))
    apex = period_nm * tan_a / (tan_a + tan_b)
    profile = np.where(t < apex, t * tan_b, (period_nm - t) * tan_a)
    rng = np.random.default_rng(seed)
    image = np.tile(profile, (rows, 1)) + rng.normal(0.0, 0.2, (rows, cols))
    return image * 1e-9


class TestFormat:
    def test_recognised_by_content_not_extension(self, tmp_path):
        path = _write_ibw(tmp_path / "scan.dat", {"HeightRetrace": _ramp()},
                          dx_m=1e-8)
        assert is_igor_wave_file(path)
        text = tmp_path / "export.txt"
        text.write_text("# Width: 2 um\n0 1 2\n", encoding="utf-8")
        assert not is_igor_wave_file(text)

    def test_a_nanoscope_file_is_not_mistaken_for_one(self, tmp_path):
        spm = tmp_path / "scan.spm"
        spm.write_bytes(b"\\*File list\r\n" + bytes(600))
        assert is_nanoscope_file(spm)
        assert not is_igor_wave_file(spm)

    @pytest.mark.parametrize("order", ["<", ">"])
    def test_either_byte_order(self, tmp_path, order):
        image = _ramp()
        path = _write_ibw(tmp_path / "s.ibw", {"HeightRetrace": image},
                          dx_m=1e-8, order=order)
        data, _ = read_ibw(path, channel="Height")
        np.testing.assert_array_equal(data, image.astype(np.float32))

    def test_lines_are_rows_and_points_are_columns(self, tmp_path):
        image = _ramp(rows=8, cols=16)
        path = _write_ibw(tmp_path / "s.ibw", {"HeightRetrace": image}, dx_m=1e-8)
        data, _ = read_ibw(path, channel="Height")
        assert data.shape == (8, 16)
        np.testing.assert_allclose(data[0, :3], [0.0, 1e-9, 2e-9], atol=1e-15)
        np.testing.assert_allclose(data[:3, 0], [0.0, 1e-7, 2e-7], rtol=1e-6)

    def test_scan_width_from_the_x_scaling(self, tmp_path):
        path = _write_ibw(tmp_path / "s.ibw", {"HeightRetrace": _ramp(cols=200)},
                          dx_m=10e-9)
        assert read_ibw(path, channel="Height")[1] == pytest.approx(2.0)

    def test_scan_width_from_the_note_without_scaling(self, tmp_path):
        path = _write_ibw(tmp_path / "s.ibw", {"HeightRetrace": _ramp()},
                          dx_m=1.0, note="ScanSize: 5e-06\rScanRate: 1\r")
        assert read_ibw(path, channel="Height")[1] == pytest.approx(5.0)

    def test_older_versions_are_refused_by_name(self, tmp_path):
        path = _write_ibw(tmp_path / "s.ibw", {"HeightRetrace": _ramp()},
                          dx_m=1e-8, version=2)
        assert is_igor_wave_file(path)
        with pytest.raises(ValueError, match="version 2"):
            read_ibw(path)

    def test_truncated_file_is_refused(self, tmp_path):
        path = _write_ibw(tmp_path / "s.ibw", {"HeightRetrace": _ramp()}, dx_m=1e-8)
        path.write_bytes(path.read_bytes()[:500])
        with pytest.raises(ValueError):
            read_ibw(path)


class TestChannels:
    LAYERS = ("HeightTrace", "HeightRetrace", "ZSensorTrace", "ZSensorRetrace",
              "AmplitudeRetrace", "PhaseRetrace")

    def _file(self, tmp_path, names=LAYERS):
        layers = {}
        for i, name in enumerate(names):
            if name.startswith("Phase"):
                layers[name] = np.full((4, 6), 37.0)  # degrees, not metres
            else:
                layers[name] = np.full((4, 6), (i + 1) * 1e-9)
        return _write_ibw(tmp_path / "s.ibw", layers, dx_m=1e-8)

    def test_every_layer_is_listed_with_its_direction(self, tmp_path):
        found = list_channels(self._file(tmp_path))
        assert [c.label for c in found] == list(self.LAYERS)
        assert {c.direction for c in found} == {"Trace", "Retrace"}

    def test_the_sensor_channel_is_preferred(self, tmp_path):
        data, _ = read_ibw(self._file(tmp_path))
        assert data[0, 0] == pytest.approx(4e-9)  # ZSensorRetrace

    def test_height_when_there_is_no_sensor(self, tmp_path):
        data, _ = read_ibw(self._file(tmp_path, ("HeightTrace", "HeightRetrace")))
        assert data[0, 0] == pytest.approx(2e-9)

    def test_direction_selects_the_layer(self, tmp_path):
        data, _ = read_ibw(self._file(tmp_path), direction="Trace")
        assert data[0, 0] == pytest.approx(3e-9)  # ZSensorTrace

    def test_a_missing_channel_names_what_the_file_has(self, tmp_path):
        with pytest.raises(ValueError, match="AmplitudeRetrace"):
            read_ibw(self._file(tmp_path), channel="Deflection")

    def test_a_channel_not_in_metres_is_refused(self, tmp_path):
        with pytest.raises(ValueError, match="not a height in metres"):
            read_ibw(self._file(tmp_path), channel="Phase")

    def test_declared_units_other_than_metres_are_refused(self, tmp_path):
        path = _write_ibw(tmp_path / "s.ibw", {"HeightRetrace": _ramp()},
                          dx_m=1e-8, units="nm")
        with pytest.raises(ValueError, match="'nm'"):
            read_ibw(path, channel="Height")


class TestPipeline:
    def test_load_afm_data_dispatches_on_content(self, tmp_path):
        image = _ramp(rows=8, cols=200)
        path = _write_ibw(tmp_path / "s.ibw", {"HeightRetrace": image}, dx_m=10e-9)
        with contextlib.redirect_stdout(io.StringIO()):
            data, scan = load_afm_data(str(path), default_scan_size=9.0)
        assert data.shape == (8, 200)
        assert scan == pytest.approx(2.0)

    def test_a_known_blaze_comes_back_through_the_analysis(self, tmp_path):
        """The physical check: a transposed or rescaled read changes the
        measured angle or period, and nothing else in this file would show it."""
        from gratinglab.metrology.analyzer import analyze_single_file
        from gratinglab.metrology.settings import AnalysisSettings

        image = _sawtooth(period_nm=400.0, blaze_deg=25.0, antiblaze_deg=70.0,
                          scan_um=4.0, cols=1024, rows=64)
        path = _write_ibw(tmp_path / "grating.ibw",
                          {"ZSensorRetrace": image, "ZSensorTrace": image},
                          dx_m=4e-6 / 1024)
        with contextlib.redirect_stdout(io.StringIO()):
            result = analyze_single_file(
                str(path), show_plots=False, settings=AnalysisSettings.from_config())
        # Measured 24.81 deg and 400.29 nm (2026-10-24). A transposed read
        # finds no grooves; a scale error moves the period by its factor.
        assert result["mean_angle"] == pytest.approx(25.0, abs=0.5)
        assert float(result["period_nm"]) == pytest.approx(400.0, abs=2.0)
