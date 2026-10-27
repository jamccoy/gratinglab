"""
Two-column line profiles: ``x z`` exports read as a one-line scan.

A Gwyddion profile export carries position and height per line, so the scan
width comes from the x column rather than from a header or a default. The
check that matters is the end-to-end one: the same surface exported as an
image and as a row-averaged profile must average to the same groove.
"""
import contextlib
import io

import numpy as np
import pytest

from scans import SYNTHETIC

import matplotlib
matplotlib.use('Agg')

from gratinglab.metrology.boundary import build_boundary_profile                # noqa: E402
from gratinglab.metrology.core.processing import load_afm_data                  # noqa: E402
from gratinglab.metrology.settings import AnalysisSettings                      # noqa: E402


def _quiet(fn, *args, **kwargs):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args, **kwargs)


def _write_profile(path, x, z, header="Horizontal profile\nx y\n[m] [m]\n"):
    with open(path, 'w', encoding='utf-8', newline='\n') as handle:
        handle.write(header)
        for a, b in zip(x, z):
            handle.write(f"{a:.9e} {b:.9e}\n")
    return str(path)


@pytest.fixture(scope='module')
def image():
    return _quiet(load_afm_data, str(SYNTHETIC))


@pytest.fixture
def profile_file(tmp_path, image):
    data, width_um = image
    z = np.nanmean(data, axis=0)
    x = np.linspace(0.0, width_um * 1e-6, len(z))
    return _write_profile(tmp_path / 'profile.txt', x, z)


def test_two_columns_load_as_one_scan_line(profile_file, image):
    data, width_um = _quiet(load_afm_data, profile_file)
    assert data.shape == (1, image[0].shape[1])
    assert width_um == pytest.approx(image[1])


def test_width_comes_from_the_x_column_not_a_default(tmp_path):
    x = np.linspace(0.0, 1.234e-6, 200)
    path = _write_profile(tmp_path / 'p.txt', x, np.zeros_like(x))
    _, width_um = _quiet(load_afm_data, path)
    assert width_um == pytest.approx(1.234)


def test_an_offset_x_column_still_gives_the_span(tmp_path):
    x = np.linspace(5e-6, 7e-6, 200)
    path = _write_profile(tmp_path / 'p.txt', x, np.zeros_like(x))
    _, width_um = _quiet(load_afm_data, path)
    assert width_um == pytest.approx(2.0)


def test_bracketed_units_are_honoured(tmp_path):
    x = np.linspace(0.0, 2000.0, 200)                  # nm
    z = 50.0 * np.sin(x / 100.0)                        # nm
    path = _write_profile(tmp_path / 'p.txt', x, z, header="x [nm] z [nm]\n")
    data, width_um = _quiet(load_afm_data, path)
    assert width_um == pytest.approx(2.0)
    assert np.ptp(data) == pytest.approx(100e-9, rel=1e-3)


def test_uneven_spacing_is_refused_rather_than_stretched(tmp_path):
    x = np.linspace(0.0, 1.0, 200) ** 2 * 2e-6
    path = _write_profile(tmp_path / 'p.txt', x, np.zeros_like(x))
    with pytest.raises(ValueError, match="evenly spaced"):
        _quiet(load_afm_data, path)


def test_image_files_still_take_the_image_path(image):
    data, width_um = image
    assert data.ndim == 2 and data.shape[0] > 1
    assert width_um == pytest.approx(2.0)


def test_profile_and_image_average_to_the_same_groove(profile_file, image):
    settings = AnalysisSettings()
    from_image = _quiet(build_boundary_profile, image[0], image[1], settings)
    data, width_um = _quiet(load_afm_data, profile_file)
    from_profile = _quiet(build_boundary_profile, data, width_um, settings)

    assert from_profile.n_used == from_image.n_used
    assert from_profile.period_nm == pytest.approx(from_image.period_nm, rel=1e-9)
    np.testing.assert_allclose(from_profile.y_norm, from_image.y_norm, atol=1e-6)
