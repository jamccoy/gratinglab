"""The Mirror tab's logic, without a window: what it computes and says."""

import numpy as np
import pytest
from matplotlib.figure import Figure

from gratinglab.gui.mirror_view import (
    MirrorSettings,
    compute,
    draw,
    is_bulk_throughout,
    lines,
    mount_markers,
    verdict,
)
from gratinglab.gui.state import FormErrors, FormState
from gratinglab.gui.theme import LIGHT
from gratinglab.mirror import Mirror

AU_ON_SI = FormState(coating="typical:Si", layers="typical:Au 10")


def text(result):
    return "".join(line.text for line in lines(result))


class TestCompute:
    def test_an_angle_scan_is_the_mirror(self):
        result = compute(AU_ON_SI, MirrorSettings(graze_count="40"))
        assert result.x[0] == pytest.approx(0.1) and result.x[-1] == pytest.approx(6.0)
        expect = Mirror("typical:Si", result.mirror.layers).reflectivity(result.x, 1.0)
        assert np.allclose(result.reflectivity, expect)
        assert list(result.depths) == ["typical:Au", "typical:Si"]
        assert result.bulk.shape == result.x.shape

    def test_a_wavelength_scan_uses_the_docks_range(self):
        form = FormState(coating="typical:Au", wavelength_start="1", wavelength_stop="3", wavelength_count="25")
        result = compute(form, MirrorSettings(scan="wavelength", graze="1.5"))
        assert result.x[0] == 1.0 and result.x[-1] == 3.0 and len(result.x) == 25
        assert result.bulk is None  # no layers, nothing to test

    def test_every_bad_field_is_reported_at_once(self):
        with pytest.raises(FormErrors) as caught:
            compute(FormState(coating=""), MirrorSettings(graze_stop="0.05", tolerance="0"))
        fields = {e.field for e in caught.value.errors}
        assert {"coating", "tolerance", "graze_stop"} <= fields

    def test_layers_take_model_names(self):
        result = compute(FormState(coating="typical:Si", layers="typical:Ni 3, typical:C 2"), MirrorSettings(graze_count="5"))
        assert [layer.material for layer in result.mirror.layers] == ["typical:Ni", "typical:C"]


class TestWhatItSays:
    def test_a_thick_layer_is_bulk_throughout(self):
        result = compute(FormState(coating="typical:Si", layers="typical:Au 200"),
                         MirrorSettings(graze_start="0.2", graze_stop="3", graze_count="20"))
        assert "reads as bulk across the scan" in text(result)
        assert verdict(result).tag == "ok" and is_bulk_throughout(result)

    def test_a_thin_layer_is_bulk_only_at_grazing(self):
        result = compute(AU_ON_SI, MirrorSettings(graze_count="60"))
        said = text(result)
        assert "is bulk at" in said and "model the stack there" in said
        assert verdict(result).tag == "warn" and not is_bulk_throughout(result)
        # 10 nm of Au at 1 nm is bulk at the smallest angles only.
        assert result.bulk[0] <= 10.0 < result.bulk[-1]

    def test_a_very_thin_layer_is_never_bulk(self):
        result = compute(FormState(coating="typical:Si", layers="typical:Au 0.5"),
                         MirrorSettings(graze_start="0.5", graze_stop="5", graze_count="10"))
        assert "is not bulk anywhere" in text(result)

    def test_a_bare_substrate_has_nothing_to_hand_off(self):
        result = compute(FormState(coating="typical:Au"), MirrorSettings(graze_count="5"))
        assert "no layer to test" in verdict(result).text and not is_bulk_throughout(result)

    def test_the_stack_is_written_surface_first(self):
        result = compute(FormState(coating="typical:Si", layers="typical:Ni 3, typical:C 2"), MirrorSettings(graze_count="5"))
        assert "vacuum / typical:C 2 nm / typical:Ni 3 nm / typical:Si (substrate)" in text(result)


class TestTheScanAndTheMount:
    def test_a_wavelength_scan_can_use_the_tabs_own_range(self):
        form = FormState(coating="typical:Au", wavelength_start="1", wavelength_stop="3")
        own = MirrorSettings(scan="wavelength", dock_wavelengths=False,
                             wavelength_start="2", wavelength_stop="4", wavelength_count="9")
        result = compute(form, own)
        assert result.x[0] == 2.0 and result.x[-1] == 4.0 and len(result.x) == 9
        assert result.markers == ()  # angles mean nothing on a wavelength axis

    def test_the_mount_is_marked_on_an_angle_scan(self):
        # Off-plane at gamma = 1.5, alpha = 25 on a 29.5 degree blaze:
        # sin(graze) = sin(gamma) cos(delta - alpha), with delta = 0 for the plane.
        markers = dict(compute(AU_ON_SI, MirrorSettings(graze_count="5")).markers)
        gamma, alpha, blaze = np.radians([1.5, 25.0, 29.5])
        assert markers["grating plane"] == pytest.approx(np.degrees(np.arcsin(np.sin(gamma) * np.cos(alpha))))
        assert markers["blaze facet"] == pytest.approx(np.degrees(np.arcsin(np.sin(gamma) * np.cos(blaze - alpha))))

    def test_a_profile_without_a_blaze_has_only_the_plane(self):
        form = FormState(coating="typical:Au", profile_kind="Sinusoidal")
        assert [label for label, _ in mount_markers(form)] == ["grating plane"]

    def test_a_form_that_is_not_a_grating_has_no_markers(self):
        assert mount_markers(FormState(period="nonsense")) == ()

    def test_markers_inside_the_scan_are_drawn(self):
        figure = Figure()
        r_axes, depth_axes = figure.subplots(2, 1)
        draw(r_axes, depth_axes, compute(AU_ON_SI, MirrorSettings(graze_count="20")), LIGHT)
        assert sorted(t.get_text().split()[0] for t in r_axes.texts) == ["blaze", "grating"]


def test_it_draws():
    figure = Figure()
    r_axes, depth_axes = figure.subplots(2, 1)
    draw(r_axes, depth_axes, compute(AU_ON_SI, MirrorSettings(graze_count="20")), LIGHT)
    labels = [line.get_label() for line in depth_axes.get_lines()]
    assert "1/e depth in typical:Au" in labels and "typical:Au needed to read as bulk" in labels
    assert r_axes.get_yscale() == "log"
