"""The Materials tab's logic without a window, then the tab in the real one."""

import dataclasses

import pytest
from matplotlib.figure import Figure

from gratinglab import materials
from gratinglab.gui.materials_view import covers_scan, draw, lines, material_name, preview
from gratinglab.gui.theme import LIGHT


def text(p):
    return "".join(line.text for line in lines(p))


class TestTheName:
    @pytest.mark.parametrize(
        "kind, name, density, expect",
        [
            ("table", "Au", "", "Au"),
            ("table", "Au", "17.5", "Au@17.5"),
            ("factors", "SiO2", "2.2", "factors:SiO2@2.2"),
            ("factors", "Au", "", "factors:Au"),
            ("typical", " C ", " 2.0 ", "typical:C@2.0"),
            ("typical", "", "2.0", ""),
        ],
    )
    def test_built_from_the_fields(self, kind, name, density, expect):
        assert material_name(kind, name, density) == expect


class TestThePreview:
    def test_a_model_covering_the_scan(self):
        p = preview("typical:Au", "1.0", "5.0")
        assert p.near_unity and covers_scan(p) is True
        said = text(p)
        assert "typical-material model" in said and "Covers the dock's scan, 1 to 5 nm" in said

    def test_a_scan_outside_the_table_is_warned_about(self):
        p = preview("typical:Au", "30", "80")  # the model stops at 41.3 nm
        assert covers_scan(p) is False
        assert any(line.tag == "warn" and "Does not cover" in line.text for line in lines(p))

    def test_caveats_for_the_scan_are_passed_on(self):
        # Gold's M edges sit near 0.56 nm; the model says so for a scan across them.
        assert any("absorption edge" in line.text for line in lines(preview("typical:Au", "0.5", "0.6")))

    def test_an_unset_scan_is_not_checked(self):
        p = preview("typical:Au", "", "nonsense")
        assert covers_scan(p) is None and "not checked" in text(p)

    def test_an_unknown_name_raises_what_lookup_raises(self):
        with pytest.raises(materials.UnknownMaterial):
            preview("Unobtainium")

    def test_a_visible_light_table_plots_n_and_k(self, tmp_path, monkeypatch):
        monkeypatch.setenv(materials.LIBRARY_ENV, str(tmp_path))
        (tmp_path / "Glass.nk").write_text(
            "wavelength (nm), n, k\n400, 1.53, 0\n700, 1.51, 0\n", encoding="utf-8")
        p = preview("Glass", "450", "650")
        assert not p.near_unity
        axes = Figure().subplots()
        draw(axes, p, LIGHT)
        assert [line.get_label() for line in axes.get_lines()] == ["n", "k"]
        assert axes.get_yscale() == "linear"

    def test_an_x_ray_material_plots_decrement_and_absorption_on_a_log_axis(self):
        axes = Figure().subplots()
        draw(axes, preview("typical:SiO2", "1", "5"), LIGHT)
        assert [line.get_label() for line in axes.get_lines()] == ["decrement δ", "absorption β"]
        assert axes.get_yscale() == "log"


class TestTheTab:
    @pytest.fixture
    def win(self, qtbot):
        pytest.importorskip("PySide6")
        from gratinglab.gui.qt.main_window import MainWindow

        window = MainWindow()
        qtbot.addWidget(window)
        with qtbot.waitSignal(window.solved, timeout=60_000):
            pass
        return window

    def test_it_previews_as_the_fields_change(self, win):
        tab = win.materials_tab
        tab.set_fields("typical", "Au")
        assert tab.coating_name() == "typical:Au"
        assert "typical-material model" in tab._panel.toPlainText()
        assert tab._use.isEnabled()

    def test_a_bad_name_disables_the_buttons_and_says_why(self, win):
        tab = win.materials_tab
        tab.set_fields("typical", "Xx")
        assert not tab._use.isEnabled() and not tab._layer.isEnabled()
        assert "Xx" in tab._panel.toPlainText()

    def test_the_buttons_write_the_dock(self, win):
        tab = win.materials_tab
        tab.set_fields("typical", "Si")
        tab._use.click()
        assert win.geometry.read_form().coating == "typical:Si"
        tab.set_fields("typical", "C", "2.0")
        tab._layer.click()
        assert win.geometry.read_form().layers == "typical:C@2.0"

    def test_the_library_is_listed_and_a_click_previews_it(self, win):
        tab = win.materials_tab
        names = [tab._list.item(i).text() for i in range(tab._list.count())]
        assert names == list(materials.available()) and "Au" in names
        tab._list.itemClicked.emit(tab._list.item(names.index("Au")))
        assert tab.coating_name() == "Au" and "CXRO/Henke export" in tab._panel.toPlainText()
        assert "Scattering factors: none installed" in tab._elements.text()

    def test_it_follows_the_docks_scan(self, win):
        tab = win.materials_tab
        tab.set_fields("typical", "Au")
        win.geometry.write_form(dataclasses.replace(
            win.geometry.read_form(), wavelength_start="30", wavelength_stop="80"))
        assert "Does not cover" in tab._panel.toPlainText()
