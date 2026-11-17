"""The coatings window: its own stack, its two tabs, and Send in both directions."""

import dataclasses

import pytest

pytest.importorskip("PySide6")

from gratinglab.gui.launcher import mode  # noqa: E402
from gratinglab.gui.qt import launcher as qt_launcher  # noqa: E402
from gratinglab.gui.qt.coatings_window import CoatingsWindow  # noqa: E402
from gratinglab.gui.state import FormState  # noqa: E402

pytestmark = pytest.mark.gui


@pytest.fixture
def coat(qtbot):
    window = CoatingsWindow()
    qtbot.addWidget(window)
    return window


def stack(window, **fields):
    window.stack.write_form(dataclasses.replace(window.stack.read_form(), **fields))


class TestTheWindow:
    def test_it_holds_materials_and_the_mirror_and_no_solver(self, coat):
        assert [coat.tabs.tabText(i) for i in range(coat.tabs.count())] == ["Materials", "Mirror"]
        assert not hasattr(coat, "solve")

    def test_its_dock_shows_the_stack_and_the_wavelengths_alone(self, coat):
        coat.show()
        shown = {key: coat.stack._fields[key].isVisible()
                 for key in ("period", "alpha", "roughness", "wavelength_start")}
        assert shown == {"period": False, "alpha": False, "roughness": True, "wavelength_start": True}
        # The hidden groups still give a whole form, at its defaults.
        assert coat.stack.read_form().period == FormState().period

    def test_the_base_material_is_called_the_substrate(self, coat):
        from PySide6.QtWidgets import QLabel

        labels = {label.text() for label in coat.stack.findChildren(QLabel)}
        assert "Substrate" in labels and "Grating" not in labels


class TestTheMirrorTab:
    def test_it_reads_the_dock_and_answers(self, coat):
        stack(coat, coating="typical:Si", layers="typical:Au 10")
        tab = coat.mirror_tab
        tab._fields["graze_count"].setText("20")
        tab.calculate()
        assert "typical:Au at 10 nm" in tab._panel.toPlainText()
        assert tab._depth_axes.get_lines()

    def test_it_waits_until_it_is_showing(self, coat, qtbot):
        tab = coat.mirror_tab
        coat.stack.set_grating_material("typical:Au")
        tab.calculate()
        before = tab._result
        assert before is not None
        coat.stack._fields["roughness"].setText("0.4")
        coat.stack.changed.emit()
        assert tab._stale and tab._result is before  # hidden: nothing computed
        coat.tabs.setCurrentWidget(tab)
        coat.show()
        qtbot.waitUntil(lambda: not tab._stale, timeout=10_000)
        assert tab._result.mirror.roughness == 0.4

    def test_it_hands_a_bulk_layer_to_the_dock_when_asked(self, coat):
        tab = coat.mirror_tab
        stack(coat, coating="typical:Si", layers="typical:Au 10")
        tab._fields["graze_count"].setText("20")
        tab.calculate()
        assert not tab._hand_off.isEnabled()  # 10 nm is not bulk past grazing
        stack(coat, layers="typical:Au 400")
        tab.calculate()
        assert tab._hand_off.isEnabled() and "reads as bulk" in tab._verdict.text()
        asked = []
        tab._confirm_hand_off = lambda *what: asked.append(what) or False
        tab._hand_off.click()
        assert asked == [("typical:Au", 1, "typical:Si")]
        assert coat.stack.read_form().layers == "typical:Au 400"  # declined: untouched
        tab._confirm_hand_off = lambda *what: True
        tab._hand_off.click()
        form = coat.stack.read_form()
        assert form.coating == "typical:Au" and form.layers == ""

    def test_bad_fields_are_reported_in_its_panel(self, coat):
        tab = coat.mirror_tab
        tab._fields["tolerance"].setText("-1")
        tab.calculate()
        assert "tolerance" in tab._panel.toPlainText()

    def test_marked_angles_are_drawn(self, coat):
        coat.stack.set_grating_material("typical:Au")
        tab = coat.mirror_tab
        tab._fields["graze_count"].setText("20")
        tab.set_marks("facet 1.5")
        tab.calculate()
        assert [t.get_text() for t in tab._r_axes.texts] == ["facet 1.5°"]


class TestSendingBetweenWindows:
    @pytest.fixture
    def both(self, qtbot, monkeypatch):
        """The efficiency and coatings windows as the launcher opens them."""
        monkeypatch.setattr(qt_launcher, "_OPEN_WINDOWS", [])
        efficiency = qt_launcher.open_window(mode("efficiency"))
        qtbot.addWidget(efficiency)
        with qtbot.waitSignal(efficiency.solved, timeout=60_000, raising=False):
            pass
        return efficiency

    def test_the_efficiency_window_sends_its_stack_and_its_mounts_angles(self, both, qtbot):
        efficiency = both
        efficiency.geometry.write_form(dataclasses.replace(
            efficiency.geometry.read_form(), coating="typical:Si",
            layers="typical:Ni 5, typical:Au 20", roughness="0.3",
            wavelength_start="2", wavelength_stop="4",
        ))
        coat = efficiency.send_to_coatings()
        qtbot.addWidget(coat)
        assert isinstance(coat, CoatingsWindow)
        form = coat.stack.read_form()
        assert (form.coating, form.layers, form.roughness) == ("typical:Si", "typical:Ni 5, typical:Au 20", "0.3")
        assert (form.wavelength_start, form.wavelength_stop) == ("2", "4")
        marks = coat.mirror_tab._fields["marks"].text()
        assert marks.startswith("grating plane 1.359") and "blaze facet 1.495" in marks
        assert efficiency.send_to_coatings() is coat  # the one window, not a second

    def test_the_coatings_window_sends_only_the_surface_back(self, both, qtbot):
        efficiency = both
        before = efficiency.geometry.read_form()
        coat = efficiency.send_to_coatings()
        qtbot.addWidget(coat)
        stack(coat, coating="typical:Au", layers="typical:C 3", roughness="0.5",
              wavelength_start="7", wavelength_stop="9")
        assert coat.send_to_efficiency() is efficiency
        after = efficiency.geometry.read_form()
        assert (after.coating, after.layers, after.roughness) == ("typical:Au", "typical:C 3", "0.5")
        # Its grating, mount and wavelengths are its own.
        assert (after.period, after.alpha, after.wavelength_start) == (
            before.period, before.alpha, before.wavelength_start)

    def test_a_library_change_reaches_every_window(self, both, qtbot, tmp_path, monkeypatch):
        from gratinglab import materials

        efficiency = both
        coat = efficiency.send_to_coatings()
        qtbot.addWidget(coat)
        library = tmp_path / "library"
        library.mkdir()
        monkeypatch.setenv(materials.LIBRARY_ENV, str(library))
        (library / "Glass.nk").write_text(
            "wavelength (nm), n, k\n400, 1.53, 0\n700, 1.51, 0\n", encoding="utf-8")
        coat.materials_tab.library.rescan()
        combo = efficiency.geometry._fields["coating"]
        assert "Glass" in [combo.itemData(i) for i in range(combo.count())]
        tab = coat.materials_tab
        assert "Glass" in [tab._list.item(i).text() for i in range(tab._list.count())]
        assert "Glass" in efficiency._setup_tab._browser.toPlainText()
