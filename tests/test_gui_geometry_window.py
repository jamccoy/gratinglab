"""The geometry window: a grating drawn from its own form, with no solver.

What it promises, and what each test holds it to: it draws without solving,
it keeps its last drawing through a half-typed field, and "Send to
efficiency window" copies the whole form into the one efficiency window
without solving there either. The form copy is checked field by field,
because a field `write_form` forgot would arrive as the efficiency window's
old value and look correct.
"""

import dataclasses

import pytest

pytest.importorskip(
    "PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"'
)

from gratinglab.gui.launcher import mode  # noqa: E402
from gratinglab.gui.qt import launcher as qt_launcher  # noqa: E402
from gratinglab.gui.state import FormState  # noqa: E402

#: Every field away from its default, so a field left unwritten shows.
#: The coating stays empty: the repository ships no optical constants.
CHANGED = FormState(
    period="400",
    profile_kind="Sinusoidal",
    blaze_angle="12",
    antiblaze_angle="60",
    depth_fraction="0.2",
    duty_cycle="0.4",
    profile_path="",
    mount="Classical",
    alpha="10",
    gamma="2.5",
    wavelength_start="2",
    wavelength_stop="6",
    wavelength_count="40",
    coating="",
    roughness="0.5",
    layers="",
)


@pytest.fixture(autouse=True)
def _no_leaked_windows():
    """Empty the launcher's registry around each test; see test_gui_launcher."""
    import shiboken6

    qt_launcher._OPEN_WINDOWS.clear()
    yield
    for window in list(qt_launcher._OPEN_WINDOWS):
        if shiboken6.isValid(window):
            window.close()
    qt_launcher._OPEN_WINDOWS.clear()


@pytest.fixture
def window(qtbot):
    window = qt_launcher.open_window(mode("geometry"))
    qtbot.addWidget(window)
    return window


def _period(window):
    return window.geometry_tab._parsed.problem.period


@pytest.mark.gui
class TestTheWindow:
    def test_it_is_the_geometry_mode_and_draws_at_once(self, window):
        from gratinglab.gui.qt.geometry_window import GeometryWindow

        assert isinstance(window, GeometryWindow)
        assert window._gratinglab_mode == "geometry"
        assert window.geometry_tab._parsed is not None, "nothing drawn on opening"

    def test_it_has_no_solver(self, window):
        assert not hasattr(window, "solve")
        assert not hasattr(window, "tabs")

    def test_enter_redraws(self, window):
        before = _period(window)
        field = window.geometry._fields["period"]
        field.setText("400")
        field.returnPressed.emit()
        assert _period(window) != before
        assert _period(window) == pytest.approx(before * 400 / 315.15)

    def test_a_half_typed_field_keeps_the_drawing(self, window):
        before = window.geometry_tab._parsed
        window.geometry._fields["period"].setText("-")
        window.refresh()
        assert window.geometry_tab._parsed is before

    def test_the_send_button_sits_under_the_form(self, window, qapp):
        """Just below the last field, not at the foot of a tall dock."""
        window.resize(1280, 820)
        window.show()
        qapp.processEvents()
        panel, button = window.geometry, window.send_button
        form_bottom = panel.mapTo(window, panel.rect().topLeft()).y() + panel.sizeHint().height()
        gap = button.mapTo(window, button.rect().topLeft()).y() - form_bottom
        assert 0 <= gap < 20, f"the button is {gap} px below the form"
        assert window.geometry_dock.height() - gap > 400, "the dock was too short to tell"

    def test_the_profile_plot_text_fits_a_narrow_window(self, window, qapp):
        """The groove profile's title and x label stay inside its figure at
        the dock's narrowest; one long title used to run off its right edge."""
        window.resize(1000, 700)
        window.show()
        qapp.processEvents()
        panel = window.geometry_tab.profile_panel
        figure, axes = panel._figure, panel._axes
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        assert figure.bbox.width < 300, "the pane was too wide to tell"
        title = figure.get_suptitle()
        assert "two periods" in title and "period" in title.splitlines()[1]
        for text in (figure._suptitle, axes.title, axes.xaxis.label):
            extent = text.get_window_extent(renderer)
            assert extent.x0 >= 0 and extent.x1 <= figure.bbox.width, text.get_text()

    def test_the_go_menu_reaches_the_efficiency_window(self, window):
        titles = [a.text() for a in window._gratinglab_switch_actions]
        assert "Compute efficiency" in titles


@pytest.mark.gui
class TestTheFormRoundTrip:
    def test_write_form_sets_every_field(self, window):
        window.geometry.write_form(CHANGED)
        assert window.geometry.read_form() == CHANGED

    def test_every_field_differs_from_its_default_but_coating_and_paths(self):
        """Otherwise the round trip above could pass by writing nothing."""
        default = FormState()
        same = {f.name for f in dataclasses.fields(FormState)
                if getattr(default, f.name) == getattr(CHANGED, f.name)}
        assert same <= {"coating", "profile_path", "layers"}


@pytest.mark.gui
class TestSending:
    def test_it_copies_the_form_and_solves_nothing(self, qtbot, window):
        efficiency = qt_launcher.open_window(mode("efficiency"))
        qtbot.addWidget(efficiency)
        with qtbot.waitSignal(efficiency.solved, timeout=60_000):
            pass  # the solve its constructor starts

        window.geometry.write_form(CHANGED)
        solves = []
        efficiency.solved.connect(lambda: solves.append(True))
        sent = window.send_to_efficiency()

        assert sent is efficiency, "the open efficiency window is reused"
        assert efficiency.geometry.read_form() == CHANGED
        qtbot.wait(300)
        assert not solves, "sending must not start a solve"
        assert len(qt_launcher._OPEN_WINDOWS) == 2

    def test_it_opens_the_efficiency_window_when_none_is_open(self, qtbot, window):
        window.geometry.write_form(CHANGED)
        sent = window.send_to_efficiency()
        qtbot.addWidget(sent)
        assert sent._gratinglab_mode == "efficiency"
        assert sent.geometry.read_form() == CHANGED

    def test_a_form_that_does_not_parse_is_not_sent(self, window, monkeypatch):
        from PySide6.QtWidgets import QMessageBox

        warned = []
        monkeypatch.setattr(QMessageBox, "warning",
                            lambda *args, **kwargs: warned.append(args))
        window.geometry._fields["period"].setText("-")

        assert window.send_to_efficiency() is None
        assert warned, "the refusal must say why"
        assert qt_launcher.window_for(mode("efficiency")) is None
