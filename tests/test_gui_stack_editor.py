"""The coating-stack editor: rows outermost first, the form's text innermost first."""

import pytest

pytest.importorskip("PySide6")

from gratinglab.gui.qt.geometry_panel import GeometryPanel  # noqa: E402
from gratinglab.gui.qt.stack_editor import StackEditor, material_color  # noqa: E402
from gratinglab.gui.state import FormState, build  # noqa: E402


@pytest.fixture
def editor(qtbot):
    widget = StackEditor()
    qtbot.addWidget(widget)
    return widget


def shown(editor):
    """The rows top to bottom, as (material, thickness)."""
    return [(row.material.currentText(), row.thickness.text()) for row in editor._rows]


class TestTheOrder:
    def test_text_is_innermost_first_and_rows_are_outermost_first(self, editor):
        editor.setText("Ni 5, Au 20, C 2")
        assert shown(editor) == [("C", "2"), ("Au", "20"), ("Ni", "5")]
        assert editor.text() == "Ni 5, Au 20, C 2"

    def test_a_new_layer_goes_on_top(self, editor, qtbot):
        editor.setText("Au 20")
        with qtbot.waitSignal(editor.changed):
            editor.add_layer("C", "2")
        assert editor.text() == "Au 20, C 2"

    def test_moving_and_removing(self, editor, qtbot):
        editor.setText("Ni 5, Au 20, C 2")
        with qtbot.waitSignal(editor.changed):
            editor._rows[0].buttons["down"].click()  # C below Au
        assert editor.text() == "Ni 5, C 2, Au 20"
        with qtbot.waitSignal(editor.changed):
            editor._rows[1].buttons["remove"].click()
        assert editor.text() == "Ni 5, Au 20"

    def test_the_ends_cannot_move_past_the_ends(self, editor):
        editor.setText("Ni 5, Au 20")
        top, bottom = editor._rows
        assert not top.buttons["up"].isEnabled() and top.buttons["down"].isEnabled()
        assert bottom.buttons["up"].isEnabled() and not bottom.buttons["down"].isEnabled()


class TestWhatReachesTheForm:
    def test_empty_is_no_layers(self, editor):
        assert editor.text() == ""
        editor.setText("")
        assert editor.layer_count() == 0

    def test_a_half_filled_row_goes_through_for_the_form_to_refuse(self, editor):
        editor.add_layer("Au", "")
        assert editor.text() == "Au"  # state.build names it as incomplete
        editor.add_layer("", "")
        assert editor.text() == "Au"  # an untouched row is not a layer

    def test_setting_text_does_not_announce_a_change(self, editor, qtbot):
        with qtbot.assertNotEmitted(editor.changed):
            editor.setText("Au 20")

    def test_enter_in_a_thickness_asks_for_a_solve(self, editor, qtbot):
        editor.setText("Au 20")
        with qtbot.waitSignal(editor.solve_requested):
            editor._rows[0].thickness.returnPressed.emit()

    def test_library_names_fill_the_dropdowns_and_typing_survives(self, editor):
        editor.setText("typical:Au 20")
        editor.set_materials(("Ni", "Pt"))
        combo = editor._rows[0].material
        assert [combo.itemText(i) for i in range(combo.count())] == ["Ni", "Pt"]
        assert combo.currentText() == "typical:Au"


def test_a_material_keeps_its_colour():
    assert material_color("Au") == material_color("Au")
    assert material_color("typical:C") != material_color("typical:Au")


class TestInTheDock:
    @pytest.fixture
    def dock(self, qtbot):
        panel = GeometryPanel()
        qtbot.addWidget(panel)
        return panel

    def test_the_form_round_trips_through_the_rows(self, dock):
        form = FormState(coating="typical:Si", layers="typical:Ni 5, typical:Au 20")
        dock.write_form(form)
        assert dock.read_form() == form
        assert build(dock.read_form()).problem.layers[1].material == "typical:Au"

    def test_editing_a_row_changes_the_form(self, dock, qtbot):
        dock.write_form(FormState(layers="typical:Au 20"))
        row = dock._fields["layers"]._rows[0]
        with qtbot.waitSignal(dock.changed):
            row.thickness.setText("35")
            row.thickness.textEdited.emit("35")
        assert dock.read_form().layers == "typical:Au 35"
