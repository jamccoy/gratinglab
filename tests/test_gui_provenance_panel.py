"""The collapsible provenance strip: its summary row (no Qt) and the widget.

``gratinglab.gui.provenance.summary`` picks the one row a collapsed panel
shows; ``gratinglab.gui.qt.provenance_panel.ProvenancePanel`` is the strip on
each solver tab.
"""

from __future__ import annotations

import pytest

from gratinglab.gui.provenance import Line, summary


class TestSummary:
    def test_a_clean_run_shows_its_first_row(self):
        lines = (
            Line("scalar 0.3.7  ·  2048 points\n"),
            Line("energy balance: ", "dim"),
            Line("ok\n", "ok"),
        )
        assert summary(lines) == (Line("scalar 0.3.7  ·  2048 points"),)

    def test_warnings_are_counted_in_the_warning_colour_only(self):
        lines = (
            Line("rcwa 0.3.7\n"),
            Line("  ⚠ first\n", "warn"),
            Line("  ⚠ second\n", "warn"),
        )
        head, count = summary(lines)
        assert head == Line("rcwa 0.3.7")
        assert count == Line("  ·  2 warnings", "warn")

    def test_one_warning_is_singular(self):
        assert summary((Line("x\n"), Line("w\n", "warn")))[-1].text.endswith("1 warning")

    def test_a_failure_is_shown_itself(self):
        lines = (Line("scalar\n"), Line("solve failed: no coating\n", "bad"))
        assert summary(lines) == (Line("solve failed: no coating", "bad"),)

    def test_a_row_takes_its_worst_run(self):
        lines = (Line("first\n"), Line("energy balance: ", "dim"), Line("NO\n", "bad"))
        assert summary(lines) == (Line("energy balance: NO", "bad"),)

    def test_nothing_gives_nothing(self):
        assert summary(()) == ()


pytest.importorskip("PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"')

from gratinglab.gui.qt.provenance_panel import ProvenancePanel  # noqa: E402
from gratinglab.gui.qt.theming import settings  # noqa: E402


@pytest.fixture
def panels(qtbot):
    made = []

    def make():
        panel = ProvenancePanel()
        qtbot.addWidget(panel)
        made.append(panel)
        return panel

    return make


@pytest.mark.gui
class TestPanel:
    def test_collapsed_by_default_with_the_text_still_there(self, panels):
        panel = panels()
        panel.show_lines((Line("scalar 0.3.7\n"), Line("more\n", "dim")))
        assert not panel.is_expanded()
        assert panel.text.isHidden()
        assert "more" in panel.text.toPlainText()
        assert "scalar 0.3.7" in panel.summary.text()

    def test_opening_one_opens_every_panel_and_is_remembered(self, panels):
        first, second = panels(), panels()
        first.toggle.click()
        assert first.is_expanded() and second.is_expanded()
        assert settings().value("panels/provenance_expanded", type=bool)
        assert panels().is_expanded(), "a panel built later starts open too"

    def test_a_failure_opens_the_panel_without_changing_the_choice(self, panels):
        panel = panels()
        panel.show_lines((Line("solve failed: no coating\n", "bad"),))
        assert panel.is_expanded()
        assert not settings().value("panels/provenance_expanded", False, type=bool)
