"""How a failed solve is worded, what the windows remember, and the File menu.

Covers ``gratinglab.gui.failures`` (no Qt), ``gratinglab.gui.qt.recent``, and
the File menu and layout memory of the efficiency window. The metrology
window's Ctrl+O and its Send button are tested beside the rest of that window
in ``tests/metrology/test_metrology_gui_qt.py``.
"""

from __future__ import annotations

import pytest

from gratinglab.gui.docs import REPOSITORY
from gratinglab.gui.failures import describe_failure, is_refusal
from gratinglab.profiles._types import ProfileRepresentationError
from gratinglab.solvers.base import UnsupportedConfiguration


class TestDescribeFailure:
    @pytest.mark.parametrize("exc", [
        UnsupportedConfiguration("rcwa cannot model this"),
        ProfileRepresentationError("no polygon for this profile"),
        ValueError("truncation must be positive"),
    ])
    def test_a_refusal_is_its_own_message_without_the_class_name(self, exc):
        assert is_refusal(exc)
        text = describe_failure(exc)
        assert str(exc) in text
        assert type(exc).__name__ not in text
        assert "bug" not in text

    def test_anything_else_keeps_the_class_and_asks_for_a_report(self):
        text = describe_failure(ZeroDivisionError("division by zero"))
        assert "ZeroDivisionError" in text
        assert f"{REPOSITORY}/issues" in text

    def test_an_empty_message_still_reads(self):
        assert describe_failure(RuntimeError()).endswith("/issues")
        assert "(no message)" in describe_failure(RuntimeError())


pytest.importorskip("PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"')

from PySide6.QtGui import QKeySequence  # noqa: E402

from gratinglab.gui.qt import recent  # noqa: E402


@pytest.mark.gui
class TestRememberedFolders:
    def test_a_folder_is_remembered_per_kind(self, tmp_path):
        scan = tmp_path / "scans" / "a.spm"
        scan.parent.mkdir()
        recent.remember("scans", scan)
        assert recent.start_in("scans") == str(scan.parent.resolve())
        assert recent.start_in("exports") == ""
        assert recent.start_in("scans", "b.csv") == str(scan.parent.resolve() / "b.csv")

    def test_a_folder_that_has_gone_falls_back(self, tmp_path):
        gone = tmp_path / "gone" / "a.json"
        gone.parent.mkdir()
        recent.remember("profiles", gone)
        gone.parent.rmdir()
        assert recent.start_in("profiles", fallback="elsewhere") == "elsewhere"
        assert recent.start_in("profiles", "x.json") == "x.json"

    def test_a_cancelled_dialog_records_nothing(self):
        recent.remember("tables", "")
        assert recent.start_in("tables") == ""


@pytest.fixture
def efficiency(qtbot):
    from gratinglab.gui.qt.main_window import MainWindow

    window = MainWindow()
    qtbot.addWidget(window)
    with qtbot.waitSignal(window.solved, timeout=60_000):
        pass  # the solve the constructor starts
    return window


@pytest.mark.gui
class TestEfficiencyWindow:
    def test_file_menu_shortcuts(self, efficiency):
        actions = efficiency._file_actions
        assert actions["open"].shortcut().matches(
            QKeySequence(QKeySequence.StandardKey.Open)
        ) == QKeySequence.SequenceMatch.ExactMatch
        assert actions["solve"].shortcut() == QKeySequence("Ctrl+Return")

    def test_solve_this_tab_solves(self, qtbot, efficiency):
        with qtbot.waitSignal(efficiency.solved, timeout=60_000):
            efficiency._file_actions["solve"].trigger()

    def test_size_is_remembered_across_windows(self, qtbot, efficiency):
        from gratinglab.gui.qt.main_window import MainWindow

        # Shown, because a geometry restored into a hidden window is not laid
        # out yet. Height only: the window's minimum width (about 1080 px) is
        # wider than the offscreen platform's 800 px screen, and its minimum height is about 550 px, so the width a
        # restore lands on depends on the screen, not on what was saved.
        efficiency.show()
        efficiency.resize(efficiency.width(), 640)
        qtbot.waitUntil(lambda: efficiency.height() == 640)
        efficiency.close()
        again = MainWindow()
        qtbot.addWidget(again)
        with qtbot.waitSignal(again.solved, timeout=60_000):
            pass
        again.show()
        qtbot.waitUntil(lambda: again.height() == 640)
