"""The Qt widget layer.

Thin on purpose, like its Tk predecessor was: the logic worth testing lives in
``gui/state.py``, ``gui/provenance.py``, ``gui/richtext.py`` and
``gui/orders.py``, all of which are pure and tested headlessly. What is checked
here is the wiring -- and above all that the window is not a second source of
truth: whatever it plots must equal what the core produces for the same inputs.

Unlike the Tk suite this replaces, these tests need no display. ``conftest.py``
sets ``QT_QPA_PLATFORM=offscreen``, under which a window always constructs, so
they run on every CI job rather than only on macOS.
"""

import numpy as np
import pytest

pytest.importorskip(
    "PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"'
)

from gratinglab.gui.scalar_options import ScalarOptionsState, build_options  # noqa: E402
from gratinglab.gui.state import FormState, build  # noqa: E402
from gratinglab.solvers import scalar  # noqa: E402

#: Generous, because it covers a whole solve plus a redraw on a loaded CI
#: runner -- but finite, because a hung worker should fail rather than hang the
#: suite. The default scalar solve takes about 70 ms.
SOLVE_TIMEOUT_MS = 10_000


@pytest.fixture
def win(qtbot):
    """A constructed window with its first solve already drawn."""
    from gratinglab.gui.qt.main_window import MainWindow

    window = MainWindow()
    qtbot.addWidget(window)  # closes and deletes it at teardown
    with qtbot.waitSignal(window.solved, timeout=SOLVE_TIMEOUT_MS):
        pass  # the solve the constructor started
    return window


def resolve(qtbot, window):
    """Press Solve and wait for the result to land."""
    with qtbot.waitSignal(window.solved, timeout=SOLVE_TIMEOUT_MS):
        window.solve()


class TestToolchain:
    """That embedding matplotlib in Qt works at all.

    Deliberately the first thing built. Everything downstream assumes a
    ``FigureCanvasQTAgg`` can be parented into a Qt window on a headless
    machine, and discovering otherwise five milestones later would be
    expensive. The toolbar is here because it is the one piece of matplotlib's
    Qt backend that reaches for icons and a style, which is where a missing
    system library shows up first.
    """

    def test_a_canvas_and_toolbar_embed_in_a_window(self, qtbot):
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg,
            NavigationToolbar2QT,
        )
        from matplotlib.figure import Figure
        from PySide6.QtWidgets import QMainWindow, QVBoxLayout, QWidget

        window = QMainWindow()
        qtbot.addWidget(window)

        figure = Figure(figsize=(4, 3), layout="constrained")
        figure.add_subplot().plot([0, 1], [0, 1])
        canvas = FigureCanvasQTAgg(figure)
        toolbar = NavigationToolbar2QT(canvas, window)

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.addWidget(toolbar)
        layout.addWidget(canvas)
        window.setCentralWidget(central)

        canvas.draw()
        assert canvas.width() > 0
        assert toolbar.actions()

    def test_embedding_does_not_hijack_the_global_backend(self, qtbot):
        """Importing the Qt backend must not switch matplotlib out from under
        the rest of the suite.

        The Tk implementation called ``matplotlib.use("TkAgg")`` from inside a
        widget method, which mutated global state for every test that ran
        afterwards. Explicit canvas embedding needs no such call, and this
        pins that it stays true.
        """
        import matplotlib

        before = matplotlib.get_backend()
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure

        canvas = FigureCanvasQTAgg(Figure())
        canvas.draw()
        assert matplotlib.get_backend() == before


class TestNotASecondSourceOfTruth:
    """The single most important property of this layer.

    It survived the toolkit change unchanged, which is the point of having
    written it this way: if the window ever computed anything itself, these
    would drift.
    """

    def test_plotted_values_equal_a_direct_solve(self, win):
        parsed = build(FormState())
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(),
        )
        expected = scalar.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )
        assert np.array_equal(win._scan.efficiency, expected.efficiency)
        assert np.array_equal(win._scan.orders, expected.orders)

    def test_changing_a_field_changes_the_result_consistently(self, qtbot, win):
        win.geometry._fields["blaze_angle"].setText("12.0")
        resolve(qtbot, win)

        parsed = build(FormState(blaze_angle="12.0"))
        options = build_options(
            parsed.problem, parsed.illumination, parsed.wavelengths,
            ScalarOptionsState(),
        )
        expected = scalar.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **options
        )
        assert np.array_equal(win._scan.efficiency, expected.efficiency)


class TestBehaviour:
    def test_solves_on_construction(self, win):
        """Opening the app shows something, not an empty window."""
        assert win._scan is not None
        assert len(win._scan) == 200

    def test_invalid_input_reports_errors_without_crashing(self, win):
        win.geometry._fields["period"].setText("not a number")
        win.solve()  # rejected synchronously; never reaches the worker
        text = win._provenance.toPlainText()
        assert "period" in text
        assert "need attention" in text

    def test_a_failed_solve_leaves_the_previous_result_intact(self, win):
        good = win._scan
        win.geometry._fields["period"].setText("")
        win.solve()
        assert win._scan is good

    def test_the_panel_shows_what_provenance_decided(self, win):
        """All this layer owes: paint the lines it was given. What each line
        says is checked in tests/test_gui_provenance.py."""
        from gratinglab.gui.provenance import provenance_lines

        expected = provenance_lines(
            win._scan, win._energy, win._parsed.lambda_over_period
        )
        shown = win._provenance.toPlainText()
        for line in expected:
            assert line.text.strip() in shown

    def test_reads_every_form_field(self, win):
        form = win.geometry.read_form()
        assert isinstance(form, FormState)
        assert form.period == "315.15"

    def test_a_geometry_field_missing_from_formstate_stops_the_window_opening(
        self, qtbot
    ):
        """`GeometryPanel._fields`' keys must exactly match FormState's field
        names -- `read_form` builds one by keyword. Previously a mismatch was
        a TypeError at solve time, visible only to whoever pressed Solve."""
        from gratinglab.gui.qt import geometry_panel as module

        original = module.GeometryPanel._build

        def drop_a_field(self):
            original(self)
            del self._fields["period"]

        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(module.GeometryPanel, "_build", drop_a_field)
            with pytest.raises(AssertionError, match="period"):
                module.GeometryPanel()

    def test_a_scalar_option_field_missing_stops_the_window_opening(self, qtbot):
        """The other half of the same guard, now that scalar's own options
        are checked separately from geometry's."""
        from gratinglab.gui.qt import main_window as module

        original = module.MainWindow._build_solver_controls

        def drop_a_field(self):
            widget = original(self)
            del self._fields["quadrature_points"]
            self._check_fields_match_scalar_options()
            return widget

        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(module.MainWindow, "_build_solver_controls", drop_a_field)
            with pytest.raises(AssertionError, match="quadrature_points"):
                module.MainWindow()

    def test_mount_change_relabels_the_angle_fields(self, win):
        win.geometry._fields["mount"].setCurrentText("Classical")
        assert "α" in win.geometry._angle_labels["alpha"].text()

    def test_classical_hides_the_second_angle(self, win):
        """Checked with isHidden, not isVisible: under the offscreen platform
        an unshown window reports every child invisible, which would make this
        pass without testing anything."""
        win.geometry._fields["mount"].setCurrentText("Classical")
        assert win.geometry._fields["gamma"].isHidden()

    def test_off_plane_shows_the_second_angle(self, win):
        win.geometry._fields["mount"].setCurrentText("Off-plane")
        assert not win.geometry._fields["gamma"].isHidden()

    def test_profile_change_hides_irrelevant_fields(self, win):
        """duty_cycle is meaningless for a sinusoid. The mapping itself is
        tested headlessly against PROFILE_FIELDS."""
        win.geometry._fields["profile_kind"].setCurrentText("Sinusoidal")
        assert win.geometry._profile_rows["duty_cycle"].isHidden()
        assert not win.geometry._profile_rows["depth_fraction"].isHidden()

    def test_from_file_reveals_the_load_button(self, win):
        win.geometry._fields["profile_kind"].setCurrentText("From file")
        assert not win.geometry._load_button.isHidden()
        assert win.geometry._profile_rows["blaze_angle"].isHidden()


class TestBackgroundSolve:
    def test_the_solve_does_not_run_on_the_ui_thread(self, qtbot, win):
        """`moveToThread` only redirects signal delivery, so calling the
        worker's method directly would quietly solve on the UI thread and
        freeze the window -- the exact failure this design exists to avoid."""
        import threading

        seen = {}

        from gratinglab.gui.qt import worker as worker_module

        original_run = worker_module.SolveWorker.run

        def record(self, token, method, geometry, options):
            seen["thread"] = threading.current_thread().ident
            return original_run(self, token, method, geometry, options)

        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(worker_module.SolveWorker, "run", record)
            resolve(qtbot, win)

        assert "thread" in seen, "the worker slot never ran"
        assert seen["thread"] != threading.current_thread().ident

    def test_solve_is_disabled_while_one_is_running(self, win):
        win.solve()
        assert not win._solve_button.isEnabled()
        assert win._cancel_button.isEnabled()

    def test_a_second_solve_while_running_is_ignored(self, win):
        """One window, one result. Superseding would stack CPU-burning
        threads; queueing would surprise the user with a stale-then-fresh
        double redraw."""
        win.solve()
        token = win._token
        win.solve()
        assert win._token == token

    def test_cancel_keeps_the_previous_result_on_screen(self, win):
        previous = win._scan
        win.solve()
        win.cancel()
        assert win._scan is previous
        assert win._solve_button.isEnabled()

    def test_cancel_says_the_calculation_is_still_finishing(self, win):
        """Cancel abandons a result; it cannot stop the CPU. Claiming
        otherwise would be the panel's first lie."""
        win.solve()
        win.cancel()
        assert "still finishing" in win._provenance.toPlainText()

    def test_a_cancelled_result_is_dropped_when_it_arrives(self, win):
        from gratinglab.gui.qt.worker import SolveResult

        win.solve()
        stale = win._token
        win.cancel()
        # would crash if accepted
        win._on_solved(stale, "scalar", SolveResult(None, None))
        assert win._scan is not None

    def test_a_solver_exception_becomes_a_message_not_a_crash(self, win):
        """A Python exception escaping a slot under PySide6 can abort the
        process rather than surfacing."""
        win._token += 1
        win._set_running(True)
        win._on_failed(win._token, "scalar", "ValueError: contrived")
        assert "solve failed" in win._provenance.toPlainText()
        assert win._solve_button.isEnabled()

    def test_no_progress_bar_for_a_fast_solve(self, win):
        """The scalar solve takes about 70 ms; flashing a bar on and off
        within a frame reads as a glitch."""
        assert win._progress.isHidden()

    def test_closing_joins_the_worker_thread(self, qtbot):
        """Otherwise Qt reports 'QThread: Destroyed while thread is still
        running' and aborts -- most visibly under offscreen CI, where windows
        are built and torn down back to back.

        Constructed without qtbot.addWidget so that closing it here is the
        only close, and the assertion is about closeEvent rather than about
        teardown order.
        """
        from gratinglab.gui.qt.main_window import MainWindow

        window = MainWindow()
        assert window._thread.isRunning()
        window.close()
        assert not window._thread.isRunning()


class TestHelpMenu:
    """Read through `win._help_menu`, never `menuBar().actions()[i].menu()`.

    Under PySide6 6.11 a QMenu reached back down from the menu bar is not
    reliably usable -- it raises "Internal C++ object already deleted" while
    the menu the window built is alive and correct. Holding the QMenuBar does
    not help; holding the QMenu does. `_build_menu` therefore keeps both, and
    these tests use what the window kept.
    """

    def _labels(self, win, *, separators=False):
        return [
            action.text()
            for action in win._help_menu.actions()
            if separators or not action.isSeparator()
        ]

    def test_a_help_menu_is_attached_to_the_menu_bar(self, win):
        titles = [action.text() for action in win._menu_bar.actions()]
        assert titles == ["&Help"], "Help is meant to be the only menu"

    def test_lists_every_registered_solver_and_general_page(self, win):
        from gratinglab.gui.docs import general_pages, theory_pages
        from gratinglab.gui.richtext import menu_label

        labels = self._labels(win)
        for page in list(general_pages()) + list(theory_pages()):
            assert any(menu_label(page.title) in label for label in labels), page.title

    def test_an_ampersand_in_a_title_is_escaped_not_eaten(self, win):
        """Qt reads `&` as a mnemonic, so 'Grating Geometry & Conventions'
        would appear as 'Grating Geometry Conventions' with C underlined."""
        assert any("&&" in label for label in self._labels(win))

    def test_the_escaping_test_is_not_vacuous(self):
        """It only means something while a real page title contains `&`."""
        from gratinglab.gui.docs import general_pages

        assert any("&" in page.title for page in general_pages())

    def test_about_reports_the_real_version(self, win, monkeypatch):
        from gratinglab import __version__

        seen = {}
        monkeypatch.setattr(
            "PySide6.QtWidgets.QMessageBox.about",
            lambda parent, title, text: seen.update(title=title, text=text),
        )
        win.show_about()
        assert seen["title"] == "About GratingLab"
        assert __version__ in seen["text"]


class TestTheoryViewer:
    def _open(self, qtbot, win, name="scalar"):
        from gratinglab.gui.docs import theory_pages

        page = next(p for p in theory_pages() if p.name == name)
        viewer = win.show_theory(page)
        qtbot.addWidget(viewer)
        return viewer

    def test_opens_a_page_with_typeset_math(self, qtbot, win):
        text = self._open(qtbot, win)._browser.toPlainText()
        assert "5. Energy is not conserved" in text
        assert "$" not in text
        assert "\\Phi_m" not in text

    def test_the_equations_are_registered_as_images(self, qtbot, win):
        """The `Text.dump()` image count this replaces, expressed against the
        document's own resources."""
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QTextDocument

        document = self._open(qtbot, win)._browser.document()
        loaded = document.resource(
            QTextDocument.ResourceType.ImageResource, QUrl("gratinglab-math:0")
        )
        assert not loaded.isNull()

    def test_a_table_survives_into_the_document(self, qtbot, win):
        """What this milestone bought: scalar.md section 5 is a table, and it
        used to render as literal pipe characters."""
        assert "<table" in self._open(qtbot, win)._browser.toHtml()

    def test_a_non_rigorous_solver_is_bannered(self, qtbot, win):
        from gratinglab.gui.docs import theory_pages

        page = next(p for p in theory_pages() if p.name == "scalar")
        assert not page.rigorous  # the case this test needs
        assert "Approximate method" in self._open(qtbot, win)._browser.toPlainText()

    def test_a_page_that_does_not_exist_yet_explains_itself(self, qtbot, win):
        from gratinglab.gui.docs import TheoryPage

        stub = TheoryPage(
            name="future", title="Future Method", available=False, path=None,
            text="No theory page has been written yet for 'future'.",
            rigorous=True,
        )
        viewer = win.show_theory(stub)  # must not raise
        qtbot.addWidget(viewer)
        text = viewer._browser.toPlainText()
        assert "No theory page has been written yet" in text
        assert "Approximate method" not in text


class TestToolkitBoundary:
    def test_only_gui_qt_may_import_a_toolkit(self):
        """A package boundary rather than a convention, so the drift that
        always happens -- someone needing 'just a QColor' in a pure module --
        fails here instead of quietly spreading.

        `app.py` is the one file at this level that's expected to *reach into*
        the toolkit -- but only deferred, inside `_require_qt()`, which is
        exactly what keeps `import gratinglab.gui.app` safe without the `gui`
        extra. So this checks module-*level* imports only (`.body`, not a full
        recursive walk) -- a recursive walk would flag that sanctioned
        deferred import as if it were the unconditional kind this test exists
        to catch.
        """
        import ast
        from pathlib import Path

        import gratinglab.gui as gui

        root = Path(gui.__file__).parent
        offenders = {}
        for path in sorted(root.glob("*.py")):
            imported: set[str] = set()
            for node in ast.parse(path.read_text()).body:
                if isinstance(node, ast.Import):
                    imported.update(a.name.split(".")[0] for a in node.names)
                elif isinstance(node, ast.ImportFrom) and node.module:
                    imported.add(node.module.split(".")[0])
            found = {"PySide6", "PyQt6", "PyQt5", "tkinter"} & imported
            if found:
                offenders[path.name] = found

        assert offenders == {}, "only gui/qt/ may import a toolkit unconditionally"


class TestEntryPoint:
    """`app.py`'s own two jobs: explain a missing toolkit, and stay safe to
    import without one.

    Ported from the Tk suite this replaces -- PySide6 is pip-installable
    (unlike Tk), so the *message* changed, but the failure mode is still real:
    an install of `[dev]` without `[gui]` hits exactly this.
    """

    def test_gives_an_explanation_not_a_traceback(self, monkeypatch):
        import builtins

        from gratinglab.gui import app as app_module

        real_import = builtins.__import__

        def refuse(name, *args, **kwargs):
            if name == "PySide6.QtWidgets" or name.startswith("PySide6"):
                raise ImportError("No module named 'PySide6'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", refuse)

        with pytest.raises(SystemExit) as excinfo:
            app_module._require_qt()

        message = str(excinfo.value)
        assert 'pip install -e ".[gui]"' in message
        assert "Everything except the GUI works without it" in message

    def test_importing_the_module_does_not_require_pyside6(self):
        """Import must be safe where PySide6 is missing; only calling
        `_require_qt`/`main` needs it."""
        import ast
        import inspect

        from gratinglab.gui import app as app_module

        top_level = {
            node.names[0].name
            for node in ast.parse(inspect.getsource(app_module)).body
            if isinstance(node, ast.Import)
        }
        assert "PySide6" not in top_level
