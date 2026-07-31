"""The widget layer.

Barely tested on purpose -- the logic worth testing lives in
``gui/state.py``. What *is* checked here is that the module degrades politely
without a toolkit, and that the app is not a second source of truth: whatever
it plots must equal what the core produces for the same inputs.
"""

import builtins
import sys

import numpy as np
import pytest

from gratinglab.gui import app as app_module
from gratinglab.gui.state import FormState, build
from gratinglab.solvers import scalar

tk = pytest.importorskip("tkinter", reason="Tk not available in this Python")


def _display_available() -> bool:
    """Whether a Tk window can actually be created.

    Importing tkinter is not enough. On a headless Linux CI runner the module
    imports cleanly and ``Tk()`` then raises ``TclError: no display name``, so
    the guard has to attempt the thing itself.
    """
    try:
        root = tk.Tk()
    except tk.TclError:
        return False
    root.destroy()
    return True


requires_display = pytest.mark.skipif(
    not _display_available(), reason="no display; Tk cannot open a window"
)


class TestMissingToolkit:
    def test_gives_an_explanation_not_a_traceback(self, monkeypatch):
        """A user without Tk should be told how to fix it."""
        real_import = builtins.__import__

        def refuse(name, *args, **kwargs):
            if name == "tkinter":
                raise ImportError("No module named '_tkinter'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", refuse)
        monkeypatch.delitem(sys.modules, "tkinter", raising=False)

        with pytest.raises(SystemExit) as excinfo:
            app_module._require_tk()

        message = str(excinfo.value)
        assert "brew install python-tk" in message
        assert "Everything except the GUI works without it" in message

    def test_importing_the_module_does_not_require_tk(self):
        """Import must be safe even where Tk is missing; only calling needs it."""
        import ast
        import inspect

        top_level = {
            node.names[0].name
            for node in ast.parse(inspect.getsource(app_module)).body
            if isinstance(node, ast.Import)
        }
        assert "tkinter" not in top_level


@pytest.fixture
def gui():
    """A constructed app with no visible window."""
    import matplotlib

    matplotlib.use("Agg")
    root = tk.Tk()
    root.withdraw()
    instance = app_module.GratingLabApp(root)
    yield instance
    root.destroy()


@requires_display
class TestNotASecondSourceOfTruth:
    def test_plotted_values_equal_a_direct_solve(self, gui):
        """The single most important property of this layer.

        If the GUI ever computed anything itself, this would drift.
        """
        parsed = build(FormState())
        expected = scalar.solve(
            parsed.problem, parsed.illumination, parsed.wavelengths, **parsed.options
        )
        assert np.array_equal(gui._scan.efficiency, expected.efficiency)
        assert np.array_equal(gui._scan.orders, expected.orders)

    def test_changing_a_field_changes_the_result_consistently(self, gui):
        gui._vars["blaze_angle"].set("12.0")
        gui.solve()

        form = FormState(blaze_angle="12.0")
        expected = scalar.solve(
            *[getattr(build(form), n) for n in ("problem", "illumination", "wavelengths")],
            **build(form).options,
        )
        assert np.array_equal(gui._scan.efficiency, expected.efficiency)


@requires_display
class TestBehaviour:
    def test_solves_on_construction(self, gui):
        """Opening the app should show something, not an empty window."""
        assert gui._scan is not None
        assert len(gui._scan) == 200

    def test_invalid_input_reports_errors_without_crashing(self, gui):
        gui._vars["period"].set("not a number")
        gui.solve()
        text = gui._provenance.get("1.0", "end")
        assert "period" in text
        assert "need attention" in text

    def test_a_failed_solve_leaves_the_previous_result_intact(self, gui):
        good = gui._scan
        gui._vars["period"].set("")
        gui.solve()
        assert gui._scan is good

    def test_energy_violation_is_surfaced(self, gui):
        """The panel's reason for existing."""
        gui._vars["antiblaze_angle"].set("90")
        gui.solve()
        text = gui._provenance.get("1.0", "end")
        assert "EXCEEDS UNITY" in text
        assert "exceeding unity" in text

    def test_specular_phase_reference_clears_the_violation(self, gui):
        gui._vars["antiblaze_angle"].set("90")
        gui._vars["phase_reference"].set("specular")
        gui.solve()
        text = gui._provenance.get("1.0", "end")
        assert "EXCEEDS UNITY" not in text

    def test_convergence_is_reported_honestly(self, gui):
        """Nothing has been shown converged, and the panel must say so."""
        assert "not yet checked" in gui._provenance.get("1.0", "end")

    def test_convergence_none_is_not_tagged_as_a_warning(self, gui):
        """The concrete regression this milestone fixes."""
        ranges = gui._provenance.tag_ranges("warn")
        warned_text = "".join(
            gui._provenance.get(ranges[i], ranges[i + 1])
            for i in range(0, len(ranges), 2)
        )
        assert "not yet checked" not in warned_text
        assert "coating" not in warned_text

    def test_normalization_is_shown_neutrally_not_as_a_warning(self, gui):
        """No coating is the default, correct mode -- it must read as a plain
        status line, not an alarm."""
        text = gui._provenance.get("1.0", "end")
        assert "normalization: relative" in text

        ranges = gui._provenance.tag_ranges("warn")
        warned_text = "".join(
            gui._provenance.get(ranges[i], ranges[i + 1])
            for i in range(0, len(ranges), 2)
        )
        assert "normalization" not in warned_text
        assert "coating" not in warned_text

    def test_mount_change_relabels_the_angle_fields(self, gui):
        gui._vars["mount"].set("Classical")
        gui._on_mount_change()
        label, _ = gui._angle_labels["alpha"]
        assert "α" in label.cget("text")

    def test_profile_change_hides_irrelevant_fields(self, gui):
        gui._vars["profile_kind"].set("Sinusoidal")
        gui._on_profile_change()
        # duty_cycle is meaningless for a sinusoid
        assert not gui._profile_rows["duty_cycle"][1].winfo_ismapped()

    def test_reads_every_form_field(self, gui):
        """A field added to the UI but missing from FormState would break here."""
        form = gui._read_form()
        assert isinstance(form, FormState)
        assert form.period == "315.15"


@requires_display
class TestHelpMenu:
    def test_a_menu_bar_is_attached(self, gui):
        assert gui.root.cget("menu") != ""

    def test_theory_viewer_shows_a_known_heading(self, gui):
        """The '##' markdown prefix is stripped -- that is the fix -- so the
        heading's *text* must survive while its raw markdown syntax does not."""
        from gratinglab.gui.docs import theory_pages

        scalar_page = next(p for p in theory_pages() if p.name == "scalar")
        before = set(gui.root.winfo_children())
        gui.show_theory(scalar_page)
        opened = [w for w in gui.root.winfo_children() if w not in before]
        assert len(opened) == 1

        text_widget = next(
            w for w in opened[0].winfo_children() if isinstance(w, tk.Text)
        )
        content = text_widget.get("1.0", "end")
        assert "5. Energy is not conserved" in content
        assert "## 5. Energy is not conserved" not in content

    def test_theory_viewer_actually_typesets_math_not_raw_latex(self, gui):
        """The whole point of this milestone: no literal '$$' or backslash
        LaTeX source should reach the widget for a page that renders cleanly."""
        from gratinglab.gui.docs import theory_pages

        scalar_page = next(p for p in theory_pages() if p.name == "scalar")
        gui.show_theory(scalar_page)
        window = [w for w in gui.root.winfo_children() if isinstance(w, tk.Toplevel)][-1]
        text_widget = next(w for w in window.winfo_children() if isinstance(w, tk.Text))

        dump = text_widget.dump("1.0", "end")
        image_count = sum(1 for kind, _, _ in dump if kind == "image")
        assert image_count > 10, "expected many rendered equations, found none"

        content = text_widget.get("1.0", "end")
        assert "$$" not in content
        assert "\\Phi_m" not in content

    def test_general_page_math_also_renders(self, gui):
        """The generalized grating equation, opened from the new Help entry."""
        from gratinglab.gui.docs import general_pages

        page = next(p for p in general_pages() if p.name == "conventions")
        gui.show_theory(page)
        window = [w for w in gui.root.winfo_children() if isinstance(w, tk.Toplevel)][-1]
        text_widget = next(w for w in window.winfo_children() if isinstance(w, tk.Text))

        dump = text_widget.dump("1.0", "end")
        assert sum(1 for kind, _, _ in dump if kind == "image") > 0
        assert "$$" not in text_widget.get("1.0", "end")

    def test_theory_viewer_banners_a_non_rigorous_solver(self, gui):
        from gratinglab.gui.docs import theory_pages

        scalar_page = next(p for p in theory_pages() if p.name == "scalar")
        assert not scalar_page.rigorous  # the case this test actually needs

        gui.show_theory(scalar_page)
        window = [w for w in gui.root.winfo_children() if isinstance(w, tk.Toplevel)][-1]
        text_widget = next(w for w in window.winfo_children() if isinstance(w, tk.Text))
        assert "Approximate method" in text_widget.get("1.0", "end")

    def test_missing_page_viewer_shows_the_explanation_not_a_crash(self, gui):
        from gratinglab.gui.docs import TheoryPage

        stub = TheoryPage(
            name="future", title="Future Method", available=False, path=None,
            text="No theory page has been written yet for 'future'.",
            rigorous=True,
        )
        gui.show_theory(stub)  # must not raise
        window = [w for w in gui.root.winfo_children() if isinstance(w, tk.Toplevel)][-1]
        text_widget = next(w for w in window.winfo_children() if isinstance(w, tk.Text))
        assert "No theory page has been written yet" in text_widget.get("1.0", "end")
        # A solver marked rigorous gets no approximate-method banner.
        assert "Approximate method" not in text_widget.get("1.0", "end")

    def test_about_reports_the_real_version(self, gui, monkeypatch):
        seen = {}
        monkeypatch.setattr(
            gui._messagebox, "showinfo", lambda title, message: seen.update(
                title=title, message=message
            )
        )
        gui.show_about()

        from gratinglab import __version__

        assert seen["title"] == "About GratingLab"
        assert __version__ in seen["message"]
        assert "BSD-3-Clause" in seen["message"]

    def test_help_menu_lists_every_registered_solver_and_general_page(self, gui):
        from gratinglab.gui.docs import general_pages, theory_pages

        help_menu_index = gui.root.nametowidget(gui.root.cget("menu"))
        # Walk the Help cascade's item labels via the menu's own introspection.
        help_cascade = None
        for i in range(help_menu_index.index("end") + 1):
            if help_menu_index.type(i) == "cascade" and help_menu_index.entrycget(
                i, "label"
            ) == "Help":
                help_cascade = help_menu_index.nametowidget(
                    help_menu_index.entrycget(i, "menu")
                )
        assert help_cascade is not None

        labels = [
            help_cascade.entrycget(i, "label")
            for i in range(help_cascade.index("end") + 1)
            if help_cascade.type(i) == "command"
        ]
        for page in theory_pages():
            assert any(page.title in label for label in labels), page.title
        for page in general_pages():
            assert any(page.title in label for label in labels), page.title
