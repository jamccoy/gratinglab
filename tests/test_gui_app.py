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
        assert "not demonstrated" in gui._provenance.get("1.0", "end")

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
