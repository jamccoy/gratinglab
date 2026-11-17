"""'Copy as Python': the scripts are run, and must give what the window gives.

A script that reads well and computes something else would be worse than no
script, so each one here is executed and its numbers compared with solving the
same form directly. The import surface the scripts rely on is checked too.
"""

import importlib
import subprocess
import sys

import numpy as np
import pytest

from gratinglab import solve
from gratinglab.compare import align, sweep
from gratinglab.gui.mirror_view import MirrorSettings, compute
from gratinglab.gui.script import efficiency_script, mirror_script
from gratinglab.gui.state import FormState, build
from gratinglab.io.ggp import write_ggp

#: A small problem every solver takes: no vertical facet, three wavelengths.
FORM = FormState(coating="typical:Au", roughness="0.3", wavelength_count="3")


def run(script: str) -> dict:
    namespace: dict = {}
    exec(compile(script, "<copied script>", "exec"), namespace)
    return namespace


def default_options(method: str, parsed) -> dict:
    """The settings that solver's tab would pass with its fields untouched."""
    module = importlib.import_module(f"gratinglab.gui.{method}_options")
    state = next(v for k, v in vars(module).items() if k.endswith("OptionsState"))
    return module.build_options(parsed.problem, parsed.illumination, parsed.wavelengths, state())


class TestTheGratingScript:
    @pytest.mark.parametrize("method", ["scalar", "rcwa", "cmethod", "integral"])
    def test_it_reproduces_the_solve(self, method, capsys):
        parsed = build(FORM)
        options = default_options(method, parsed)
        script = efficiency_script(parsed, [method], {method: options})
        got = run(script)["scan"]
        expect = solve(parsed.problem, parsed.illumination, parsed.wavelengths, method=method, **options)
        assert np.array_equal(got.orders, expect.orders)
        assert np.allclose(got.efficiency, expect.efficiency, rtol=1e-12, atol=0.0)
        assert "order +1: peak" in capsys.readouterr().out

    def test_the_problem_is_what_the_form_built(self):
        parsed = build(FormState(coating="typical:Si", layers="typical:Ni 5, typical:Au 20",
                                 roughness="0.3", wavelength_count="4"))
        got = run(efficiency_script(parsed))
        assert got["problem"] == parsed.problem
        assert got["illumination"] == parsed.illumination
        assert np.array_equal(got["wavelengths"], parsed.wavelengths)

    def test_defaults_are_left_out_and_the_problem_is_spread_over_lines(self):
        script = efficiency_script(build(FORM))
        assert "n_grooves" not in script and "substrate" not in script and "layers" not in script
        assert "problem = Problem(\n    period=315.15,\n" in script

    def test_several_methods_are_compared(self, capsys):
        parsed = build(FORM)
        options = {"scalar": default_options("scalar", parsed)}
        got = run(efficiency_script(parsed, ["scalar", "integral"], options))
        expect = align(sweep(parsed.problem, parsed.illumination, parsed.wavelengths,
                             ["scalar", "integral"], options))
        assert got["comparison"].summary("integral", "scalar") == expect.summary("integral", "scalar")
        assert "integral against scalar" in capsys.readouterr().out

    def test_a_measured_groove_is_loaded_again_not_written_out(self, tmp_path):
        path = write_ggp(tmp_path / "groove.ggp", t=[0.0, 0.4, 1.0], y=[0.0, 0.25, 0.0])
        form = FormState(profile_kind="From file", profile_path=str(path), wavelength_count="3")
        parsed = build(form)
        script = efficiency_script(parsed, ["scalar"], profile_path=form.profile_path)
        assert "read_profile_file(" in script and "groove.ggp" in script
        got = run(script)
        assert got["problem"] == parsed.problem
        expect = solve(parsed.problem, parsed.illumination, parsed.wavelengths)
        assert np.allclose(got["scan"].efficiency, expect.efficiency, rtol=1e-12, atol=0.0)

    def test_an_uneven_scan_is_written_out(self):
        parsed = build(FORM)
        uneven = type(parsed)(parsed.problem, parsed.illumination, np.array([1.0, 1.5, 4.0]))
        got = run(efficiency_script(uneven))
        assert list(got["wavelengths"]) == [1.0, 1.5, 4.0]


class TestTheMirrorScript:
    @pytest.mark.parametrize("settings", [
        MirrorSettings(graze_count="7"),
        MirrorSettings(scan="wavelength", graze="1.5", polarization="s"),
    ])
    def test_it_reproduces_the_tab(self, settings, capsys):
        form = FormState(coating="typical:Si", layers="typical:Ni 5, typical:Au 20",
                         roughness="0.3", wavelength_count="6")
        result = compute(form, settings)
        got = run(mirror_script(result))
        assert np.allclose(got["reflectivity"], result.reflectivity, rtol=1e-12, atol=0.0)
        assert np.allclose(got["depth"]["typical:Au"], result.depths["typical:Au"], rtol=1e-12)
        assert got["needed"] == pytest.approx(float(result.bulk[0]))
        assert "needed to read as bulk" in capsys.readouterr().out

    def test_a_bare_substrate_has_no_layer_to_test(self):
        script = mirror_script(compute(FormState(coating="typical:Au"), MirrorSettings(graze_count="5")))
        assert "bulk_thickness" not in script and "ConformalLayer" not in script
        run(script)


class TestTheTopLevelImport:
    def test_a_script_needs_one_import_line(self):
        import gratinglab

        for name in ("Problem", "ConformalLayer", "Blazed", "Lamellar", "Sinusoidal", "Trapezoidal",
                     "Illumination", "solve", "get_solver", "available_solvers", "sweep", "align",
                     "Mirror", "materials"):
            assert getattr(gratinglab, name) is not None, name
            assert name in gratinglab.__all__ and name in dir(gratinglab)

    def test_an_unknown_name_is_still_an_attribute_error(self):
        import gratinglab

        with pytest.raises(AttributeError, match="no attribute 'solv'"):
            gratinglab.solv

    def test_importing_the_package_does_not_load_the_solvers(self):
        """The lazy names keep `import gratinglab` as cheap as it was: the
        solvers pull in scipy's linear algebra, which a script that only
        reads a material should not pay for."""
        code = (
            "import sys, gratinglab; "
            "loaded = [m for m in ('gratinglab.solvers', 'gratinglab.compare', 'gratinglab.mirror') if m in sys.modules]; "
            "assert not loaded, loaded; "
            "from gratinglab import solve; assert 'gratinglab.solvers' in sys.modules; print('ok')"
        )
        done = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, encoding="utf-8"
        )
        assert done.returncode == 0 and done.stdout.strip() == "ok", done.stderr


@pytest.mark.gui
class TestTheMenuActions:
    @pytest.fixture
    def win(self, qtbot):
        pytest.importorskip("PySide6")
        from gratinglab.gui.qt.main_window import MainWindow

        window = MainWindow()
        qtbot.addWidget(window)
        with qtbot.waitSignal(window.solved, timeout=60_000, raising=False):
            pass
        return window

    def test_a_solver_tab_copies_that_solver_with_its_settings(self, win):
        from PySide6.QtWidgets import QApplication

        win._tab_widget.setCurrentWidget(win.tabs["scalar"])
        win.copy_as_python()
        script = QApplication.clipboard().text()
        assert "method='scalar'" in script and "quadrature_points=" in script
        assert "Copied" in win.statusBar().currentMessage()
        run(script)

    def test_the_geometry_tab_copies_the_problem_alone(self, win):
        win._tab_widget.setCurrentWidget(win.geometry_tab)
        script = win.python_script()
        assert "solve(" not in script and "problem = Problem(" in script

    def test_the_compare_tab_copies_the_ticked_solvers(self, win):
        win._tab_widget.setCurrentWidget(win.compare_tab)
        win.compare_tab.set_selected(["scalar", "integral"])
        script = win.python_script()
        assert "methods = ['scalar', 'integral']" in script or "methods = ['integral', 'scalar']" in script

    def test_a_form_that_does_not_parse_copies_nothing(self, win, monkeypatch):
        from PySide6.QtWidgets import QApplication, QMessageBox

        shown = []
        monkeypatch.setattr(QMessageBox, "warning", lambda *args: shown.append(args[2]))
        QApplication.clipboard().setText("untouched")
        win.geometry._fields["period"].setText("-")
        win.copy_as_python()
        assert QApplication.clipboard().text() == "untouched"
        assert "period" in shown[0]

    def test_the_coatings_window_copies_its_mirror(self, qtbot):
        import dataclasses

        from PySide6.QtWidgets import QApplication

        from gratinglab.gui.qt.coatings_window import CoatingsWindow

        coat = CoatingsWindow()
        qtbot.addWidget(coat)
        coat.copy_as_python()  # no substrate yet
        assert "Nothing copied" in coat.statusBar().currentMessage()
        coat.stack.write_form(dataclasses.replace(
            coat.stack.read_form(), coating="typical:Si", layers="typical:Au 20"))
        coat.mirror_tab._fields["graze_count"].setText("5")
        coat.copy_as_python()
        script = QApplication.clipboard().text()
        assert "mirror = Mirror(substrate='typical:Si'" in script
        run(script)
