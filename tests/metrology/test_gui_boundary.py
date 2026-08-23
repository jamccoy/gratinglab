"""The toolkit boundary, for the metrology window.

`gui/qt` is the only part of this GUI allowed to import a Qt binding. That is a
package boundary rather than a convention precisely so this test can enforce it.
The drift it prevents is the ordinary kind - someone needing "just a QColor" in
a module that was pure, and the headless tests quietly becoming un-runnable.

Was an AST scan of the source; now a subprocess import per module, via
`endstation.testing.boundary`. Stronger in two ways: it catches a *transitive*
toolkit import that a source scan cannot see, and it discovers the modules
rather than trusting a list, so a new pure module is covered the day it lands.

Needs no Qt installed: importing a pure module is the whole test, and a pure
module is one that does not need it.
"""
from endstation.testing.boundary import assert_qt_free, modules_in

PURE = modules_in("gratinglab.metrology.gui",
                  exclude=("gratinglab.metrology.gui.qt",))


def test_no_toolkit_outside_gui_qt():
    assert_qt_free(
        PURE,
        reason="only gratinglab/metrology/gui/qt may import a Qt binding.",
    )


def test_pure_modules_exist():
    """Guard against the walk silently finding nothing and passing vacuously"""
    assert "gratinglab.metrology.gui.state" in PURE, PURE


def test_the_gui_package_loads_neither_matplotlib_nor_a_toolkit():
    """The console script's fast path, guarded at both ends.

    This assertion has been through three states, and the history is the point.
    It first claimed `gratinglab.metrology.gui` did not reach matplotlib and
    passed for the wrong reason -- it read the import *statements* in
    `gui/__init__.py`, which indeed do not mention it. The subprocess check
    showed the truth: matplotlib arrived transitively through the parent
    package, whose `__init__` re-exported workflows that reach plotting code.
    So the assertion was inverted to match reality, with a note to update it if
    the parent ever got lazier.

    It has. `gratinglab.metrology` now hands its re-exports out through a
    module-level `__getattr__` (PEP 562), so the plotting stack is imported by
    the first caller who wants an analysis rather than by anyone who touches
    the package. `gratinglab-metrology-gui` exists partly to report a missing
    PySide6 quickly and clearly, and it no longer pays for matplotlib -- about
    a second -- before it can.

    What was never a defect, and still is not, is the matplotlib *requirement*:
    it is hard, that is what the `metrology` extra declares, and the parent
    still checks for it eagerly. Only the import moved, not the check --
    `test_the_parent_still_reports_a_missing_matplotlib_eagerly` holds that end
    down, and it is the one to read if this test starts failing.
    """
    from endstation.testing.boundary import assert_qt_free, leaked_imports

    assert leaked_imports("gratinglab.metrology.gui", ["matplotlib"]) == (), (
        "importing the GUI package should no longer pull in the plotting stack; "
        "something in the parent package's import path went eager again"
    )

    assert_qt_free(
        ["gratinglab.metrology.gui"],
        reason=("The console script touches this first and must be able to "
                "report a missing toolkit rather than fail on importing one."),
    )


def test_the_parent_still_reports_a_missing_matplotlib_eagerly():
    """Laziness moved the import, and must not have moved the check with it.

    The failure this guards against is quiet: if `require_matplotlib()` were
    ever deferred alongside the re-exports, a user missing the `metrology`
    extra would import the package cleanly and only find out four modules deep,
    in a traceback naming a package they never asked for. That is exactly the
    experience the eager check exists to prevent.

    Hides matplotlib from a fresh interpreter rather than uninstalling it, so
    this runs in an environment that has the extra -- which is every
    environment that can run the rest of this suite.
    """
    import subprocess
    import sys

    # The finder *raises* rather than returning None: None from a meta_path
    # finder means "not mine, ask the next one", which would leave the real
    # matplotlib perfectly importable and this test passing vacuously.
    script = """
import sys

class Block:
    def find_spec(self, name, path=None, target=None):
        if name == "matplotlib" or name.startswith("matplotlib."):
            raise ModuleNotFoundError(name)
        return None

sys.meta_path.insert(0, Block())
sys.modules.pop("matplotlib", None)

try:
    import gratinglab.metrology
except ModuleNotFoundError as err:
    print(err)
else:
    raise SystemExit("importing the parent package did not check for matplotlib")
"""
    probe = subprocess.run([sys.executable, "-c", script],
                           capture_output=True, text=True)
    assert probe.returncode == 0, probe.stderr
    assert "matplotlib" in probe.stdout
    assert ".[metrology]" in probe.stdout, probe.stdout


def test_the_parents_lazy_re_exports_still_resolve():
    """PEP 562 is easy to get subtly wrong; this is the blunt check that it isn't.

    A typo in the name-to-submodule table costs nothing at import time and
    raises `AttributeError` at the call site, which reads like the function was
    removed. Asking for each of them here turns that into a test failure.
    """
    import gratinglab.metrology as metrology

    for name in metrology.__all__:
        assert getattr(metrology, name) is not None, name


def test_gui_package_names_the_extra_in_its_error():
    """A missing toolkit should tell the user what to install"""
    from gratinglab.metrology.gui import QT_MISSING_MESSAGE
    assert 'PySide6' in QT_MISSING_MESSAGE
    assert '.[gui]' in QT_MISSING_MESSAGE
