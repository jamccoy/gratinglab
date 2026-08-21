"""
The toolkit boundary.

`gui/qt` is the only part of the GUI allowed to import a Qt binding. That is a
package boundary rather than a convention precisely so this test can enforce it.
The drift it prevents is the ordinary kind - someone needing "just a QColor" in a
module that was pure, and the headless tests quietly becoming un-runnable.

Needs no Qt installed: it reads source, it does not import widgets.
"""
import ast
import pathlib

GUI_ROOT = pathlib.Path(__file__).parent.parent.parent / 'src' / 'gratinglab' / 'metrology' / 'gui'
TOOLKITS = ('PySide6', 'PySide2', 'PyQt5', 'PyQt6')


def _imported_modules(path):
    tree = ast.parse(path.read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                yield alias.name
        elif isinstance(node, ast.ImportFrom) and node.module:
            yield node.module


def _pure_modules():
    """Every GUI module outside gui/qt"""
    qt_dir = GUI_ROOT / 'qt'
    for path in GUI_ROOT.rglob('*.py'):
        if qt_dir not in path.parents and path.parent != qt_dir:
            yield path


def test_no_toolkit_outside_gui_qt():
    offenders = []
    for path in _pure_modules():
        for module in _imported_modules(path):
            if module.split('.')[0] in TOOLKITS:
                offenders.append(f"{path.relative_to(GUI_ROOT.parent)} imports {module}")
    assert not offenders, (
        "only gratinglab/metrology/gui/qt may import a Qt binding:\n  "
        + "\n  ".join(offenders))


def test_pure_modules_exist():
    """Guard against the walk silently finding nothing and passing vacuously"""
    names = {p.name for p in _pure_modules()}
    assert 'state.py' in names, f"expected state.py among {names}"


def test_matplotlib_is_not_imported_at_gui_package_import():
    """
    Importing gratinglab.metrology.gui must stay cheap and display-free.

    It is what the `gratinglab-metrology-gui` console script touches first, and what raises the
    friendly message when PySide6 is missing; pulling matplotlib in at that point
    would make a missing-toolkit error slower and noisier than it needs to be.
    """
    for module in _imported_modules(GUI_ROOT / '__init__.py'):
        assert not module.startswith('matplotlib'), \
            f"gui/__init__.py imports {module} at module level"


def test_gui_package_names_the_extra_in_its_error():
    """A missing toolkit should tell the user what to install"""
    from gratinglab.metrology.gui import QT_MISSING_MESSAGE
    assert 'PySide6' in QT_MISSING_MESSAGE
    assert '.[gui]' in QT_MISSING_MESSAGE
