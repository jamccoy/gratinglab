"""One way in: the chooser that asks which window you want.

The package has always held both halves -- `gratinglab.metrology` measures a
scan and the solvers compute an efficiency from the result -- but the desktop
layer shipped them as two console scripts and so two application bundles, one
of them still named after the repository the metrology was merged from. The
code was one package and the Dock said otherwise, which cost a user who knew
the project well twenty minutes looking for the AFM window inside the
efficiency one.

This module is that single front door. It is deliberately thin: it owns the
list of windows and nothing else, and each window is named as a
``"module:attribute"`` string that is imported only once chosen, so starting
the chooser does not pay for both windows.

Like :mod:`gratinglab.gui.app`, its PySide6 import is deferred into
:func:`_require_qt`, which keeps this module importable on an install that
never asked for the ``gui`` extra. The widget lives in
:mod:`gratinglab.gui.qt.launcher`, because only ``gui/qt`` may touch a
toolkit and a test enforces it.

Run with ``gratinglab-app``. The two direct scripts still work and are
unchanged, for a shortcut or a habit.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass

from . import exit_without_qt

__all__ = ["Mode", "MODES", "mode", "main"]


@dataclass(frozen=True)
class Mode:
    """One window the chooser can open.

    ``window`` is a ``"module:attribute"`` path rather than the class itself,
    so this module stays free of Qt and of the cost of importing either
    window. ``script`` is the console script that opens the same window
    directly, which the chooser shows so that a reader learns the shortcut.
    """

    key: str
    title: str
    blurb: str
    window: str
    script: str


#: The metrology window comes first because it is the half in regular use.
#: Adding a third window is one entry here and nothing else.
MODES = (
    Mode(
        key="metrology",
        title="Measure a scan",
        blurb="Import an AFM scan, fit a blaze angle to every groove with "
              "row-group statistics, and average the grooves into a boundary "
              "profile to export.",
        window="gratinglab.metrology.gui.qt.main_window:MainWindow",
        script="gratinglab-metrology-gui",
    ),
    Mode(
        key="efficiency",
        title="Compute efficiency",
        blurb="Set the grating and the illumination, then run the scalar and "
              "integral solvers on that one problem and compare what they "
              "return.",
        window="gratinglab.gui.qt.main_window:MainWindow",
        script="gratinglab-gui",
    ),
)


def mode(key: str) -> Mode:
    """One mode by key. Raises KeyError, since a missing key is a typo."""
    for candidate in MODES:
        if candidate.key == key:
            return candidate
    raise KeyError(f"no such mode {key!r}; have {[m.key for m in MODES]}")


def _require_qt():
    """Import PySide6 with an explanation instead of a traceback."""
    exit_without_qt()
    from PySide6.QtWidgets import QApplication

    return QApplication


def main(argv=None) -> int:
    """Console-script entry point for ``gratinglab-app``."""
    QApplication = _require_qt()
    from .qt.launcher import show_chooser

    app = QApplication.instance() or QApplication(list(argv or sys.argv))
    # Fusion for the same reason the metrology window asks for it: the native
    # macOS style gives a QFormLayout row heights this layout does not expect.
    # `install` sets it, and applies the saved theme to both windows.
    from .qt.theming import install

    install(app)

    # show_chooser, not LauncherWindow(), so that "Main Menu" in either
    # window brings this same chooser back instead of a second copy.
    show_chooser()
    return app.exec()


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
