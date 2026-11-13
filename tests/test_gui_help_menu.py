"""The Help menu, now in every window rather than only the efficiency one.

The pages the menu shows are the project's normative documents: the angle
conventions every reported number is in, and one page per solver. They were
reachable only from the efficiency window, which left the AFM window -- the
half in regular use, and the one whose own analysis `theory/metrology.md`
describes -- with no route to any of them.

What each page says is tested in `test_gui_docs.py` and `test_gui_qt.py`;
what is asserted here is only that every window offers them, and that nothing
has to know which window it is in to do so.
"""

import pytest

from gratinglab.gui.docs import general_pages, theory_pages
from gratinglab.gui.richtext import menu_label

pytest.importorskip(
    "PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"'
)

from gratinglab.gui.qt.help_menu import HELP_MENU_TITLE  # noqa: E402

#: Every window that carries a menu bar, as "module:attribute".
WINDOWS = {
    "efficiency": "gratinglab.gui.qt.main_window:MainWindow",
    "metrology": "gratinglab.metrology.gui.qt.main_window:MainWindow",
}


def labels_of(menu, *, separators=False):
    return [
        action.text()
        for action in menu.actions()
        if separators or not action.isSeparator()
    ]


@pytest.mark.gui
class TestEveryWindowOffersIt:
    @pytest.fixture(params=sorted(WINDOWS), ids=sorted(WINDOWS))
    def window(self, request, qtbot):
        from importlib import import_module

        module_name, _, attribute = WINDOWS[request.param].partition(":")
        window = getattr(import_module(module_name), attribute)()
        qtbot.addWidget(window)
        yield window
        window.close()

    def test_the_menu_bar_carries_help_last(self, window):
        """Last, so the bar reads the same way in both: View, then Help."""
        titles = [action.text() for action in window._menu_bar.actions()]
        assert titles[-1] == HELP_MENU_TITLE
        assert "&View" in titles

    def test_about_comes_first(self, window):
        assert labels_of(window._help_menu)[0] == "About GratingLab"

    def test_every_page_is_listed(self, window):
        """The same pages in both windows: the menu is built from
        `gratinglab.gui.docs`, so a page added there appears in both at once
        and neither window holds its own list."""
        shown = labels_of(window._help_menu)
        for page in (*general_pages(), *theory_pages()):
            assert any(menu_label(page.title) in label for label in shown), page.title

    def test_about_reports_the_real_version(self, window, monkeypatch):
        from gratinglab import __version__

        seen = {}
        monkeypatch.setattr(
            "PySide6.QtWidgets.QMessageBox.about",
            lambda parent, title, text: seen.update(title=title, text=text),
        )
        window._help_menu.actions()[0].trigger()
        assert seen["title"] == "About GratingLab"
        assert __version__ in seen["text"]


@pytest.mark.gui
class TestOpeningAPageFromTheAfmWindow:
    """The route that did not exist, driven end to end.

    One window and one page rather than the full matrix: the viewer itself is
    tested in `test_gui_qt.py`, and what is new here is that the AFM window
    can reach it at all.
    """

    def test_the_conventions_page_opens(self, qtbot):
        from gratinglab.gui.qt.theory_viewer import TheoryViewer
        from gratinglab.metrology.gui.qt.main_window import MainWindow

        window = MainWindow()
        qtbot.addWidget(window)
        try:
            wanted = menu_label(general_pages()[0].title)
            action = next(
                a for a in window._help_menu.actions() if a.text() == wanted
            )
            action.trigger()

            viewers = window.findChildren(TheoryViewer)
            assert len(viewers) == 1, "the menu entry opened no page"
            qtbot.addWidget(viewers[0])
            assert viewers[0].isVisible()
        finally:
            window.close()
