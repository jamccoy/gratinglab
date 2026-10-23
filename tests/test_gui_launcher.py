"""The chooser: one front door onto the two windows.

The failure this guards against is the one that prompted it. Both halves have
always been one package, but the desktop layer shipped two scripts and two
icons, one named after the repository the metrology was merged from, so a user
who knew the project well went looking for the AFM window inside the efficiency
one. A chooser only helps if it keeps naming both halves, and if clicking
really opens the window rather than appearing to.

The mode table is pure and tested without a toolkit; only the widget needs Qt.
"""

import tomllib
from importlib import import_module
from pathlib import Path

import pytest

from gratinglab.gui.launcher import MODES, Mode, mode

REPO = Path(__file__).resolve().parent.parent


class TestTheModeTable:
    """Pure: no Qt, no window, no display."""

    def test_both_halves_are_offered(self):
        keys = [m.key for m in MODES]
        assert keys == ['metrology', 'efficiency'], (
            "the metrology window comes first, being the half in regular use")

    def test_every_mode_is_fully_described(self):
        for chosen in MODES:
            assert chosen.title and not chosen.title.endswith('.')
            assert len(chosen.blurb) > 40, f"{chosen.key} has no real blurb"
            assert ':' in chosen.window, chosen.window
            assert chosen.script.startswith('gratinglab-')

    def test_keys_and_titles_are_unique(self):
        assert len({m.key for m in MODES}) == len(MODES)
        assert len({m.title for m in MODES}) == len(MODES)

    def test_lookup_by_key(self):
        assert mode('metrology').script == 'gratinglab-metrology-gui'
        with pytest.raises(KeyError, match='no such mode'):
            mode('nope')

    def test_no_mode_names_the_predecessor_project(self):
        """
        The AFM half is not a separate application and must not read as one.

        Its bundle used to be called "AFM Blaze Meas" with the bundle id
        org.jamccoy.afmblazemeas, inherited from the repository it was merged
        from. That naming is what made the half in regular use the half nobody
        could find, so it does not come back through the chooser.
        """
        text = ' '.join(m.title + m.blurb + m.window + m.script
                        for m in MODES).lower()
        assert 'afm blaze meas' not in text
        assert 'afmblazemeas' not in text
        assert 'jamccoy' not in text

    def test_each_script_is_a_real_entry_point(self):
        """A blurb naming a command that does not exist is worse than none."""
        pyproject = tomllib.loads(
            (REPO / 'pyproject.toml').read_text(encoding='utf-8'))
        declared = pyproject['project']['scripts']
        for chosen in MODES:
            assert chosen.script in declared, chosen.script
        assert 'gratinglab-app' in declared, "the front door is not installed"
        assert declared['gratinglab-app'] == 'gratinglab.gui.launcher:main'


def test_the_launcher_module_needs_no_toolkit():
    """
    `gratinglab.gui.launcher` sits outside `gui/qt`, so it must stay pure.

    This mirrors the rule the rest of the package follows, and it is why the
    mode table names each window as a string instead of importing it.
    """
    from endstation.testing.boundary import assert_qt_free

    assert_qt_free(
        ['gratinglab.gui.launcher'],
        reason='only gratinglab/gui/qt may import a Qt binding.',
    )


pytest.importorskip(
    "PySide6", reason='Qt not installed; pip install -e ".[dev,gui]"'
)

from PySide6.QtWidgets import QLabel  # noqa: E402

from gratinglab.gui.qt import launcher as qt_launcher  # noqa: E402


class TestTheWindowsItNames:
    """The paths in the table have to resolve, without building a window."""

    @pytest.mark.parametrize('chosen', MODES, ids=lambda m: m.key)
    def test_the_window_path_resolves_to_a_class(self, chosen):
        module_name, _, attribute = chosen.window.partition(':')
        found = getattr(import_module(module_name), attribute, None)
        assert isinstance(found, type), (
            f"{chosen.key} names {chosen.window}, which is not a class")


class TestTheChooser:
    def test_there_is_a_button_per_mode(self, qtbot):
        window = qt_launcher.LauncherWindow()
        qtbot.addWidget(window)

        assert set(window.buttons) == {m.key for m in MODES}
        for chosen in MODES:
            assert window.buttons[chosen.key].text() == chosen.title

    def test_both_halves_are_named_in_the_window(self, qtbot):
        """
        Whatever else changes, the word AFM has to be on this window.

        Its absence from the efficiency window is the whole reason the chooser
        exists.
        """
        window = qt_launcher.LauncherWindow()
        qtbot.addWidget(window)

        shown = ' '.join(label.text() for label in window.findChildren(QLabel))
        shown += ' '.join(b.text() for b in window.buttons.values())
        assert 'AFM' in shown

    def test_clicking_opens_that_mode_and_closes_the_chooser(self, qtbot,
                                                             monkeypatch):
        opened = []
        monkeypatch.setattr(qt_launcher, 'open_window', opened.append)

        window = qt_launcher.LauncherWindow()
        qtbot.addWidget(window)
        window.show()

        window.buttons['metrology'].click()

        assert [m.key for m in opened] == ['metrology']
        assert not window.isVisible(), "the chooser should get out of the way"

    def test_the_chooser_closes_only_after_the_window_is_up(self, qtbot,
                                                           monkeypatch):
        """
        Order matters: Qt quits when the last window closes.

        Closing the chooser first would end the process before the window it
        was asked for existed, which is a blank launch rather than an error.
        """
        order = []
        monkeypatch.setattr(qt_launcher, 'open_window',
                            lambda chosen: order.append('opened'))

        window = qt_launcher.LauncherWindow()
        qtbot.addWidget(window)
        window.show()
        original_close = window.close

        def recording_close():
            order.append('closed')
            return original_close()

        monkeypatch.setattr(window, 'close', recording_close)
        window.buttons['efficiency'].click()

        assert order == ['opened', 'closed']

    def test_an_opened_window_is_kept_alive(self, qtbot):
        """
        A top-level Qt window with no Python reference is collected and
        disappears a moment after it appears, which reads as a dead button.

        Checked with a bare QWidget rather than either real window, because
        what is under test is the bookkeeping, not the windows.
        """
        before = len(qt_launcher._OPEN_WINDOWS)
        stub = Mode(key='stub', title='Stub', blurb='x' * 50,
                    window='PySide6.QtWidgets:QWidget',
                    script='gratinglab-app')

        opened = qt_launcher.open_window(stub)
        qtbot.addWidget(opened)
        try:
            assert len(qt_launcher._OPEN_WINDOWS) == before + 1
            assert qt_launcher._OPEN_WINDOWS[-1] is opened
        finally:
            qt_launcher._OPEN_WINDOWS.remove(opened)


class TestItReallyOpensTheRealWindows:
    """
    The integration the rest of this file stubs out.

    Everything above monkeypatches `open_window`, which checks the wiring and
    not the thing most likely to break: constructing either real window inside
    the chooser's own QApplication. Both windows used to create their own, and
    both set QT_API before matplotlib's shim loads, so in-process construction
    is where an ordering problem would show up.

    Slow -- it builds a full window per mode -- and worth it, because the
    failure it catches is a button that does nothing.
    """

    @pytest.mark.parametrize('chosen', MODES, ids=lambda m: m.key)
    def test_the_window_opens_in_this_process(self, qtbot, chosen):
        before = len(qt_launcher._OPEN_WINDOWS)

        window = qt_launcher.open_window(chosen)
        qtbot.addWidget(window)
        try:
            module_name, _, attribute = chosen.window.partition(':')
            expected = getattr(import_module(module_name), attribute)
            assert isinstance(window, expected)
            assert window.isVisible()
            assert qt_launcher._OPEN_WINDOWS[-1] is window
        finally:
            qt_launcher._OPEN_WINDOWS.remove(window)
            window.close()
        assert len(qt_launcher._OPEN_WINDOWS) == before

    def test_the_chooser_opens_the_metrology_window_for_real(self, qtbot):
        """
        The click path end to end, for the half that was unfindable.

        No monkeypatch: the button is clicked and a metrology window has to
        exist afterwards.
        """
        from gratinglab.metrology.gui.qt.main_window import (
            MainWindow as MetrologyWindow,
        )

        chooser = qt_launcher.LauncherWindow()
        qtbot.addWidget(chooser)
        chooser.show()

        before = len(qt_launcher._OPEN_WINDOWS)
        chooser.buttons['metrology'].click()
        try:
            assert len(qt_launcher._OPEN_WINDOWS) == before + 1
            opened = qt_launcher._OPEN_WINDOWS[-1]
            qtbot.addWidget(opened)
            assert isinstance(opened, MetrologyWindow)
            assert not chooser.isVisible()
        finally:
            opened = qt_launcher._OPEN_WINDOWS.pop()
            opened.close()
