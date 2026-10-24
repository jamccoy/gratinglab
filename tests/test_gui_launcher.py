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


@pytest.fixture(autouse=True)
def _no_leaked_windows():
    """
    Clear the module-level registries around every test.

    `_OPEN_WINDOWS` and `_CHOOSER` exist to outlive the chooser, which is
    exactly what makes them a hazard here: a widget left in either one
    survives qtbot's teardown as a Python object whose C++ half is gone, and
    the next test to touch it fails somewhere unrelated.
    """
    qt_launcher._OPEN_WINDOWS.clear()
    qt_launcher._CHOOSER.clear()
    yield
    for window in list(qt_launcher._OPEN_WINDOWS) + list(qt_launcher._CHOOSER):
        window.close()
    qt_launcher._OPEN_WINDOWS.clear()
    qt_launcher._CHOOSER.clear()


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
                    window='PySide6.QtWidgets:QMainWindow',
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


class TestGettingBackAndSwitching:
    """
    The chooser was a one-way door: pick a mode and there was no route to the
    other one or back. Each window now carries a Go menu.
    """

    def _stub_mode(self, key='stub'):
        return Mode(key=key, title=f'Stub {key}', blurb='x' * 50,
                    window='PySide6.QtWidgets:QMainWindow',
                    script='gratinglab-app')

    def test_every_opened_window_gets_the_menu(self, qtbot):
        window = qt_launcher.open_window(self._stub_mode())
        qtbot.addWidget(window)

        menu = window._gratinglab_switch_menu
        assert menu.title() == qt_launcher.SWITCH_MENU_TITLE
        labels = [a.text() for a in window._gratinglab_switch_actions]
        assert labels[-1] == 'Main Menu'

    def test_the_menu_offers_the_other_modes_not_the_current_one(self, qtbot):
        """
        A window must not offer a route to itself.

        Checked with the real metrology mode, since that is the pairing that
        matters: from the AFM window the menu should name the efficiency one.
        """
        window = qt_launcher.open_window(mode('metrology'))
        qtbot.addWidget(window)

        labels = [a.text() for a in window._gratinglab_switch_actions]
        assert 'Compute efficiency' in labels
        assert 'Measure a scan' not in labels, "it offered a route to itself"

    def test_switching_reuses_the_window_rather_than_opening_a_second(self, qtbot):
        """
        The failure this prevents looks like losing work.

        The metrology window holds a loaded scan and its settings. A switch
        that produced an empty second copy would read as the scan being gone.
        """
        chosen = self._stub_mode()
        first = qt_launcher.open_window(chosen)
        qtbot.addWidget(first)

        again = qt_launcher.open_window(chosen)
        assert again is first
        assert len(qt_launcher._OPEN_WINDOWS) == 1

    def test_a_closed_window_comes_back_with_its_state(self, qtbot):
        """Closing is not destroying while a reference is held."""
        chosen = self._stub_mode()
        window = qt_launcher.open_window(chosen)
        qtbot.addWidget(window)
        window.close()
        assert not window.isVisible()

        again = qt_launcher.open_window(chosen)
        assert again is window
        assert again.isVisible()

    def test_main_menu_brings_back_the_same_chooser(self, qtbot):
        first = qt_launcher.show_chooser()
        qtbot.addWidget(first)
        first.close()

        again = qt_launcher.show_chooser()
        assert again is first, "a second chooser would be a second main menu"
        assert again.isVisible()
        assert len(qt_launcher._CHOOSER) == 1

    def test_the_round_trip(self, qtbot):
        """
        Chooser, into a mode, back to the chooser, into the other mode.

        The journey the Go menu exists for, with nothing stubbed but the
        windows themselves.
        """
        chooser = qt_launcher.show_chooser()
        qtbot.addWidget(chooser)

        first = qt_launcher.open_window(self._stub_mode('one'))
        qtbot.addWidget(first)
        chooser.close()

        back = qt_launcher.show_chooser()
        assert back is chooser and back.isVisible()

        second = qt_launcher.open_window(self._stub_mode('two'))
        qtbot.addWidget(second)

        assert first is not second
        assert len(qt_launcher._OPEN_WINDOWS) == 2

    def test_the_main_menu_action_actually_shows_it(self, qtbot):
        """Triggering the menu action, not calling the function behind it."""
        window = qt_launcher.open_window(self._stub_mode())
        qtbot.addWidget(window)

        assert not qt_launcher._CHOOSER
        back = window._gratinglab_switch_actions[-1]
        back.trigger()

        assert len(qt_launcher._CHOOSER) == 1
        qtbot.addWidget(qt_launcher._CHOOSER[0])
        assert qt_launcher._CHOOSER[0].isVisible()

    def test_the_other_mode_action_opens_that_mode(self, qtbot):
        from gratinglab.metrology.gui.qt.main_window import (
            MainWindow as MetrologyWindow,
        )

        window = qt_launcher.open_window(mode('efficiency'))
        qtbot.addWidget(window)

        switch = window._gratinglab_switch_actions[0]
        assert switch.text() == 'Measure a scan'
        switch.trigger()

        opened = qt_launcher._OPEN_WINDOWS[-1]
        qtbot.addWidget(opened)
        assert isinstance(opened, MetrologyWindow)
