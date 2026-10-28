"""Themes: the palettes as data, the order colours, and switching live.

The contrast thresholds are WCAG 2.1's: 4.5:1 for anything read as text and
3:1 for plotted lines. They are asserted for every theme, so a palette edit
that makes a colour unreadable fails here rather than on somebody's screen.
"""

from pathlib import Path

import pytest

from gratinglab.gui import diagram, provenance
from gratinglab.gui.efficiency_plot import draw_efficiency
from gratinglab.gui.theme import (
    LIGHT,
    MIN_LINE_CONTRAST,
    THEMES,
    contrast,
    order_colors,
    provenance_colors,
)

ALL = list(THEMES.values())

#: The committed synthetic scan, for the AFM window.
SAMPLE = str(Path(__file__).parent / "metrology" / "fixtures" / "synthetic_blazed_scan.txt")


class TestPalettes:
    @pytest.mark.parametrize("theme", ALL, ids=lambda t: t.key)
    def test_text_is_readable(self, theme):
        for name in ("text", "muted", "link", "warn", "bad", "ok"):
            ratio = contrast(getattr(theme, name), theme.base)
            assert ratio >= 4.5, f"{theme.key}.{name} is {ratio:.2f}:1 on base"
        assert contrast(theme.highlight_text, theme.highlight) >= 4.5
        assert contrast(theme.ink, theme.axes) >= 4.5

    def test_light_is_the_look_before_themes(self):
        """Choosing Light must change nothing but the order colours."""
        assert provenance_colors(LIGHT) == provenance.TAG_COLORS
        assert dict(LIGHT.diagram) == diagram.TAG_COLORS

    @pytest.mark.parametrize("theme", ALL, ids=lambda t: t.key)
    def test_every_diagram_tag_is_covered(self, theme):
        assert set(theme.diagram) == set(diagram.TAG_COLORS)


class TestOrderColors:
    ORDERS = range(-8, 12)

    @pytest.mark.parametrize("theme", ALL, ids=lambda t: t.key)
    def test_every_order_is_distinct_and_legible(self, theme):
        colors = order_colors(self.ORDERS, theme)
        assert len(set(colors.values())) == len(colors) == 20
        worst = min(contrast(c, theme.axes) for c in colors.values() if c != theme.muted)
        assert worst >= MIN_LINE_CONTRAST

    def test_the_old_cycle_collision_is_gone(self):
        """The default ten-colour cycle drew m=-3 and m=+7 alike."""
        colors = order_colors(self.ORDERS, LIGHT)
        assert colors[-3] != colors[7]
        assert colors[-1] != colors[9]

    def test_sign_reads_from_hue(self):
        import colorsys

        def hue(color):
            value = color.lstrip("#")
            rgb = (int(value[i:i + 2], 16) / 255 for i in (0, 2, 4))
            return colorsys.rgb_to_hls(*rgb)[0] * 360

        colors = order_colors(self.ORDERS, LIGHT)
        assert all(170 <= hue(colors[m]) <= 290 for m in range(-8, 0))
        assert all(hue(colors[m]) <= 60 or hue(colors[m]) >= 300 for m in range(1, 12))
        assert colors[0] == LIGHT.muted

    def test_a_colour_depends_only_on_the_set_of_orders(self):
        assert order_colors(self.ORDERS, LIGHT) == order_colors(
            reversed(list(self.ORDERS)), LIGHT
        )


class TestEfficiencyPlot:
    def test_legend_sits_outside_the_axes(self):
        import numpy as np
        from matplotlib.figure import Figure

        from gratinglab.illumination import Illumination
        from gratinglab.problem import Problem
        from gratinglab.profiles import Blazed
        from gratinglab.solvers import scalar

        scan = scalar.solve(
            Problem(period=315.15, profile=Blazed(blaze_angle=29.5, antiblaze_angle=70.5)),
            Illumination.offplane(graze=1.5, azimuth=25.0),
            np.linspace(1.0, 5.0, 40),
        )
        figure = Figure()
        axes = figure.add_subplot(1, 1, 1)
        visible = {int(m) for m in scan.orders}
        colors = draw_efficiency(axes, scan, visible, THEMES["dark"])
        legend = axes.get_legend()
        assert legend.get_bbox_to_anchor().transformed(axes.transAxes.inverted()).x0 > 1.0
        lines = {line.get_label(): line.get_color() for line in axes.get_lines()}
        assert all(lines[f"m={m:+d}"] == colors[m] for m in visible)
        assert axes.get_facecolor()[:3] != (1.0, 1.0, 1.0)


@pytest.fixture
def light_again():
    """Every test that switches theme hands Light back, since the manager
    lives on the QApplication and outlasts the test."""
    yield
    from gratinglab.gui.qt import theming

    found = theming.manager()
    if found is not None:
        found.set_choice("light")


class TestSwitching:
    @pytest.fixture
    def win(self, qtbot, light_again):
        from gratinglab.gui.qt.main_window import MainWindow

        window = MainWindow()
        qtbot.addWidget(window)
        with qtbot.waitSignal(window.solved, timeout=60_000):
            pass
        return window

    def test_the_suite_runs_in_light(self, win):
        from gratinglab.gui.qt import theming

        assert theming.current() is LIGHT

    def test_switching_restyles_figures_pages_and_palette(self, win):
        from matplotlib.colors import to_hex

        from gratinglab.gui.qt import theming

        dark = THEMES["dark"]
        win.set_theme("dark")
        assert theming.current() is dark
        scalar = win.tabs["scalar"]
        assert to_hex(scalar._figure.get_facecolor()) == dark.figure
        assert to_hex(scalar._axes.get_facecolor()) == dark.axes
        from PySide6.QtWidgets import QApplication

        assert QApplication.palette().window().color().name() == dark.window
        assert dark.muted in scalar._provenance.toHtml()
        checked = [a.data() for a in win._theme_actions.actions() if a.isChecked()]
        assert checked == ["dark"]

    def test_order_list_shows_each_lines_colour(self, win):
        scalar = win.tabs["scalar"]
        items = [scalar._order_list.item(i) for i in range(scalar._order_list.count())]
        assert items and all(not item.icon().isNull() for item in items)

    def test_the_choice_is_remembered(self, win, tmp_path, monkeypatch):
        from PySide6.QtCore import QSettings

        from gratinglab.gui.qt import theming

        ini = tmp_path / "settings.ini"
        monkeypatch.setenv(theming.SETTINGS_ENV, str(ini))
        win.set_theme("solarized-dark")
        assert QSettings(str(ini), QSettings.Format.IniFormat).value(
            "appearance/theme"
        ) == "solarized-dark"

    def test_both_windows_offer_the_theme_and_stay_in_step(self, win, qtbot):
        """The AFM window had no theme menu; now both do, kept in step."""
        from matplotlib.colors import to_hex

        from gratinglab.metrology.gui.qt.main_window import MainWindow as AfmWindow

        afm = AfmWindow()
        qtbot.addWidget(afm)
        afm.load(SAMPLE)
        dark = THEMES["dark"]

        afm._theme_menu.actions_by_key()["dark"].trigger()

        assert [a.data() for a in win._theme_actions.actions() if a.isChecked()] == ["dark"]
        for figure in (afm.canvas.fig, afm.importer.figure, afm.boundary.figure):
            assert to_hex(figure.get_facecolor()) == dark.figure
            assert all(to_hex(ax.get_facecolor()) == dark.axes for ax in figure.axes)
        afm.close()

    def test_system_resolves_to_light_or_dark(self, win):
        from gratinglab.gui.qt import theming

        win.set_theme("system")
        assert theming.current().key in ("light", "dark")

    def test_an_unknown_theme_is_refused(self, win):
        from gratinglab.gui.qt import theming

        with pytest.raises(ValueError, match="unknown theme"):
            theming.manager().set_choice("sepia")
