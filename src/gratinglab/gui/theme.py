"""Colour themes for the desktop windows, as data.

Nothing here imports Qt. A theme is a set of named colours: the window
palette, the plot surfaces, the three verdict colours the provenance panel
uses, and the diagram's tag colours. The Qt layer (``gui.qt.theming``) turns
one into a ``QPalette`` and matplotlib settings; everything that draws asks
the active theme for a colour by meaning, never by value, which is the
arrangement ``provenance.TAG_COLORS`` and ``diagram.TAG_COLORS`` already
had.

``LIGHT`` is the look the windows had before themes existed, value for
value, so choosing it changes nothing but the order colours.

Order colours are by sign rather than from a cycle. The plots used
matplotlib's ten-colour cycle, so a scan with sixteen orders drew m = -3
and m = +7 in the same orange. Here negative orders take cool hues and
positive orders warm ones, each running away from m = 0 as |m| grows, with
alternate orders a little lighter and darker so neighbours separate; m = 0
is neutral. The colour of an order depends only on the set of orders in the
scan, so hiding one leaves every other colour where it was.
"""

from __future__ import annotations

import colorsys
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

__all__ = [
    "THEMES",
    "LIGHT",
    "DARK",
    "SYSTEM",
    "MIN_LINE_CONTRAST",
    "Theme",
    "contrast",
    "order_colors",
    "provenance_colors",
    "rc_params",
    "set_titles",
    "style_axes",
    "style_figure",
]

#: The choice that follows the operating system's light or dark setting.
SYSTEM = "system"


@dataclass(frozen=True)
class Theme:
    """One theme. Every field is a ``#rrggbb`` colour except the first three."""

    key: str
    label: str
    dark: bool
    # Window palette.
    window: str
    base: str
    alt_base: str
    text: str
    muted: str
    button: str
    highlight: str
    highlight_text: str
    link: str
    # Plot surfaces.
    figure: str
    axes: str
    ink: str
    grid: str
    # Verdicts, as provenance uses them.
    warn: str
    bad: str
    ok: str
    # Diagram tags, as ``diagram.TAG_COLORS`` names them.
    diagram: Mapping[str, str] = field(default_factory=dict)
    # Order colours: lightness and saturation of the hue ramps.
    order_lightness: float = 0.42
    order_saturation: float = 0.75


LIGHT = Theme(
    key="light", label="Light", dark=False,
    window="#efefef", base="#ffffff", alt_base="#f7f7f7", text="#1a1a1a",
    muted="#666666", button="#efefef", highlight="#2d5aa8",
    highlight_text="#ffffff", link="#2d5aa8",
    figure="#ffffff", axes="#ffffff", ink="#1a1a1a", grid="#b0b0b0",
    warn="#a5370d", bad="#b00020", ok="#0a7d28",
    diagram={
        "surface": "#1f3b73", "facet": "#2d5aa8", "incident": "#8fa3c8",
        "zero": "#666666", "order": "#4dc4ff", "blaze": "#ffc233",
        "normal": "#0a7d28", "cone": "#999999", "evanescent": "#b00020",
        "axis": "#888888",
    },
    order_lightness=0.46, order_saturation=0.90,
)

DARK = Theme(
    key="dark", label="Dark", dark=True,
    window="#1f2125", base="#17181b", alt_base="#222429", text="#e6e6e6",
    muted="#9aa0a6", button="#2a2d32", highlight="#2f6fd6",
    highlight_text="#ffffff", link="#7fb0ff",
    figure="#17181b", axes="#17181b", ink="#e6e6e6", grid="#4a4f55",
    warn="#ff9f5a", bad="#ff6b81", ok="#5ad17a",
    diagram={
        "surface": "#8fb0ff", "facet": "#b5c9ff", "incident": "#8fa3c8",
        "zero": "#b0b0b0", "order": "#4dc4ff", "blaze": "#ffc233",
        "normal": "#5ad17a", "cone": "#8a8f96", "evanescent": "#ff6b81",
        "axis": "#9aa0a6",
    },
    order_lightness=0.66, order_saturation=0.80,
)

# Solarized, Ethan Schoonover's palette: the base tones and accents are his.
_S = {
    "base03": "#002b36", "base02": "#073642", "base01": "#586e75",
    "base00": "#657b83", "base0": "#839496", "base1": "#93a1a1",
    "base2": "#eee8d5", "base3": "#fdf6e3", "yellow": "#b58900",
    "orange": "#cb4b16", "red": "#dc322f", "magenta": "#d33682",
    "violet": "#6c71c4", "blue": "#268bd2", "cyan": "#2aa198",
    "green": "#859900",
}
_SOLARIZED_DIAGRAM = {
    "surface": _S["blue"], "facet": _S["violet"], "incident": _S["base1"],
    "zero": _S["base00"], "order": _S["cyan"], "blaze": _S["yellow"],
    "normal": _S["green"], "cone": _S["base1"], "evanescent": _S["red"],
    "axis": _S["base1"],
}

SOLARIZED_LIGHT = Theme(
    key="solarized-light", label="Solarized Light", dark=False,
    window=_S["base2"], base=_S["base3"], alt_base="#f5efdc", text=_S["base02"],
    muted=_S["base01"], button=_S["base2"], highlight="#1f6fae",
    highlight_text=_S["base3"], link="#1f6fae",
    figure=_S["base3"], axes=_S["base3"], ink=_S["base02"], grid="#d9d2bf",
    warn="#b4410f", bad="#c42a27", ok="#5f6e00",
    diagram=_SOLARIZED_DIAGRAM,
    order_lightness=0.46, order_saturation=0.85,
)

SOLARIZED_DARK = Theme(
    key="solarized-dark", label="Solarized Dark", dark=True,
    window=_S["base02"], base=_S["base03"], alt_base="#04313d", text=_S["base2"],
    muted=_S["base1"], button=_S["base02"], highlight="#1f6fae",
    highlight_text=_S["base3"], link="#5fb3f0",
    figure=_S["base03"], axes=_S["base03"], ink=_S["base2"], grid="#1c4a56",
    warn="#e8743b", bad="#f06260", ok="#a4bf1d",
    diagram={**_SOLARIZED_DIAGRAM, "surface": "#5fb3f0", "facet": "#9a9ee8",
             "zero": _S["base1"], "evanescent": "#f06260"},
    order_lightness=0.62,
)

HIGH_CONTRAST = Theme(
    key="high-contrast", label="High contrast", dark=True,
    window="#000000", base="#1a1a1a", alt_base="#262626", text="#ffffff",
    muted="#a8a8a8", button="#2a2a2a", highlight="#ffff00",
    highlight_text="#000000", link="#00ffff",
    figure="#000000", axes="#000000", ink="#ffffff", grid="#5a5a5a",
    warn="#ffb000", bad="#ff5c5c", ok="#3cff6e",
    diagram={
        "surface": "#ffffff", "facet": "#00ffff", "incident": "#c0c0c0",
        "zero": "#ffffff", "order": "#00ffff", "blaze": "#ffff00",
        "normal": "#3cff6e", "cone": "#a0a0a0", "evanescent": "#ff5c5c",
        "axis": "#c0c0c0",
    },
    order_lightness=0.66, order_saturation=1.0,
)

#: Every theme by key, in menu order.
THEMES: dict[str, Theme] = {
    t.key: t for t in (LIGHT, DARK, SOLARIZED_LIGHT, SOLARIZED_DARK, HIGH_CONTRAST)
}

#: Hue ramps in degrees, from |m| = 1 outwards. Cool for negative orders,
#: warm for positive, so the sign is readable from the colour alone. The warm
#: ramp starts at orange, not yellow: on a light background the contrast floor
#: darkens yellow into olive and brown, and a run of browns is not legible.
_NEGATIVE_HUES = (185.0, 285.0)
_POSITIVE_HUES = (32.0, -48.0)
#: Alternate orders step this far in lightness, so adjacent ones separate.
_ALTERNATE_STEP = 0.10


def _hex(red: float, green: float, blue: float) -> str:
    return "#" + "".join(f"{round(255 * max(0.0, min(1.0, c))):02x}" for c in (red, green, blue))


def order_colors(orders: Iterable[int], theme: Theme) -> dict[int, str]:
    """A colour for every order in ``orders``, by sign and by rank in |m|.

    Pass every order the scan has, not only the visible ones, so that hiding
    an order leaves the others' colours alone.
    """
    unique = {int(m) for m in orders}
    colors: dict[int, str] = {}
    if 0 in unique:
        colors[0] = theme.muted
    for ranked, (start, stop) in (
        (sorted((m for m in unique if m < 0), reverse=True), _NEGATIVE_HUES),
        (sorted(m for m in unique if m > 0), _POSITIVE_HUES),
    ):
        count = len(ranked)
        for i, order in enumerate(ranked):
            hue = start + (stop - start) * (i / (count - 1) if count > 1 else 0.0)
            step = _ALTERNATE_STEP if i % 2 else -_ALTERNATE_STEP / 2
            lightness = theme.order_lightness + (step if theme.dark else -step)
            colors[order] = _legible(hue, lightness, theme)
    return colors


#: Lines must stand off the plot background by at least this WCAG ratio, the
#: 3:1 that WCAG 2.1 asks of graphical objects.
MIN_LINE_CONTRAST = 3.0


def _legible(hue: float, lightness: float, theme: Theme) -> str:
    """The colour at this hue, moved away from the background until legible.

    A fixed lightness cannot serve every hue: yellow at a lightness that
    suits blue on white is barely visible. So each colour starts at the
    ramp's lightness and steps darker (light theme) or lighter (dark theme)
    until it clears :data:`MIN_LINE_CONTRAST`.
    """
    direction = 1.0 if theme.dark else -1.0
    color = theme.ink
    for _ in range(40):
        color = _hex(
            *colorsys.hls_to_rgb((hue % 360.0) / 360.0, lightness, theme.order_saturation)
        )
        if contrast(color, theme.axes) >= MIN_LINE_CONTRAST:
            return color
        lightness = max(0.0, min(1.0, lightness + 0.02 * direction))
    return color


def provenance_colors(theme: Theme) -> dict[str, str]:
    """``provenance.TAG_COLORS`` for this theme."""
    return {"warn": theme.warn, "bad": theme.bad, "ok": theme.ok, "dim": theme.muted}


def _luminance(color: str) -> float:
    def channel(c: float) -> float:
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    value = color.lstrip("#")
    red, green, blue = (int(value[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
    return 0.2126 * channel(red) + 0.7152 * channel(green) + 0.0722 * channel(blue)


def contrast(foreground: str, background: str) -> float:
    """WCAG 2 contrast ratio between two ``#rrggbb`` colours, 1 to 21."""
    a, b = sorted((_luminance(foreground), _luminance(background)), reverse=True)
    return (a + 0.05) / (b + 0.05)


def rc_params(theme: Theme) -> dict[str, object]:
    """matplotlib settings for this theme, for ``rcParams.update``.

    They govern what is created after they are set. A figure that already
    exists also needs :func:`style_figure`, and each panel restyles its axes
    when it redraws.
    """
    return {
        "figure.facecolor": theme.figure,
        "figure.edgecolor": theme.figure,
        "savefig.facecolor": theme.figure,
        "savefig.edgecolor": theme.figure,
        "axes.facecolor": theme.axes,
        "axes.edgecolor": theme.ink,
        "axes.labelcolor": theme.ink,
        "axes.titlecolor": theme.ink,
        "text.color": theme.ink,
        "xtick.color": theme.ink,
        "ytick.color": theme.ink,
        "xtick.labelcolor": theme.ink,
        "ytick.labelcolor": theme.ink,
        "grid.color": theme.grid,
        "legend.facecolor": theme.axes,
        "legend.edgecolor": theme.grid,
        "legend.labelcolor": theme.ink,
        "axes3d.xaxis.panecolor": _pane(theme),
        "axes3d.yaxis.panecolor": _pane(theme),
        "axes3d.zaxis.panecolor": _pane(theme),
    }


def _pane(theme: Theme) -> tuple[float, float, float, float]:
    """A 3D pane tinted a little away from the axes colour."""
    value = theme.axes.lstrip("#")
    red, green, blue = (int(value[i:i + 2], 16) / 255.0 for i in (0, 2, 4))
    shift = 0.06 if theme.dark else -0.04
    return (*(max(0.0, min(1.0, c + shift)) for c in (red, green, blue)), 1.0)


def style_axes(axes, theme: Theme) -> None:
    """Restyle one existing axes, 2D or 3D, after a clear."""
    axes.set_facecolor(theme.axes)
    axes.tick_params(colors=theme.ink, labelcolor=theme.ink, which="both")
    for spine in getattr(axes, "spines", {}).values():
        spine.set_edgecolor(theme.ink)
    axes.xaxis.label.set_color(theme.ink)
    axes.yaxis.label.set_color(theme.ink)
    axes.title.set_color(theme.ink)
    zaxis = getattr(axes, "zaxis", None)
    if zaxis is not None:
        zaxis.label.set_color(theme.ink)
        for axis in (axes.xaxis, axes.yaxis, zaxis):
            axis.set_pane_color(_pane(theme))
            axis.line.set_color(theme.ink)


#: Plot titles: what the panel is, at this size and weight, left-aligned.
TITLE_SIZE = 11
TITLE_WEIGHT = "semibold"
#: The figures a panel reports beside its title, smaller and in ``muted``.
DETAIL_SIZE = 9


def set_titles(axes, title: str, detail: str | None = None, theme: Theme | None = None) -> None:
    """Title a panel: its name on the left, the numbers it reports on the right.

    One style for every plot, so a title says what the panel is and never has
    to carry a measurement, a setting or an explanation on a second line. The
    detail is where those go, quieter than the title. Pass ``theme`` from a
    caller that has one; otherwise the plain ``LIGHT`` colours are used.
    """
    theme = theme or LIGHT
    axes.set_title(title, loc="left", fontsize=TITLE_SIZE, fontweight=TITLE_WEIGHT,
                   color=theme.ink)
    axes.set_title(detail or "", loc="right", fontsize=DETAIL_SIZE, color=theme.muted)


def style_figure(figure, theme: Theme) -> None:
    """Restyle an existing figure and every axes on it."""
    figure.set_facecolor(theme.figure)
    for axes in figure.axes:
        style_axes(axes, theme)
