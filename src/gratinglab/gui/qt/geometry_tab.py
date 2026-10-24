"""The Grating Geometry tab: what the grating and the beams look like.

Like every widget here, this tab is thin. :mod:`gratinglab.gui.diagram`
decides what to draw (every ray direction, cutoff, colour tag and caption),
and this class only positions the result. It computes no angle, makes no
comparison against 1.0 and picks no colour, so a bug in this file can
misplace a line but cannot draw wrong physics.

This is not a solver tab: it has no ``name`` attribute, no
``solve_requested`` and no ``build_options``. `MainWindow._solve_active_tab`
reads ``getattr(widget,
"name", None)`` and therefore skips it, in the same way that it skips
:class:`~.setup_tab.SetupTab`. Nothing here needs a solver, since the shape of
a groove and the directions in which light leaves it are geometry, and the
geometry is known as soon as the form parses.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QPalette
from PySide6.QtWidgets import (
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .. import diagram as diagram_module
from .. import diagram3d as diagram3d_module
from .. import provenance
from ..theme import provenance_colors, style_axes, style_figure
from . import theming
from .profile_plot_panel import ProfilePlotPanel
from .scene3d_panel import Scene3DPanel

__all__ = ["GeometryTab"]

#: Named glyph -> matplotlib marker. The diagram module names what a point
#: *is*; this is the only place that decides what it looks like.
_GLYPHS = {
    "out-of-page": "o",
    "strike": "o",
}

_PANEL_TITLES = {
    "main": "Down the groove axis ĝ",
    "ladder": "sin β: every order, propagating or not",
}


class GeometryTab(QWidget):
    """The grating, the beams, and what the drawing cannot show."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._parsed = None
        self._diagram = None
        self._scene = None
        self._needs_paint = False
        self._build()

    def showEvent(self, event) -> None:  # noqa: N802 - Qt's spelling
        """Paint on becoming visible, if a redraw happened while hidden.

        See `_redraw`: the figure is the expensive part, and rendering it into
        a tab that is not visible would serve no purpose.
        """
        super().showEvent(event)
        if self._needs_paint and self._diagram is not None:
            self._paint(self._diagram)
            self._needs_paint = False

    # -- construction ------------------------------------------------

    def _build(self) -> None:
        """Lay out the main view, the supporting panels and the controls.

        The main view is on the left, the supporting panels on the right and
        the controls beneath. The controls sit in a strip along the bottom
        instead of the column on the right that they occupied before. Once
        the profile plot moved into this tab, the tab became limited by width
        instead of height, and a 320 px column was using the scarce dimension.
        A ~170 px strip uses the dimension that is now available and returns
        the width to the plots.
        """
        hero = self._build_hero()
        self.profile_panel = ProfilePlotPanel()

        side = QSplitter(Qt.Orientation.Vertical)
        side.addWidget(self._build_plots())
        side.addWidget(self.profile_panel)
        side.setStretchFactor(0, 1)
        side.setStretchFactor(1, 1)
        side.setSizes([360, 229])
        side.setMinimumWidth(280)
        self.profile_panel.setMinimumHeight(150)

        body = QSplitter(Qt.Orientation.Horizontal)
        body.addWidget(hero)
        body.addWidget(side)
        body.setStretchFactor(0, 2)
        body.setStretchFactor(1, 1)
        body.setSizes([578, 290])

        outer = QSplitter(Qt.Orientation.Vertical)
        outer.addWidget(body)
        outer.addWidget(self._build_controls())
        outer.setStretchFactor(0, 1)
        outer.setStretchFactor(1, 0)
        outer.setSizes([593, 170])

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(outer)

    def _build_hero(self) -> QWidget:
        """The 3D conical view, which shows what a dispersion-plane drawing
        cannot."""
        self._scene3d = Scene3DPanel()
        self._scene3d.setMinimumWidth(380)
        self._scene3d.setMinimumHeight(320)
        return self._scene3d

    def _build_plots(self) -> QWidget:
        """The 2D panels, beside the 3D one: the dispersion-plane
        cross-section (where the groove facets are legible) and the sin β
        ladder (where an evanescent order still appears)."""
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg,
            NavigationToolbar2QT,
        )
        from matplotlib.figure import Figure

        # Two panels since M13-I: the γ sliver retired when the 3D view began
        # drawing the cone itself. Its row goes to the cross-section, which is
        # the one that has angles to be read.
        self._figure = Figure(figsize=(4, 6), layout="constrained")
        grid = self._figure.add_gridspec(2, 1, height_ratios=[3.0, 1.0])
        self._axes = {
            "main": self._figure.add_subplot(grid[0, 0]),
            "ladder": self._figure.add_subplot(grid[1, 0]),
        }
        self._canvas = FigureCanvasQTAgg(self._figure)
        style_figure(self._figure, theming.current())
        theming.subscribe(self._on_theme)
        self._caption_lines: tuple = ()

        panel = QWidget()
        panel.setMinimumHeight(200)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(NavigationToolbar2QT(self._canvas, panel))
        layout.addWidget(self._canvas)
        return panel

    def _build_controls(self) -> QWidget:
        strip = QWidget()
        strip.setMinimumHeight(150)
        strip.setMaximumHeight(240)
        row = QHBoxLayout(strip)
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self._build_wavelength_group(), 1)
        row.addWidget(self._build_caption_group(), 3)
        return strip

    def _build_wavelength_group(self) -> QWidget:
        """The wavelength being drawn, chosen from the scan.

        Laid out as two tight rows -- the value and the blaze jump, then the
        slider between the scan's two ends -- rather than three stacked
        widgets with the box's spare height between them, which left most of
        the box empty and narrowed the captions beside it.
        """
        group = QGroupBox("Wavelength")
        layout = QVBoxLayout(group)

        self._wavelength_label = QLabel("—")
        bigger = self._wavelength_label.font()
        bigger.setPointSizeF(bigger.pointSizeF() * 1.25)
        bigger.setBold(True)
        self._wavelength_label.setFont(bigger)
        self._blaze_button = QPushButton("Jump to blaze λ")
        self._blaze_button.clicked.connect(self._jump_to_blaze)
        top = QHBoxLayout()
        top.addWidget(self._wavelength_label)
        top.addStretch(1)
        top.addWidget(self._blaze_button)
        layout.addLayout(top)

        self._slider = QSlider(Qt.Orientation.Horizontal)
        self._slider.setMinimum(0)
        self._slider.setMaximum(0)
        self._slider.valueChanged.connect(self._on_slider_moved)
        self._scan_start = QLabel("")
        self._scan_stop = QLabel("")
        for end in (self._scan_start, self._scan_stop):
            end.setForegroundRole(QPalette.ColorRole.Mid)
        bottom = QHBoxLayout()
        bottom.addWidget(self._scan_start)
        bottom.addWidget(self._slider, 1)
        bottom.addWidget(self._scan_stop)
        layout.addLayout(bottom)
        layout.addStretch(1)
        return group

    def _build_caption_group(self) -> QWidget:
        group = QGroupBox("What this shows")
        layout = QVBoxLayout(group)
        self._captions = QTextEdit()
        self._captions.setReadOnly(True)
        layout.addWidget(self._captions)
        return group

    # -- the contract MainWindow uses --------------------------------

    def show_geometry(self, parsed) -> None:
        """Draw a parsed geometry, which needs no solver and no solve."""
        self._parsed = parsed
        count = len(parsed.wavelengths)

        # Keep the slider's position across a redraw where the grid is the
        # same size; a wavelength that jumped on every keystroke would make
        # the picture unreadable while typing.
        index = self._slider.value() if self._slider.maximum() == count - 1 else count // 2

        self._scan_start.setText(f"{parsed.wavelengths[0]:.4g}")
        self._scan_stop.setText(f"{parsed.wavelengths[-1]:.4g} nm")

        self._slider.blockSignals(True)
        self._slider.setMaximum(max(0, count - 1))
        self._slider.setValue(min(index, count - 1))
        self._slider.blockSignals(False)

        self._sync_blaze_button()
        self.profile_panel.draw(parsed)
        self._redraw()

    def show_pending(self, errors) -> None:
        """Keep the last good drawing while the form does not yet parse.

        A form in the middle of an edit is not a form with an error: ``"2."``
        and ``"-"`` each fail to parse for exactly one keystroke, and blanking
        the canvas or turning a label red for them would present ordinary
        typing as a problem. The picture therefore stays, and a `dim` line
        states that it is one edit behind. The line is never `warn` or `bad`,
        which the panel reserves for something that is wrong.

        The errors are counted but not listed, because naming a half-typed
        field would also present it as an error.
        """
        if self._parsed is None:
            return  # nothing good to keep; the tab has never drawn
        count = len(errors)
        note = provenance.Line(
            # Trailing newline: `to_html` concatenates runs, so a line that
            # does not end one runs into the caption beneath it.
            f"waiting for {count} field{'' if count == 1 else 's'} to finish; "
            "the drawing below is from the last complete geometry\n",
            "dim",
        )
        self._set_captions((note,) + self._scene.captions + self._diagram.captions)

    # -- reactions ----------------------------------------------------

    def _on_slider_moved(self, _value: int) -> None:
        """Redraw when the slider moves.

        It never solves, because the geometry is already known and nothing
        here depends on a solver."""
        self._redraw()

    def _jump_to_blaze(self) -> None:
        if self._parsed is None:
            return
        index, _ = diagram_module.blaze_jump(
            self._parsed.wavelengths, self._parsed.problem, self._parsed.illumination
        )
        if index is not None:
            self._slider.setValue(index)

    def _sync_blaze_button(self) -> None:
        """Enable the button when it can act, and otherwise disable it with its reason.

        A disabled control that states why it is disabled is preferable to an
        enabled one that silently does nothing; this project has made and
        documented the latter mistake once before (see `docs/theory/scalar.md`
        §5).
        """
        index, reason = diagram_module.blaze_jump(
            self._parsed.wavelengths, self._parsed.problem, self._parsed.illumination
        )
        self._blaze_button.setEnabled(index is not None)
        self._blaze_button.setToolTip(reason)

    # -- drawing ------------------------------------------------------

    def _redraw(self) -> None:
        """Rebuild the diagram; paint it only if anyone can see it.

        Building takes ~1 ms of numpy and always runs, so `_diagram` is always
        current and a caller (or a test) can read it without a window.
        Painting takes ~50 ms of matplotlib, so when the tab is hidden it
        waits for `showEvent`, since rendering into a tab that is not visible
        would be wasted time.
        """
        if self._parsed is None:
            return
        wavelength = float(self._parsed.wavelengths[self._slider.value()])
        self._wavelength_label.setText(f"λ = {wavelength:.4g} nm")

        self._diagram = diagram_module.build(
            self._parsed.problem, self._parsed.illumination, wavelength
        )
        self._scene = diagram3d_module.build_scene(
            self._parsed.problem, self._parsed.illumination, wavelength
        )
        # Cheap, and it is the part a reader needs even before the picture
        # arrives, so it is never deferred. The 3D scene's captions carry what
        # the projection cannot show, so both sets appear.
        self._set_captions(self._scene.captions + self._diagram.captions)

        if self.isVisible():
            self._paint(self._diagram)
            self._needs_paint = False
        else:
            self._needs_paint = True

    def _set_captions(self, lines) -> None:
        self._caption_lines = tuple(lines)
        self._captions.setHtml(
            provenance.to_html(self._caption_lines, provenance_colors(theming.current()))
        )

    def _on_theme(self, theme) -> None:
        style_figure(self._figure, theme)
        if self._caption_lines:
            self._set_captions(self._caption_lines)
        if self._diagram is None:
            self._canvas.draw_idle()
        elif self.isVisible():
            self._paint(self._diagram)
        else:
            self._needs_paint = True

    def _paint(self, drawing) -> None:
        """Position what `diagram` decided, choosing no angle, cutoff or colour."""
        theme = theming.current()
        colors = {**diagram_module.TAG_COLORS, **theme.diagram}
        self._scene3d.show_scene(self._scene)

        for name, axes in self._axes.items():
            axes.clear()
            style_axes(axes, theme)
            axes.set_title(_PANEL_TITLES[name], fontsize=9)

        for path in drawing.paths:
            axes = self._axes[path.panel]
            axes.plot(
                path.x, path.y,
                color=colors[path.tag],
                lw=1.6,
                ls="--" if path.dashed else "-",
                label=path.label,
            )
            if path.fill_to is not None:
                axes.fill_between(path.x, path.fill_to, path.y,
                                  color=colors[path.tag], alpha=0.12)

        labels: list = []
        for arrow in drawing.arrows:
            axes = self._axes[arrow.panel]
            axes.annotate(
                "",
                xy=(arrow.x1, arrow.y1),
                xytext=(arrow.x0, arrow.y0),
                arrowprops={
                    "arrowstyle": "-|>",
                    "color": colors[arrow.tag],
                    "lw": 1.5,
                    "ls": "--" if arrow.dashed else "-",
                    "shrinkA": 0,
                    "shrinkB": 0,
                },
            )
            if arrow.label:
                # Placed just beyond the tip along the arrow's own direction,
                # then pushed further out if it collides (`_declutter`).
                dx, dy = arrow.x1 - arrow.x0, arrow.y1 - arrow.y0
                length = float(np.hypot(dx, dy)) or 1.0
                direction = (dx / length, dy / length)
                text = axes.annotate(
                    arrow.label,
                    xy=(arrow.x1, arrow.y1),
                    xytext=(_LABEL_GAP * direction[0], _LABEL_GAP * direction[1]),
                    textcoords="offset points",
                    fontsize=7,
                    color=colors[arrow.tag],
                    ha="center",
                    va="center",
                )
                labels.append((axes, text, direction))

        for marker in drawing.markers:
            axes = self._axes[marker.panel]
            axes.plot(
                marker.x, marker.y,
                _GLYPHS[marker.glyph],
                color=colors[marker.tag],
                ms=4 if marker.glyph == "out-of-page" else 5,
                mfc="none" if marker.glyph == "out-of-page" else colors[marker.tag],
            )
            if marker.label:
                axes.annotate(
                    marker.label,
                    xy=(marker.x, marker.y),
                    xytext=(0, 7),
                    textcoords="offset points",
                    fontsize=6,
                    color=colors[marker.tag],
                    ha="center",
                )

        for name, ((x0, x1), (y0, y1)) in drawing.limits.items():
            axes = self._axes[name]
            axes.set_xlim(x0, x1)
            axes.set_ylim(y0, y1)

        # Equal aspect on the panel whose angles are meant to be true. The
        # ladder is a number line and has no aspect to preserve.
        self._axes["main"].set_aspect("equal", adjustable="box")
        self._axes["main"].set_xlabel("position along d̂ (nm)", fontsize=8)
        self._axes["main"].set_ylabel("height (nm)", fontsize=8)
        self._axes["ladder"].set_yticks([])
        for axes in self._axes.values():
            axes.tick_params(labelsize=7)
        self._axes["main"].grid(alpha=0.15)

        _declutter(self._canvas, labels)
        self._canvas.draw_idle()


#: Points between an arrow's tip and its label, and the step a colliding
#: label is pushed out by, along its own arrow.
_LABEL_GAP = 7.0
_LABEL_STEP = 5.0
_LABEL_MAX = 45.0


def _declutter(canvas, labels) -> None:
    """Push colliding arrow labels outward until each one is clear.

    The order rays fan out from one strike point and their tips sit a few
    degrees apart on the same arc, so labels placed at the tips stacked on
    one another. A fixed alternation cannot know which pairs collide, so
    this measures: labels are taken in fan order (by angle) and each one is
    pushed further along its own arrow until its box clears every label
    already placed and the panel title, up to a limit. Moving along the
    arrow keeps every label pointing at the ray it names.
    """
    if not labels:
        return
    renderer = canvas.get_renderer()
    placed: dict[object, list] = {}
    for axes, text, (ux, uy) in sorted(
        labels, key=lambda item: np.arctan2(item[2][1], item[2][0])
    ):
        blockers = placed.setdefault(axes, [axes.title.get_window_extent(renderer)])
        reach = _LABEL_GAP
        while True:
            box = text.get_window_extent(renderer).expanded(1.08, 1.15)
            if reach >= _LABEL_MAX or not any(box.overlaps(b) for b in blockers):
                break
            reach += _LABEL_STEP
            text.set_position((reach * ux, reach * uy))
        blockers.append(text.get_window_extent(renderer))
