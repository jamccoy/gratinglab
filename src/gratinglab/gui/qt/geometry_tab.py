"""The Grating Geometry tab: what the grating and the beams actually look like.

Thin, like every widget here. :mod:`gratinglab.gui.diagram` decides what to
draw -- every ray direction, every cutoff, every colour tag, every caption --
and this class only positions the result. It computes no angle, makes no
comparison against 1.0, and picks no colour, so a bug in this file misplaces a
line rather than drawing wrong physics.

Not a solver tab: no ``name`` attribute, no ``solve_requested``, no
``build_options``. `MainWindow._solve_active_tab` reads ``getattr(widget,
"name", None)`` and therefore skips it, the same way it skips
:class:`~.setup_tab.SetupTab`. Nothing here needs a solver -- a groove's shape
and the directions light leaves in are geometry, and geometry is already known
the moment the form parses.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
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
    "ladder": "sin β  —  every order, propagating or not",
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

        See `_redraw`: the figure is the expensive part and there is no point
        rendering it into a tab nobody is looking at.
        """
        super().showEvent(event)
        if self._needs_paint and self._diagram is not None:
            self._paint(self._diagram)
            self._needs_paint = False

    # -- construction ------------------------------------------------

    def _build(self) -> None:
        """Hero on the left, supporting panels on the right, controls beneath.

        Controls sit in a bottom *strip* rather than the right *column* they
        had before. Once the profile plot moved in here the tab became
        width-constrained rather than height-constrained, and a 320 px column
        was spending the scarce dimension; a ~170 px strip spends the one we
        now have to spare, handing that width back to the plots.
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
        """The 3D conical view -- the one thing a dispersion-plane drawing
        cannot show."""
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
        row.addWidget(self._build_wavelength_group(), 2)
        row.addWidget(self._build_caption_group(), 3)
        return strip

    def _build_wavelength_group(self) -> QWidget:
        group = QGroupBox("Wavelength")
        layout = QVBoxLayout(group)

        self._wavelength_label = QLabel("—")
        layout.addWidget(self._wavelength_label)

        self._slider = QSlider(Qt.Orientation.Horizontal)
        self._slider.setMinimum(0)
        self._slider.setMaximum(0)
        self._slider.valueChanged.connect(self._on_slider_moved)
        layout.addWidget(self._slider)

        self._blaze_button = QPushButton("Jump to blaze λ")
        self._blaze_button.clicked.connect(self._jump_to_blaze)
        row = QHBoxLayout()
        row.addWidget(self._blaze_button)
        layout.addLayout(row)
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
        """Draw a parsed geometry. Needs no solver and no solve."""
        self._parsed = parsed
        count = len(parsed.wavelengths)

        # Keep the slider's position across a redraw where the grid is the
        # same size; a wavelength that jumped on every keystroke would make
        # the picture unreadable while typing.
        index = self._slider.value() if self._slider.maximum() == count - 1 else count // 2

        self._slider.blockSignals(True)
        self._slider.setMaximum(max(0, count - 1))
        self._slider.setValue(min(index, count - 1))
        self._slider.blockSignals(False)

        self._sync_blaze_button()
        self.profile_panel.draw(parsed)
        self._redraw()

    def show_pending(self, errors) -> None:
        """The form does not parse yet. Keep the last good drawing.

        A form mid-edit is not a form with an error: ``"2."`` and ``"-"`` each
        fail to parse for exactly one keystroke, and blanking the canvas or
        reddening a label for them would style ordinary typing as a problem.
        So the picture stays, and a `dim` line says plainly that it is one
        edit behind -- never `warn`, never `bad`, which are what the panel
        reserves for something actually wrong.

        Errors are counted, not listed: naming a half-typed field would be the
        same accusation in smaller type.
        """
        if self._parsed is None:
            return  # nothing good to keep; the tab has never drawn
        count = len(errors)
        note = provenance.Line(
            # Trailing newline: `to_html` concatenates runs, so a line that
            # does not end one runs into the caption beneath it.
            f"waiting for {count} field{'' if count == 1 else 's'} to finish — "
            "the drawing below is from the last complete geometry\n",
            "dim",
        )
        self._captions.setHtml(
            provenance.to_html(
                (note,) + self._scene.captions + self._diagram.captions
            )
        )

    # -- reactions ----------------------------------------------------

    def _on_slider_moved(self, _value: int) -> None:
        """Moving the slider redraws. It never solves -- the geometry is
        already known, and nothing here depends on a solver."""
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
        """Enabled when it can act; disabled *with its reason* when it cannot.

        A disabled control that explains itself is honest. An enabled one that
        silently does nothing is the mistake this project already made once
        and named -- see `docs/theory/scalar.md` §5.
        """
        index, reason = diagram_module.blaze_jump(
            self._parsed.wavelengths, self._parsed.problem, self._parsed.illumination
        )
        self._blaze_button.setEnabled(index is not None)
        self._blaze_button.setToolTip(reason)

    # -- drawing ------------------------------------------------------

    def _redraw(self) -> None:
        """Rebuild the diagram; paint it only if anyone can see it.

        Building is ~1 ms of numpy and always runs, so `_diagram` is always
        current and a caller (or a test) can read it without a window. Painting
        is ~50 ms of matplotlib, so it waits for `showEvent` when the tab is
        hidden -- rendering a figure into a tab nobody is looking at is time
        spent for no one.
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
        self._captions.setHtml(
            provenance.to_html(self._scene.captions + self._diagram.captions)
        )

        if self.isVisible():
            self._paint(self._diagram)
            self._needs_paint = False
        else:
            self._needs_paint = True

    def _paint(self, drawing) -> None:
        """Position what `diagram` decided. No trig, no cutoffs, no colours."""
        colors = diagram_module.TAG_COLORS
        self._scene3d.show_scene(self._scene)

        for name, axes in self._axes.items():
            axes.clear()
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
                axes.annotate(
                    arrow.label,
                    xy=(arrow.x1, arrow.y1),
                    fontsize=7,
                    color=colors[arrow.tag],
                    ha="center",
                    va="bottom",
                )

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

        self._canvas.draw_idle()
