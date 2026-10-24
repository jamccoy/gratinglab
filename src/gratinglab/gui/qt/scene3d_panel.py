"""The 3D conical-diffraction view.

Like every widget here, this panel is thin. :mod:`gratinglab.gui.diagram3d`
decides every direction, surface and caption, and this panel positions them.
It computes no angle and picks no colour.

These axes have no navigation toolbar. The zoom of `NavigationToolbar2QT`
rescales the data limits of an `Axes3D`, which breaks the cubic bounding box
that `set_box_aspect((1, 1, 1))` depends on, and the angles this view is meant
to show would then silently become incorrect. A control that can make the
picture misleading must not be offered; the same rule led to the removal of
the obliquity factor. The 2D figures keep their toolbars, because there pan
and zoom cannot distort an angle that the figure claims to show true.

Rotation needs no wiring: `Axes3D.__init__` calls `mouse_init()` and connects
its own mouse handlers.
"""

from __future__ import annotations

import numpy as np
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .. import diagram3d
from ..theme import style_axes, style_figure
from . import theming

#: How far the cube is zoomed inside its canvas. matplotlib's default leaves
#: a wide margin around a 3D box; this fills the panel without clipping the
#: tick labels at the default window size.
_BOX_ZOOM = 1.2

#: Where each kind of label sits around its anchor. At grazing incidence the
#: cone-axis tip, the m=0 tip and the blaze direction's tip are all within a
#: degree or two of one point, and 3D text cannot move out of its
#: neighbours' way, so each is put in a different quadrant around the point
#: instead: the geometry is not moved to make room.
_LABEL_ANCHORS = {
    "cone": ("right", "bottom"),
    "zero": ("left", "top"),
    "blaze": ("left", "bottom"),
}

#: How far along its ray a label sits, as a fraction of the ray. The cone
#: axis is labelled part-way along its shaft rather than at its tip, which
#: is exactly where the order tips cluster; inside 1.0 also keeps the label
#: within the frame, which hides anything beyond the tips.
_LABEL_REACH = {"cone": 0.72}

__all__ = ["Scene3DPanel"]

#: Named glyph -> matplotlib marker. The scene names what a point *is*; this is
#: the only place that decides what it looks like.
_GLYPHS = {"out-of-page": "o", "strike": "o"}


def _display(v) -> tuple:
    r"""Physical :math:`(\hat{d}, \hat{n}, -\hat{g})` to the frame matplotlib
    draws.

    `Axes3D` always draws its own z axis vertically, while the axis a reader
    expects to point up is the grating normal, since a grating lies flat. The
    scene is in the frame of conventions.md §3, where :math:`\hat{n}` is the
    second component, so the widget permutes to
    :math:`(\hat{d}, \hat{g}, \hat{n})`.

    This ordering is used instead of :math:`(\hat{d}, -\hat{g}, \hat{n})`
    because :math:`\hat{d} \times \hat{g} = \hat{n}`; the alternative is
    left-handed, and a left-handed set would mirror the whole picture.

    The permutation is a choice of view with no physical content, since no
    length or angle changes under a permutation of orthonormal axes, and for
    that reason it lives here and not in the scene.
    """
    x, y, z = v
    return x, -z, y


def _view_angles(physical) -> tuple[float, float]:
    """A physical look-direction to matplotlib's ``(elev, azim)``, degrees.

    The inverse of :func:`gratinglab.gui.diagram3d.view_direction`, after the
    display permutation. It concerns the view and not the physics, so it
    lives here. The round trip through `view_direction` can be asserted
    directly, and a test uses it to show that each preset points where its
    name says.
    """
    dx, dy, dz = _display(np.asarray(physical, dtype=float))
    return float(np.degrees(np.arcsin(-dz))), float(np.degrees(np.arctan2(-dy, -dx)))


class Scene3DPanel(QWidget):
    """One `Axes3D`, with the controls needed to keep its angles true."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._scene = None
        self._focus = "whole"  # "whole" | "rim"
        self._lod = "fine"
        self._artists: dict[str, list] = {"fine": [], "coarse": []}
        self._bounds: dict = {}
        self._lod_of: dict = {}
        self._build()

    def _build(self) -> None:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(6, 5), layout="constrained")
        self._axes = self._figure.add_subplot(projection="3d")
        self._canvas = FigureCanvasQTAgg(self._figure)
        style_figure(self._figure, theming.current())
        theming.subscribe(self._on_theme)

        # Rotation needs no wiring -- Axes3D connects its own handlers -- but
        # the level-of-detail swap does.
        self._canvas.mpl_connect("button_press_event", self._on_press)
        self._canvas.mpl_connect("button_release_event", self._on_release)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._canvas)
        layout.addWidget(self._build_controls())

        self._apply_preset("Oblique")

    def _build_controls(self) -> QWidget:
        strip = QWidget()
        row = QHBoxLayout(strip)
        row.setContentsMargins(4, 0, 4, 0)

        self._preset_box = QComboBox()
        self._preset_box.addItems(list(diagram3d.PRESET_VIEWS))
        self._preset_box.currentTextChanged.connect(self._apply_preset)
        row.addWidget(QLabel("View:"))
        row.addWidget(self._preset_box)

        self._focus_button = QPushButton("Focus on cone rim")
        self._focus_button.setCheckable(True)
        self._focus_button.toggled.connect(self._on_focus_toggled)
        row.addWidget(self._focus_button)

        self._magnification = QLabel("")
        self._magnification.setStyleSheet("color: #666")
        row.addWidget(self._magnification)

        row.addStretch(1)
        save = QPushButton("Save image…")
        save.clicked.connect(self.save_image)
        row.addWidget(save)
        return strip

    # -- reactions ----------------------------------------------------

    def _apply_preset(self, name: str) -> None:
        """Point the camera along a preset's physical direction."""
        if name not in diagram3d.PRESET_VIEWS:
            return
        elev, azim = _view_angles(diagram3d.PRESET_VIEWS[name])
        self._axes.view_init(elev=elev, azim=azim)

        # The rim preset exists *for* the rim; separating the two would make
        # the one useful combination a two-step discovery.
        if name == diagram3d.RIM_PRESET and self._focus_button.isEnabled():
            self._focus_button.setChecked(True)

        self._canvas.draw_idle()

    def _on_focus_toggled(self, checked: bool) -> None:
        self._focus = "rim" if checked else "whole"
        if self._scene is not None:
            self._apply_limits(self._scene)
            self._canvas.draw_idle()

    def _on_press(self, event) -> None:
        """Coarsen the meshes for the duration of a rotate drag.

        `get_navigate_mode() is None` is the exact predicate matplotlib itself
        uses at axes3d.py:1787 to decide whether a drag will rotate. There is
        no toolbar on these axes, so it is always true, but the test documents
        the coupling and remains correct if a toolbar is added later.
        """
        if (
            event.inaxes is self._axes
            and event.button == 1
            and self._axes.get_navigate_mode() is None
        ):
            self._set_lod("coarse")

    def _on_release(self, _event) -> None:
        self._set_lod("fine")  # unconditional: fine is always the safe direction

    def _set_lod(self, lod: str) -> None:
        if lod == self._lod or self._scene is None:
            return
        self._lod = lod
        # Through _apply_limits, so the frame filter and the level-of-detail
        # choice cannot disagree about what is visible.
        self._apply_limits(self._scene)
        self._canvas.draw_idle()

    def save_image(self) -> None:
        """Save the figure, the one toolbar function worth offering on its own."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save 3D view", "geometry.png", "PNG (*.png);;PDF (*.pdf)"
        )
        if not path:
            return
        self._figure.savefig(path, dpi=200)
        QMessageBox.information(self, "Saved", f"Wrote {path}")

    # -- the contract GeometryTab uses --------------------------------

    def _on_theme(self, theme) -> None:
        style_figure(self._figure, theme)
        if self._scene is not None:
            self._paint(self._scene)
        else:
            self._canvas.draw_idle()

    def show_scene(self, scene) -> None:
        """Draw a scene.

        The panel keeps the scene, so a caller can read back what was drawn
        without re-deriving it."""
        self._scene = scene
        self._sync_focus_button(scene)
        self._paint(scene)

    def _sync_focus_button(self, scene) -> None:
        """Enable the button when it can act, and otherwise disable it with its reason.

        At γ = 90° the rim is already the whole scene, so focusing on it would
        zoom out. A disabled control that states why it is disabled is
        preferable to an enabled one that silently does nothing.
        """
        self._focus_button.setEnabled(scene.rim_available)
        self._focus_button.setToolTip(scene.rim_reason)
        if not scene.rim_available and self._focus_button.isChecked():
            self._focus_button.setChecked(False)

    def _apply_limits(self, scene) -> None:
        """Frame either the whole scene or the cone rim.

        Both are cubes, and both get ``set_box_aspect((1, 1, 1))``, so moving
        between them is a uniform scale plus a translation, i.e. a similarity
        transform. Every angle in the picture is preserved exactly, and only
        the part of the geometry that fills the frame changes. For this reason
        the magnification can be stated without any disclaimer.
        """
        use_rim = self._focus == "rim" and scene.rim_available
        (x0, x1), (y0, y1), (z0, z1) = scene.rim_limits if use_rim else scene.limits

        # Limits permute with the coordinates, exactly as the points do.
        self._axes.set_xlim(x0, x1)
        self._axes.set_ylim(-z1, -z0)
        self._axes.set_zlim(y0, y1)
        # A cube, so no axis is stretched relative to another -- without this
        # every angle in the picture is a lie.
        # The zoom scales the whole box uniformly, so it preserves every
        # angle just as the cube does; it only stops the default margin from
        # leaving the scene a small cube in a large canvas.
        self._axes.set_box_aspect((1.0, 1.0, 1.0), zoom=_BOX_ZOOM)

        self._magnification.setText(
            f"cone rim, {scene.rim_magnification:.0f}×" if use_rim else ""
        )
        # Display ranges, not the physical ones: `_bounds` holds permuted
        # coordinates, and comparing the two frames would hide everything.
        self._hide_artists_outside_the_frame((x0, x1), (-z1, -z0), (y0, y1))

    def _hide_artists_outside_the_frame(self, xr, yr, zr) -> None:
        """Hide primitives with no point inside the current box.

        `Axes3D` does not clip its artists to the axes limits, so at 27x the
        grating patch, which lies entirely outside the rim frame, would be
        smeared across the edges instead of disappearing. This is a purely
        geometric test on coordinates that are already computed; it hides
        nothing that would have been visible and decides no physics.
        """
        for artist, bounds in self._bounds.items():
            (ax0, ax1), (ay0, ay1), (az0, az1) = bounds
            inside = (
                ax1 >= min(xr) and ax0 <= max(xr)
                and ay1 >= min(yr) and ay0 <= max(yr)
                and az1 >= min(zr) and az0 <= max(zr)
            )
            lod_ok = self._lod_of.get(artist, self._lod) == self._lod
            artist.set_visible(inside and lod_ok)

    # -- drawing ------------------------------------------------------

    def _paint(self, scene) -> None:
        theme = theming.current()
        colors = {**diagram3d.TAG_COLORS, **theme.diagram}
        axes = self._axes
        # view_init survives cla(); preserve it explicitly anyway so a redraw
        # never yanks the camera out from under a user mid-inspection.
        elev, azim = axes.elev, axes.azim
        axes.clear()
        style_axes(axes, theme)

        # Both levels of detail are built once and swapped by visibility, so a
        # drag never resamples anything -- see `_set_lod`.
        self._artists = {"fine": [], "coarse": []}
        self._bounds = {}
        self._lod_of = {}

        def record(artist, xs, ys, zs, lod=None):
            """Remember an artist's extent so the frame filter can reach it."""
            self._bounds[artist] = (
                (float(np.min(xs)), float(np.max(xs))),
                (float(np.min(ys)), float(np.max(ys))),
                (float(np.min(zs)), float(np.max(zs))),
            )
            if lod is not None:
                self._lod_of[artist] = lod
                self._artists[lod].append(artist)

        for surface in scene.surfaces:
            sx, sy, sz = _display((surface.x, surface.y, surface.z))
            artist = axes.plot_surface(
                sx, sy, sz,
                color=colors[surface.tag],
                alpha=surface.alpha,
                linewidth=0,
                antialiased=False,
                shade=False,
            )
            record(artist, sx, sy, sz, surface.lod)

        for curve in scene.curves:
            cx, cy, cz = _display((curve.x, curve.y, curve.z))
            (artist,) = axes.plot(
                cx, cy, cz,
                color=colors[curve.tag],
                lw=1.0,
                ls="--" if curve.dashed else "-",
                alpha=0.7,
            )
            record(artist, cx, cy, cz, curve.lod)

        for ray in scene.rays:
            tail = _display(ray.origin)
            head = _display(ray.head)
            xs = [tail[0], head[0]]
            ys = [tail[1], head[1]]
            zs = [tail[2], head[2]]
            (artist,) = axes.plot(
                xs, ys, zs,
                color=colors[ray.tag],
                lw=1.8,
                ls="--" if ray.dashed else "-",
            )
            record(artist, xs, ys, zs)
            if ray.label:
                ha, va = _LABEL_ANCHORS.get(ray.tag, ("left", "bottom"))
                # A space on the side facing the point, so neighbouring
                # labels anchored on the same tip do not run together.
                padded = f"{ray.label} " if ha == "right" else f" {ray.label}"
                reach = _LABEL_REACH.get(ray.tag, 1.0)
                at = [t + reach * (h - t) for t, h in zip(tail, head)]
                text = axes.text(
                    *at, padded, fontsize=7, color=colors[ray.tag],
                    ha=ha, va=va,
                )
                record(text, [at[0]], [at[1]], [at[2]])

        for point in scene.points:
            px, py, pz = _display((point.x, point.y, point.z))
            (artist,) = axes.plot(
                [px], [py], [pz],
                _GLYPHS[point.glyph],
                color=colors[point.tag],
                ms=4,
                mfc="none" if point.glyph == "out-of-page" else colors[point.tag],
            )
            record(artist, [px], [py], [pz])

        self._apply_limits(scene)
        axes.set_title(scene.title, fontsize=9)
        axes.set_xlabel("d̂  (dispersion)", fontsize=8)
        axes.set_ylabel("ĝ  (grooves)", fontsize=8)
        axes.set_zlabel("n̂  (normal)", fontsize=8)
        axes.tick_params(labelsize=6, pad=0)
        # Five ticks a side rather than nine: the scale is in arbitrary
        # units (only directions carry meaning), so the ticks are there to
        # orient, and nine of them crowded each other.
        from matplotlib.ticker import MaxNLocator

        for axis in (axes.xaxis, axes.yaxis, axes.zaxis):
            axis.set_major_locator(MaxNLocator(5))
        axes.view_init(elev=elev, azim=azim)

        self._canvas.draw_idle()
