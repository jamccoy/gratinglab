"""The 3D conical-diffraction view.

Thin, like every widget here. :mod:`gratinglab.gui.diagram3d` decides every
direction, every surface and every caption; this positions them. It computes no
angle and picks no colour.

**No navigation toolbar on this axes, deliberately.** `NavigationToolbar2QT`'s
zoom rescales an `Axes3D`'s data limits, which breaks the cube bounding box
`set_box_aspect((1, 1, 1))` depends on -- and then the angles this view exists
to show would silently stop being true. A control that can make the picture lie
must not be offered, the same rule that removed the obliquity factor. The 2D
figures keep their toolbars, where pan and zoom cannot distort an angle that
isn't being claimed.

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

__all__ = ["Scene3DPanel"]

#: Named glyph -> matplotlib marker. The scene names what a point *is*; this is
#: the only place that decides what it looks like.
_GLYPHS = {"out-of-page": "o", "strike": "o"}


def _display(v) -> tuple:
    r"""Physical :math:`(\hat{d}, \hat{n}, -\hat{g})` to the frame matplotlib
    draws.

    `Axes3D` always puts *its* z axis vertical, and the axis a reader expects
    to point up is the grating normal -- a grating lies flat. The scene is in
    conventions.md §3's frame, where :math:`\hat{n}` is the *second*
    component, so the widget permutes to
    :math:`(\hat{d}, \hat{g}, \hat{n})`.

    That ordering and not :math:`(\hat{d}, -\hat{g}, \hat{n})` because
    :math:`\hat{d} \times \hat{g} = \hat{n}` -- the alternative is left-handed,
    and a left-handed set would mirror the whole picture.

    A view choice, not a physical one: no length or angle changes under a
    permutation of orthonormal axes, which is why it lives here rather than in
    the scene.
    """
    x, y, z = v
    return x, -z, y


def _view_angles(physical) -> tuple[float, float]:
    """A physical look-direction to matplotlib's ``(elev, azim)``, degrees.

    The inverse of :func:`gratinglab.gui.diagram3d.view_direction`, after the
    display permutation. View math, not physics, so it lives here -- and the
    round trip through `view_direction` is directly assertable, which is what
    a test uses to prove a preset points where its name says.
    """
    dx, dy, dz = _display(np.asarray(physical, dtype=float))
    return float(np.degrees(np.arcsin(-dz))), float(np.degrees(np.arctan2(-dy, -dx)))


class Scene3DPanel(QWidget):
    """One `Axes3D`, and what it takes to keep it honest."""

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
        """Point the camera along a preset's *physical* direction."""
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
        no toolbar on this axes, so it is always true -- it documents the
        coupling and survives someone adding one.
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
        """The one toolbar function worth having, offered on its own."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save 3D view", "geometry.png", "PNG (*.png);;PDF (*.pdf)"
        )
        if not path:
            return
        self._figure.savefig(path, dpi=200)
        QMessageBox.information(self, "Saved", f"Wrote {path}")

    # -- the contract GeometryTab uses --------------------------------

    def show_scene(self, scene) -> None:
        """Draw a scene. The panel keeps it so a caller can read back what was
        drawn without re-deriving it."""
        self._scene = scene
        self._sync_focus_button(scene)
        self._paint(scene)

    def _sync_focus_button(self, scene) -> None:
        """Enabled when it can act; disabled *with its reason* when it cannot.

        At γ = 90° the rim already is the whole scene, so focusing on it would
        zoom out. A disabled control that explains itself is honest; an
        enabled one that silently does nothing is not.
        """
        self._focus_button.setEnabled(scene.rim_available)
        self._focus_button.setToolTip(scene.rim_reason)
        if not scene.rim_available and self._focus_button.isChecked():
            self._focus_button.setChecked(False)

    def _apply_limits(self, scene) -> None:
        """Frame either the whole scene or the cone rim.

        Both are cubes, and both get ``set_box_aspect((1, 1, 1))``, so moving
        between them is a uniform scale plus a translation -- a similarity
        transform. Every angle in the picture survives it exactly; only which
        part of the geometry fills the frame changes. That is why the
        magnification can simply be stated rather than disclaimed.
        """
        use_rim = self._focus == "rim" and scene.rim_available
        (x0, x1), (y0, y1), (z0, z1) = scene.rim_limits if use_rim else scene.limits

        # Limits permute with the coordinates, exactly as the points do.
        self._axes.set_xlim(x0, x1)
        self._axes.set_ylim(-z1, -z0)
        self._axes.set_zlim(y0, y1)
        # A cube, so no axis is stretched relative to another -- without this
        # every angle in the picture is a lie.
        self._axes.set_box_aspect((1.0, 1.0, 1.0))

        self._magnification.setText(
            f"cone rim — {scene.rim_magnification:.0f}×" if use_rim else ""
        )
        # Display ranges, not the physical ones: `_bounds` holds permuted
        # coordinates, and comparing the two frames would hide everything.
        self._hide_artists_outside_the_frame((x0, x1), (-z1, -z0), (y0, y1))

    def _hide_artists_outside_the_frame(self, xr, yr, zr) -> None:
        """Hide primitives with no point inside the current box.

        `Axes3D` does not clip its artists to the axes limits, so at 27x the
        grating patch -- which lies entirely outside the rim frame -- smears
        across the edges instead of disappearing. This is a purely geometric
        test on coordinates already computed; it hides nothing that would have
        been visible, and decides no physics.
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
        colors = diagram3d.TAG_COLORS
        axes = self._axes
        # view_init survives cla(); preserve it explicitly anyway so a redraw
        # never yanks the camera out from under a user mid-inspection.
        elev, azim = axes.elev, axes.azim
        axes.clear()

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
                text = axes.text(*head, ray.label, fontsize=7, color=colors[ray.tag])
                record(text, [head[0]], [head[1]], [head[2]])

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
        axes.tick_params(labelsize=6)
        axes.view_init(elev=elev, azim=azim)

        self._canvas.draw_idle()
