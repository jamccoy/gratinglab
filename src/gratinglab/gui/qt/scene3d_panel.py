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
from PySide6.QtWidgets import QLabel, QVBoxLayout, QWidget

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


class Scene3DPanel(QWidget):
    """One `Axes3D`, and what it takes to keep it honest."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._scene = None
        self._build()

    def _build(self) -> None:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(6, 5), layout="constrained")
        self._axes = self._figure.add_subplot(projection="3d")
        self._canvas = FigureCanvasQTAgg(self._figure)

        self._hint = QLabel("Drag to rotate.")
        self._hint.setStyleSheet("color: #666")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._canvas)
        layout.addWidget(self._hint)

        elev, azim = diagram3d.PRESET_VIEWS["oblique"]
        self._axes.view_init(elev=elev, azim=azim)

    # -- the contract GeometryTab uses --------------------------------

    def show_scene(self, scene) -> None:
        """Draw a scene. The panel keeps it so a caller can read back what was
        drawn without re-deriving it."""
        self._scene = scene
        self._paint(scene)

    # -- drawing ------------------------------------------------------

    def _paint(self, scene) -> None:
        colors = diagram3d.TAG_COLORS
        axes = self._axes
        # view_init survives cla(); preserve it explicitly anyway so a redraw
        # never yanks the camera out from under a user mid-inspection.
        elev, azim = axes.elev, axes.azim
        axes.clear()

        for surface in scene.surfaces:
            if surface.lod != "fine":
                continue
            sx, sy, sz = _display((surface.x, surface.y, surface.z))
            axes.plot_surface(
                sx, sy, sz,
                color=colors[surface.tag],
                alpha=surface.alpha,
                linewidth=0,
                antialiased=False,
                shade=False,
            )

        for curve in scene.curves:
            if curve.lod != "fine":
                continue
            cx, cy, cz = _display((curve.x, curve.y, curve.z))
            axes.plot(
                cx, cy, cz,
                color=colors[curve.tag],
                lw=1.0,
                ls="--" if curve.dashed else "-",
                alpha=0.7,
            )

        for ray in scene.rays:
            tail = _display(ray.origin)
            head = _display(ray.head)
            axes.plot(
                [tail[0], head[0]], [tail[1], head[1]], [tail[2], head[2]],
                color=colors[ray.tag],
                lw=1.8,
                ls="--" if ray.dashed else "-",
            )
            if ray.label:
                axes.text(*head, ray.label, fontsize=7, color=colors[ray.tag])

        for point in scene.points:
            px, py, pz = _display((point.x, point.y, point.z))
            axes.plot(
                [px], [py], [pz],
                _GLYPHS[point.glyph],
                color=colors[point.tag],
                ms=4,
                mfc="none" if point.glyph == "out-of-page" else colors[point.tag],
            )

        # Limits permute with the coordinates. The box is a cube, so the
        # permutation only relabels which range goes on which axis.
        (x0, x1), (y0, y1), (z0, z1) = scene.limits
        axes.set_xlim(x0, x1)
        axes.set_ylim(-z1, -z0)
        axes.set_zlim(y0, y1)
        # A cube, so no axis is stretched relative to another -- without this
        # every angle in the picture is a lie.
        axes.set_box_aspect((1.0, 1.0, 1.0))

        axes.set_title(scene.title, fontsize=9)
        axes.set_xlabel("d̂  (dispersion)", fontsize=8)
        axes.set_ylabel("ĝ  (grooves)", fontsize=8)
        axes.set_zlabel("n̂  (normal)", fontsize=8)
        axes.tick_params(labelsize=6)
        axes.view_init(elev=elev, azim=azim)

        self._canvas.draw_idle()
