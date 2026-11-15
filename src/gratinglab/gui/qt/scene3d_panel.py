"""The 3D conical-diffraction view, painted with QPainter.

:mod:`gratinglab.gui.diagram3d` decides every direction, surface and caption;
:mod:`gratinglab.gui.camera3d` decides where they land on screen; this module
paints them in the active theme's colours. It computes no angle.

It replaced a matplotlib `Axes3D` on 2026-11-14. The scene was fine; the
chart was the problem: a fixed cube with panes and ticks left the scene a
tenth of the panel, `Axes3D` sorts whole artists rather than faces so rays
and grooves crossed in the wrong order, its text could not move out of a
neighbour's way, and a drag cost ~50 ms a frame. This is the arrangement the
gratinglab app uses (its own orthographic camera over a 2D canvas), ported.

The view is orthographic and every change of framing is a similarity
transform, so zoom, pan, roll and the rim slider all keep every angle true.

Gestures on the canvas, all reaching the one `Camera`:

- drag: orbit; Shift+drag or right/middle drag: pan; Alt+drag: roll;
- wheel: zoom about the cursor; double-click or 0: recentre;
- arrows orbit 5° (Shift 1°), + and - zoom, ] and [ roll.

Called by :class:`gratinglab.gui.qt.geometry_tab.GeometryTab` through
`show_scene`, and it reads nothing else.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from .. import camera3d, diagram3d, label_layout
from ..camera3d import Camera
from . import recent, theming

__all__ = ["Scene3DPanel", "SceneCanvas", "LIGHT"]

#: Where the light comes from, in the world frame: above the grating, toward
#: the viewer's default side. A drawing choice. It only has to separate the
#: two faces of a groove the way a lamp would; nothing physical reads it.
LIGHT = np.array([0.45, 1.0, 0.55]) / np.linalg.norm([0.45, 1.0, 0.55])
#: Shading = AMBIENT + (1 - AMBIENT) * max(0, n.L): a face turned away from the
#: light still reads as a surface rather than a hole.
AMBIENT = 0.35
#: Each ray and curve is cut into this many pieces so that it sorts against
#: the grooves piece by piece. Rays leave the grating surface, so a whole-ray
#: sort would put each one entirely in front of or behind the patch.
SEGMENTS = 12
#: Line widths in pixels, by what the line is.
WIDTHS = {"order": 2.0, "zero": 2.0, "blaze": 2.2, "incident": 1.8,
          "normal": 1.4, "axis": 1.1, "cone": 1.1}
#: Cone generators drawn over the mantle for the sense of a ruled surface.
#: Many more read as hatching at rim zoom.
GENERATORS = 12
#: Orbit per pixel of drag, and per arrow press (fine with Shift).
DRAG_DEG, KEY_DEG, FINE_DEG = 0.4, 5.0, 1.0
#: Zoom per wheel notch and per +/- press (fine with Shift).
WHEEL_STEP, KEY_ZOOM, FINE_ZOOM = 1.15, 1.25, 1.05
#: Pixels reserved round the fit for labels, which do not scale.
LABEL_PAD = 34.0


def _mix(a: QColor, b: QColor, t: float) -> QColor:
    t = float(np.clip(t, 0.0, 1.0))
    return QColor.fromRgbF(
        a.redF() + (b.redF() - a.redF()) * t,
        a.greenF() + (b.greenF() - a.greenF()) * t,
        a.blueF() + (b.blueF() - a.blueF()) * t,
    )


class SceneCanvas(QWidget):
    """The painted view. Owns the camera; the panel owns the controls."""

    #: Emitted after any change of framing, so the panel's readout follows.
    framed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMinimumSize(200, 160)
        self.setToolTip(
            "Drag to orbit, Shift+drag to pan, Alt+drag to roll, wheel to zoom, "
            "double-click to recentre."
        )
        self.scene = None
        self.camera = Camera()
        self.rim_offered = True
        self.lod = "fine"
        self.frame: camera3d.Frame | None = None
        #: Label rectangles placed on the last paint, for the tests.
        self.label_boxes: list[tuple[str, QRectF]] = []
        self._drag: tuple[QPointF, Qt.MouseButton, Qt.KeyboardModifier] | None = None

    # -- framing ------------------------------------------------------

    def set_camera(self, camera: Camera) -> None:
        self.camera = camera
        self.update()
        self.framed.emit()

    def _fit_points(self, scene) -> np.ndarray:
        """What the frame is fitted to: the patch, the rim, and every ray
        except the incident one's far end, which is as long as the rest of the
        scene again and would halve it. Read from the scene, so the fit
        follows whatever is drawn."""
        pts = []
        for s in scene.surfaces:
            if s.lod == self.lod and s.tag == "surface":
                pts.append(np.column_stack([s.x.ravel(), s.y.ravel(), s.z.ravel()]))
        for c in scene.curves:
            if c.lod == self.lod and c.tag == "cone":
                pts.append(np.column_stack([c.x, c.y, c.z]))
        for r in scene.rays:
            pts.append(np.asarray([r.head]))
            if r.tag != "incident":
                pts.append(np.asarray([r.origin]))
        return np.vstack(pts)

    def _frame_for(self, width: float, height: float) -> camera3d.Frame:
        return camera3d.fit(
            self.camera, self._fit_points(self.scene), self.scene.gamma,
            (0.0, 0.0, width, height), pad=LABEL_PAD, rim_offered=self.rim_offered,
        )

    # -- painting -----------------------------------------------------

    def paintEvent(self, _event) -> None:
        painter = QPainter(self)
        try:
            self.render_into(painter, self.width(), self.height())
        finally:
            painter.end()

    def render_into(self, painter: QPainter, width: float, height: float) -> None:
        """Paint one frame; `save_image` calls this at a larger size."""
        theme = theming.current()
        colors = {**diagram3d.TAG_COLORS, **theme.diagram}
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(QRectF(0, 0, width, height), QColor(theme.axes))
        self.label_boxes = []
        if self.scene is None:
            return
        frame = self.frame = self._frame_for(width, height)
        scene = self.scene

        # --- the cone mantle, behind everything: it is a faint fill and
        # sorting it against the rays would only cost time.
        mantle = [s for s in scene.surfaces if s.tag == "cone" and s.lod == self.lod]
        for s in mantle:
            rim = np.column_stack([s.x[:, 1], s.y[:, 1], s.z[:, 1]])
            apex = np.array([[s.x[0, 0], s.y[0, 0], s.z[0, 0]]])
            xy, _ = camera3d.project(frame, np.vstack([apex, rim]))
            color = QColor(colors["cone"])
            color.setAlphaF(s.alpha)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(color)
            painter.drawPolygon(QPolygonF([QPointF(*p) for p in xy]))
            step = max(1, (len(rim) - 1) // GENERATORS)
            gen = QColor(colors["cone"])
            gen.setAlphaF(0.22)
            painter.setPen(QPen(gen, 1.0))
            for p in xy[1::step][:GENERATORS]:
                painter.drawLine(QPointF(*xy[0]), QPointF(*p))

        # --- everything else, sorted far to near, piece by piece.
        items: list[tuple[float, int, object]] = []
        self._collect_patch(frame, colors, items)
        self._collect_lines(frame, colors, items)
        items.sort(key=lambda it: (-it[0], it[1]))
        for _depth, _rank, draw in items:
            draw(painter)

        self._paint_points(painter, frame, colors)
        self._paint_labels(painter, frame, colors, theme, width, height)
        self._paint_title(painter, theme, width)

    def _collect_patch(self, frame, colors, items) -> None:
        """The grating: one quad per sample interval of the ruled surface,
        shaded by how squarely it faces the light, so a blaze facet and its
        anti-blaze face come out in two tones and a sinusoid in a gradient."""
        lit, shade = QColor(colors["facet"]), QColor(colors["surface"])
        # On the way to the rim the patch fades and then goes: at the rim it
        # lies far outside the frame, and drawn at that magnification it is a
        # sawtooth wall behind the fan. The same schedule as the app's.
        rim = self.camera.rim if self.rim_offered else 0.0
        fade = 1.0 - float(np.clip((rim - 0.25) / 0.30, 0.0, 1.0))
        if fade <= 0.0:
            return
        for s in self.scene.surfaces:
            if s.tag != "surface" or s.lod != self.lod:
                continue
            a = np.column_stack([s.x[:, 0], s.y[:, 0], s.z[:, 0]])
            b = np.column_stack([s.x[:, 1], s.y[:, 1], s.z[:, 1]])
            xa, da = camera3d.project(frame, a)
            xb, db = camera3d.project(frame, b)
            along = b[0] - a[0]
            for i in range(len(a) - 1):
                across = a[i + 1] - a[i]
                normal = np.cross(along, across)
                length = np.linalg.norm(normal)
                if length == 0.0:
                    continue
                normal /= length
                if normal[1] < 0:
                    normal = -normal  # the face a viewer above the grating sees
                light = AMBIENT + (1 - AMBIENT) * max(0.0, float(normal @ LIGHT))
                # Hue from the theme's two surface colours, brightness from
                # the light, so the two faces of a groove separate in every
                # theme, light or dark.
                color = _mix(shade, lit, (light - AMBIENT) / (1 - AMBIENT))
                h, s_, v, _a = color.getHsvF()
                color = QColor.fromHsvF(h, s_, float(np.clip(v * light, 0.0, 1.0)))
                color.setAlphaF(max(s.alpha, 0.9) * fade)
                quad = QPolygonF([QPointF(*xa[i]), QPointF(*xa[i + 1]),
                                  QPointF(*xb[i + 1]), QPointF(*xb[i])])
                depth = float((da[i] + da[i + 1] + db[i] + db[i + 1]) / 4)

                def draw(p, quad=quad, color=color):
                    # A cosmetic pen in the fill colour closes the hairline
                    # seams antialiasing leaves between neighbouring quads.
                    p.setPen(QPen(color, 1.0))
                    p.setBrush(color)
                    p.drawPolygon(quad)

                items.append((depth, 0, draw))

    def _collect_lines(self, frame, colors, items) -> None:
        segments = []
        for c in self.scene.curves:
            if c.lod != self.lod and c.tag == "cone":
                continue
            pts = np.column_stack([c.x, c.y, c.z])
            segments.append((pts, c.tag, c.dashed, 0.75 if c.tag == "cone" else 0.45))
        for r in self.scene.rays:
            t = np.linspace(0.0, 1.0, SEGMENTS + 1)[:, None]
            pts = np.asarray(r.origin) + t * (np.asarray(r.head) - np.asarray(r.origin))
            segments.append((pts, r.tag, r.dashed, 1.0))
        for pts, tag, dashed, alpha in segments:
            xy, depth = camera3d.project(frame, pts)
            color = QColor(colors[tag])
            color.setAlphaF(alpha)
            pen = QPen(color, WIDTHS.get(tag, 1.2))
            pen.setCapStyle(Qt.PenCapStyle.FlatCap)
            if dashed:
                pen.setStyle(Qt.PenStyle.DashLine)
            for i in range(len(xy) - 1):
                line = (QPointF(*xy[i]), QPointF(*xy[i + 1]))

                def draw(p, line=line, pen=pen):
                    p.setPen(pen)
                    p.drawLine(*line)

                items.append((float((depth[i] + depth[i + 1]) / 2), 1, draw))

    def _paint_points(self, painter, frame, colors) -> None:
        for point in self.scene.points:
            (xy,), _ = camera3d.project(frame, [(point.x, point.y, point.z)])
            color = QColor(colors[point.tag])
            painter.setPen(QPen(color, 1.4))
            painter.setBrush(Qt.BrushStyle.NoBrush if point.glyph == "out-of-page" else color)
            painter.drawEllipse(QPointF(*xy), 3.5, 3.5)

    # -- labels ---------------------------------------------------------

    def _labels(self):
        """``(world point, text, tag)`` for everything named.

        The scene names the axes, the incident ray, the blaze direction and
        m = 0. It left the other orders unnamed because matplotlib's 3D text
        could not move out of a neighbour's way; this view can, so each order
        is named here from the ray's own ``order``, the scene's number."""
        out = []
        for r in self.scene.rays:
            text = r.label
            if text is None and r.order is not None:
                text = f"{r.order:+d}"
            if not text:
                continue
            tail, head = np.asarray(r.origin), np.asarray(r.head)
            if r.tag == "incident":
                # k_i ends AT the strike point, where everything starts; its
                # name belongs at the far end, where it comes from.
                tail, head = head, tail
            out.append((tail, head, text, r.tag))
        return out

    def _paint_labels(self, painter, frame, colors, theme, width, height) -> None:
        """Name everything, covering no other label and no other tip.

        Labels whose tips crowd together go in a column beside the crowd,
        laid out by :mod:`gratinglab.gui.label_layout`: at grazing incidence
        every order's tip lands in one knot, and labels set round each tip
        filled it. The rest go at the first of several spots round their tip
        that is free, orders last, since there are many of them."""
        font = QFont(painter.font())
        font.setPointSizeF(8.5)
        painter.setFont(font)
        metrics = QFontMetricsF(font)
        bounds = QRectF(2, 2, width - 4, height - 4)

        labels = []
        for tail, world, text, tag in self._labels():
            (t_xy, xy), _ = camera3d.project(frame, [tail, world])
            size = (metrics.horizontalAdvance(text) + 4, metrics.height())
            labels.append((t_xy, (float(xy[0]), float(xy[1])), text, tag, size))
        tips = [lab[1] for lab in labels]
        markers = []
        for point in self.scene.points:
            (xy,), _ = camera3d.project(frame, [(point.x, point.y, point.z)])
            markers.append((float(xy[0]), float(xy[1])))
        keep_clear = tips + markers
        # Every labelled line as drawn on screen, tail to tip, for a column
        # to keep out of.
        lines = [((float(lab[0][0]), float(lab[0][1])), lab[1]) for lab in labels]

        # --- crowds first, each as one column, while there is room for it.
        placed: list[QRectF] = []
        done: set[int] = set()
        frame_rect = bounds.getRect()
        on_screen = [i for i, tip in enumerate(tips) if bounds.contains(QPointF(*tip))]
        for crowd in label_layout.crowds([tips[i] for i in on_screen]):
            members = [on_screen[k] for k in crowd]
            rects = label_layout.column(
                members, tips, [lab[4] for lab in labels], frame_rect,
                keep_clear, [b.getRect() for b in placed], lines,
            )
            if rects is None:
                continue  # no room for a column; each takes its own chance below
            for i, rect in rects.items():
                box = QRectF(*rect)
                self._draw_label(painter, colors, labels[i], box, leader=True)
                placed.append(box)
                done.add(i)

        # --- everything else, round its own tip.
        rest = sorted((i for i in range(len(labels)) if i not in done),
                      key=lambda i: labels[i][3] in ("order", "zero"))
        for i in rest:
            t_xy, (x, y), text, tag, (w, h) = labels[i]
            others = [p for k, p in enumerate(keep_clear) if k != i]
            # First choice: just past the tip, continuing the line, where the
            # label reads as the line's own and covers none of it.
            dx, dy = x - t_xy[0], y - t_xy[1]
            norm = float(np.hypot(dx, dy))
            ux, uy = (dx / norm, dy / norm) if norm > 1e-6 else (1.0, 0.0)
            reach = 6 + 0.5 * abs(ux) * w + 0.5 * abs(uy) * h
            beyond = (ux * reach - w / 2, uy * reach - h / 2)
            spots = [beyond, (6, -h - 2), (6, 2), (-w - 6, -h - 2), (-w - 6, 2),
                     (-w / 2, -h - 8), (-w / 2, 8), (14, -h / 2), (-w - 14, -h / 2)]
            # Then rings further out, joined to the tip by a leader line. A
            # label never covers another: one that cannot be placed at all is
            # left out, which the tests would notice for an order.
            for ring in (1, 2, 3, 4, 6):
                r = ring * (h + 6) + 10
                for ang in np.radians(np.arange(0, 360, 30)):
                    spots.append((np.cos(ang) * r - w / 2, -np.sin(ang) * r - h / 2))
            box, leader = None, False
            for n, (dx, dy) in enumerate(spots):
                candidate = QRectF(x + dx, y + dy, w, h)
                if bounds.contains(candidate) \
                        and not any(candidate.intersects(b) for b in placed) \
                        and not any(label_layout.covers(candidate.getRect(), p) for p in others):
                    box, leader = candidate, n >= 9
                    break
            if box is None or not bounds.contains(QPointF(x, y)):
                continue  # nowhere free, or its tip is off screen
            self._draw_label(painter, colors, labels[i], box, leader)
            placed.append(box)

    def _draw_label(self, painter, colors, label, box, leader) -> None:
        """One label in its box, joined to its tip when it is set apart."""
        _tail, (x, y), text, tag, _size = label
        if leader:
            lead = QColor(colors[tag])
            lead.setAlphaF(0.5)
            painter.setPen(QPen(lead, 0.8))
            # To the nearest point of the box, which keeps clear of the tip.
            painter.drawLine(QPointF(x, y), QPointF(
                min(max(x, box.left()), box.right()), min(max(y, box.top()), box.bottom())))
        self.label_boxes.append((text, box))
        painter.setPen(QColor(colors[tag]))
        painter.drawText(box, Qt.AlignmentFlag.AlignCenter, text)

    def _paint_title(self, painter, theme, width) -> None:
        font = QFont(painter.font())
        font.setPointSizeF(9.5)
        painter.setFont(font)
        painter.setPen(QColor(theme.muted))
        painter.drawText(QRectF(10, 6, width - 20, 20), Qt.AlignmentFlag.AlignLeft, self.scene.title)

    # -- gestures -----------------------------------------------------

    def mousePressEvent(self, event) -> None:
        self._drag = (event.position(), event.button(), event.modifiers())
        self.lod = "coarse"
        self.setFocus()

    def mouseMoveEvent(self, event) -> None:
        if self._drag is None:
            return
        last, button, mods = self._drag
        pos = event.position()
        dx, dy = pos.x() - last.x(), pos.y() - last.y()
        self._drag = (pos, button, mods)
        cam = self.camera
        if mods & Qt.KeyboardModifier.AltModifier:
            cam = cam.rolled(dx * DRAG_DEG)
        elif button != Qt.MouseButton.LeftButton or mods & Qt.KeyboardModifier.ShiftModifier:
            cam = cam.panned(dx, dy)
        else:
            cam = cam.orbited(-dx * DRAG_DEG, dy * DRAG_DEG)
        self.set_camera(cam)

    def mouseReleaseEvent(self, _event) -> None:
        self._drag = None
        self.lod = "fine"  # unconditional: fine is always the safe direction
        self.update()

    def mouseDoubleClickEvent(self, _event) -> None:
        self.set_camera(self.camera.recentred())

    def wheelEvent(self, event) -> None:
        if self.frame is None:
            return
        notches = event.angleDelta().y() / 120.0
        at = event.position()
        self.set_camera(camera3d.zoom_about(
            self.camera, self.frame, WHEEL_STEP ** notches, (at.x(), at.y())))

    def keyPressEvent(self, event) -> None:
        fine = bool(event.modifiers() & Qt.KeyboardModifier.ShiftModifier)
        deg = FINE_DEG if fine else KEY_DEG
        zoom = FINE_ZOOM if fine else KEY_ZOOM
        centre = (self.width() / 2, self.height() / 2)
        key, cam = event.key(), self.camera
        moves = {
            Qt.Key.Key_Left: lambda: cam.orbited(deg, 0),
            Qt.Key.Key_Right: lambda: cam.orbited(-deg, 0),
            Qt.Key.Key_Up: lambda: cam.orbited(0, deg),
            Qt.Key.Key_Down: lambda: cam.orbited(0, -deg),
            Qt.Key.Key_BracketRight: lambda: cam.rolled(deg),
            Qt.Key.Key_BracketLeft: lambda: cam.rolled(-deg),
            Qt.Key.Key_0: cam.recentred,
            Qt.Key.Key_Home: cam.recentred,
        }
        if key in (Qt.Key.Key_Plus, Qt.Key.Key_Equal, Qt.Key.Key_Minus) and self.frame:
            factor = zoom if key != Qt.Key.Key_Minus else 1 / zoom
            self.set_camera(camera3d.zoom_about(cam, self.frame, factor, centre))
        elif key in moves:
            self.set_camera(moves[key]())
        else:
            super().keyPressEvent(event)


class Scene3DPanel(QWidget):
    """The painted view with its controls: view, rim focus, recentre, save."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._scene = None
        self._canvas = SceneCanvas(self)
        self._canvas.framed.connect(self._sync_readout)
        theming.subscribe(self._on_theme)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._canvas, 1)
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

        row.addWidget(QLabel("Zoom to rim:"))
        self._rim_slider = QSlider(Qt.Orientation.Horizontal)
        self._rim_slider.setRange(0, 1000)
        self._rim_slider.setMinimumWidth(120)
        self._rim_slider.valueChanged.connect(self._on_rim)
        row.addWidget(self._rim_slider, 1)

        self._magnification = QLabel("")
        self._magnification.setMinimumWidth(48)
        row.addWidget(self._magnification)

        recentre = QPushButton("Recentre")
        recentre.setToolTip("Undo pan, zoom and roll; keep the view and the rim focus. (0)")
        recentre.clicked.connect(lambda: self._canvas.set_camera(self._canvas.camera.recentred()))
        row.addWidget(recentre)

        save = QPushButton("Save image…")
        save.clicked.connect(self.save_image)
        row.addWidget(save)
        return strip

    # -- reactions ----------------------------------------------------

    def _apply_preset(self, name: str) -> None:
        """Look along a preset's physical direction, with the free transform
        cleared. The rim preset exists for the rim, so it also zooms there."""
        if name not in diagram3d.PRESET_VIEWS:
            return
        yaw, pitch = camera3d.orbit_for(diagram3d.PRESET_VIEWS[name])
        camera = replace(self._canvas.camera.recentred(), yaw=yaw, pitch=pitch)
        if name == diagram3d.RIM_PRESET and self._rim_slider.isEnabled():
            self._rim_slider.blockSignals(True)
            self._rim_slider.setValue(1000)
            self._rim_slider.blockSignals(False)
            camera = replace(camera, rim=1.0)
        self._canvas.set_camera(camera)

    def _on_rim(self, value: int) -> None:
        self._canvas.set_camera(replace(self._canvas.camera, rim=value / 1000.0))

    def _sync_readout(self) -> None:
        if self._scene is None:
            return
        frame = self._canvas._frame_for(self._canvas.width(), self._canvas.height())
        self._magnification.setText(f"{frame.magnification:.1f}×")

    def _on_theme(self, _theme) -> None:
        self._canvas.update()

    def show_scene(self, scene) -> None:
        """Draw a scene. The panel keeps it, so a caller can read back what
        was drawn without re-deriving it."""
        self._scene = scene
        self._canvas.scene = scene
        # Disabled with its reason rather than enabled and doing nothing: at
        # gamma = 90 deg the rim is the whole scene.
        self._canvas.rim_offered = scene.rim_available
        self._rim_slider.setEnabled(scene.rim_available)
        self._rim_slider.setToolTip(scene.rim_reason or "Zoom from the whole scene to the cone rim.")
        if not scene.rim_available and self._rim_slider.value():
            self._rim_slider.setValue(0)
        self._canvas.update()
        self._sync_readout()

    def save_image(self) -> None:
        """Save the view as drawn, at twice the size on screen."""
        path, _ = QFileDialog.getSaveFileName(
            self, "Save 3D view", recent.start_in("exports", "geometry.png"), "PNG (*.png)",
        )
        if not path:
            return
        recent.remember("exports", path)
        self.render_image(2.0).save(path)
        QMessageBox.information(self, "Saved", f"Wrote {path}")

    def render_image(self, scale: float = 1.0):
        """The view as a QImage, at ``scale`` times the canvas's size."""
        from PySide6.QtGui import QImage

        w, h = self._canvas.width(), self._canvas.height()
        image = QImage(int(w * scale), int(h * scale), QImage.Format.Format_ARGB32)
        painter = QPainter(image)
        try:
            painter.scale(scale, scale)
            self._canvas.render_into(painter, w, h)
        finally:
            painter.end()
        return image
