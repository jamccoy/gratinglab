"""The Mirror tab: the coating stack as a flat mirror, before any grooves.

Reads the stack from the geometry dock's Surface group and its own fields for
the scan, and shows the reflectivity, how deep the field reaches into each
material, and whether the outermost layer is thick enough to treat as bulk.
The answer to that is one line above the plots; when it is yes, a button
makes the layer the grating material and drops the layers beneath it, which
is the simplification the tab exists to justify. Everything it computes and
says is in :mod:`gratinglab.gui.mirror_view`; this is widgets.

There is no Calculate button. The tab recomputes a moment after any of its
fields or the dock changes, and only while it is showing: a scan takes about a
second in the GUI thread, which is too long to spend on every keystroke made
in another tab. It is not a solver and has no worker.

Built by :class:`~gratinglab.gui.qt.main_window.MainWindow`, which hands it
the dock.
"""

from __future__ import annotations

from endstation.qt.canvas import restyle_toolbar
from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .. import provenance
from ..mirror_view import (
    MirrorResult,
    MirrorSettings,
    compute,
    draw,
    is_bulk_throughout,
    lines,
    verdict,
)
from ..state import FormErrors
from ..theme import provenance_colors, style_figure
from . import theming
from .sizing import fit_width_to_contents

__all__ = ["MirrorTab"]

#: How long after the last edit the tab recomputes, ms.
_SETTLE_MS = 300

_ANGLE_FIELDS = ("graze_start", "graze_stop", "graze_count", "wavelength")
_OWN_WAVELENGTHS = ("wavelength_start", "wavelength_stop", "wavelength_count")


class MirrorTab(QWidget):
    """Reflectivity and penetration of the dock's coating stack."""

    def __init__(self, dock, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        #: The geometry dock: read for the stack and the mount, written by
        #: the hand-off button.
        self._dock = dock
        self._result: MirrorResult | None = None
        self._lines: tuple = ()
        self._fields: dict[str, QLineEdit] = {}
        self._stale = True
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(_SETTLE_MS)
        self._timer.timeout.connect(self._settled)
        self._build()
        dock.changed.connect(self.schedule)

    # -- construction ----------------------------------------------------

    def _build(self) -> None:
        controls = QWidget()
        column = QVBoxLayout(controls)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(self._build_scan_group())
        self._hand_off = QPushButton("Model the top layer as bulk")
        self._hand_off.setToolTip(
            "Make the outermost layer the dock's Grating material and remove\n"
            "the layers. Offered when that layer reads as bulk across the scan,\n"
            "so what lies beneath it does not change the reflectivity."
        )
        self._hand_off.setEnabled(False)
        self._hand_off.clicked.connect(self._model_as_bulk)
        column.addWidget(self._hand_off)
        column.addStretch(1)

        scroller = QScrollArea()
        scroller.setWidget(controls)
        scroller.setWidgetResizable(True)
        fit_width_to_contents(scroller)

        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(self._build_plot())
        right.addWidget(self._build_panel())
        right.setStretchFactor(0, 1)
        right.setSizes([560, 160])

        body = QSplitter(Qt.Orientation.Horizontal)
        body.addWidget(scroller)
        body.addWidget(right)
        body.setStretchFactor(1, 1)
        body.setSizes([260, 700])

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(body)

    def _build_scan_group(self) -> QGroupBox:
        defaults = MirrorSettings()
        group = QGroupBox("Scan")
        form = QFormLayout(group)

        self._scan = QComboBox()
        self._scan.addItem("Graze angle", "angle")
        self._scan.addItem("Wavelength", "wavelength")
        self._scan.currentIndexChanged.connect(self._on_scan_change)
        form.addRow("Scan", self._scan)

        def entry(label: str, key: str, tip: str) -> None:
            field = QLineEdit(getattr(defaults, key))
            field.setToolTip(tip)
            field.textEdited.connect(self.schedule)
            field.returnPressed.connect(self.calculate)
            self._fields[key] = field
            form.addRow(label, field)

        entry("From (deg)", "graze_start", "Graze angle from the surface; 0 is grazing.")
        entry("To (deg)", "graze_stop", "Graze angle from the surface.")
        entry("Points", "graze_count", "Angles in the scan, at most 1000.")
        entry("At wavelength (nm)", "wavelength", "The one wavelength an angle scan uses.")
        entry("At graze (deg)", "graze", "The one graze angle a wavelength scan uses.")

        self._dock_wavelengths = QCheckBox("Use the grating's wavelengths")
        self._dock_wavelengths.setChecked(defaults.dock_wavelengths)
        self._dock_wavelengths.setToolTip("Scan the dock's Wavelengths, or the three fields below.")
        self._dock_wavelengths.toggled.connect(self._on_scan_change)
        form.addRow(self._dock_wavelengths)
        entry("From (nm)", "wavelength_start", "Shortest wavelength of this tab's own scan.")
        entry("To (nm)", "wavelength_stop", "Longest wavelength of this tab's own scan.")
        entry("Points", "wavelength_count", "Wavelengths in the scan, at most 1000.")

        self._polarization = QComboBox()
        for label, value in (("Unpolarized", "unpolarized"), ("s", "s"), ("p", "p")):
            self._polarization.addItem(label, value)
        self._polarization.setToolTip("s and p are local to the mirror's plane of incidence.")
        self._polarization.currentIndexChanged.connect(self.schedule)
        form.addRow("Polarization", self._polarization)

        entry("Bulk within", "tolerance",
              "The outermost layer reads as bulk once the layers beneath it\n"
              "change the reflectivity by less than this, relative.")
        self._log = QCheckBox("Log reflectivity")
        self._log.setChecked(True)
        self._log.toggled.connect(self._redraw)
        form.addRow(self._log)
        self._set_field_states()
        return group

    def _build_plot(self) -> QWidget:
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(8, 5), layout="constrained")
        grid = self._figure.add_gridspec(2, 1, height_ratios=(1, 1))
        self._r_axes = self._figure.add_subplot(grid[0])
        self._depth_axes = self._figure.add_subplot(grid[1], sharex=self._r_axes)
        self._canvas = FigureCanvasQTAgg(self._figure)
        style_figure(self._figure, theming.current())
        theming.subscribe(self._on_theme)

        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        # The answer, above the plots that support it.
        self._verdict = QLabel("")
        self._verdict.setContentsMargins(8, 6, 8, 2)
        self._verdict.setTextFormat(Qt.TextFormat.RichText)
        bold = self._verdict.font()
        bold.setBold(True)
        self._verdict.setFont(bold)
        layout.addWidget(self._verdict)
        self._toolbar = NavigationToolbar2QT(self._canvas, panel)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, 1)  # the plots take the spare height
        return panel

    def _build_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 4, 8, 8)
        heading = QLabel("Mirror")
        bold = heading.font()
        bold.setBold(True)
        heading.setFont(bold)
        self._panel = QTextEdit()
        self._panel.setReadOnly(True)
        layout.addWidget(heading)
        layout.addWidget(self._panel)
        return panel

    # -- what the window and tests ask ---------------------------------------

    def settings(self) -> MirrorSettings:
        values = {key: field.text() for key, field in self._fields.items()}
        return MirrorSettings(
            scan=self._scan.currentData(),
            polarization=self._polarization.currentData(),
            dock_wavelengths=self._dock_wavelengths.isChecked(),
            **values,
        )

    def schedule(self, *_args) -> None:
        """Something changed: recompute once the edits settle, if showing."""
        self._stale = True
        if self.isVisible():
            self._timer.start()

    def calculate(self) -> None:
        """Compute now, from the dock and the tab's fields."""
        self._timer.stop()
        self._stale = False
        try:
            self._result = compute(self._dock.read_form(), self.settings())
        except FormErrors as exc:
            self._result = None
            self._paint(provenance.error_lines(exc.errors))
        except (ValueError, KeyError) as exc:
            # A material that does not cover the scan, or one not installed.
            self._result = None
            self._paint((provenance.Line(f"{exc}\n", "bad"),))
        self._redraw()

    # -- internals -------------------------------------------------------

    def showEvent(self, event) -> None:  # noqa: N802 - Qt's name
        super().showEvent(event)
        if self._stale:
            self._timer.start(0)

    def _settled(self) -> None:
        self._timer.setInterval(_SETTLE_MS)  # showEvent may have shortened it
        if self.isVisible():
            self.calculate()

    def _on_scan_change(self, *_args) -> None:
        self._set_field_states()
        self.schedule()

    def _set_field_states(self) -> None:
        angle = self._scan.currentData() == "angle"
        for key in _ANGLE_FIELDS:
            self._fields[key].setEnabled(angle)
        self._fields["graze"].setEnabled(not angle)
        self._dock_wavelengths.setEnabled(not angle)
        own = not angle and not self._dock_wavelengths.isChecked()
        for key in _OWN_WAVELENGTHS:
            self._fields[key].setEnabled(own)

    def _model_as_bulk(self) -> None:
        if self._result is None or not is_bulk_throughout(self._result):
            return
        mirror = self._result.mirror
        top = mirror.layers[-1].material
        if self._confirm_hand_off(top, len(mirror.layers), mirror.substrate):
            self._dock.replace_stack(top)

    def _confirm_hand_off(self, top: str, layers: int, substrate: str) -> bool:
        """Its own method so a test can answer it. The dock has no undo."""
        removed = "the layer" if layers == 1 else f"all {layers} layers"
        answer = QMessageBox.question(
            self,
            "Model the top layer as bulk",
            f"Make {top} the grating material, replacing {substrate}, and remove "
            f"{removed} from the stack? This cannot be undone.",
        )
        return answer == QMessageBox.StandardButton.Yes

    def _paint(self, lines_) -> None:
        self._lines = tuple(lines_)
        self._panel.setHtml(provenance.to_html(self._lines, provenance_colors(theming.current())))

    def _on_theme(self, theme) -> None:
        restyle_toolbar(self._toolbar, theme)
        style_figure(self._figure, theme)
        self._redraw()
        if self._lines:
            self._paint(self._lines)

    def _redraw(self, *_args) -> None:
        result = self._result
        self._hand_off.setEnabled(result is not None and is_bulk_throughout(result))
        if result is not None:
            draw(self._r_axes, self._depth_axes, result, theming.current(), self._log.isChecked())
            self._paint(lines(result))
            self._verdict.setText(provenance.to_html((verdict(result),), provenance_colors(theming.current())))
        else:
            for axes in (self._r_axes, self._depth_axes):
                axes.clear()
            style_figure(self._figure, theming.current())
            self._verdict.setText("")
        self._canvas.draw_idle()
