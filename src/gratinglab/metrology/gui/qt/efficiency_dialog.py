"""
Diffraction efficiency of the groove currently previewed on the Boundary tab.

The point of the two halves being one package, made visible. Everything here
was reachable before by exporting a .ggp, opening the other window and loading
it back -- except the period, which the file cannot carry and which a user then
had to retype from the metrics panel. This path keeps it, because
`BoundaryProfile.to_problem` carries the number the scan measured.

Deliberately a *preview*, not a second efficiency application. It offers one
mount, one solver and no convergence check, and says as much: a result whose
convergence was never demonstrated reports itself as not defensible, and this
dialog does not demonstrate it. Anything load-bearing belongs in the main
window, where the convergence harness and the comparison view live.

Synchronous, following the Boundary tab's own reasoning rather than the main
window's: a scalar solve over a few hundred wavelengths is well under a second
(3.4 s at the extreme this dialog allows), so a worker thread would buy a
progress bar nobody waits on. The wait cursor is the whole concession.
"""
from __future__ import annotations

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFormLayout, QGroupBox,
    QHBoxLayout, QLabel, QPushButton, QSpinBox, QVBoxLayout, QWidget,
)
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from gratinglab import materials
from gratinglab.illumination import Illumination
from gratinglab.solvers import UnsupportedConfiguration, get_solver

__all__ = ["EfficiencyDialog"]

NO_COATING = "(perfect conductor)"


class EfficiencyDialog(QDialog):
    """Solve the previewed boundary profile and plot the orders."""

    def __init__(self, profile, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._profile = profile
        self._scan = None

        self.setWindowTitle("Efficiency of the measured groove")
        self.resize(940, 620)
        self._build_ui()
        self._show_placeholder(
            "Set a mount and press Compute.\n\n"
            f"Period {profile.period_nm:.2f} nm, measured from "
            f"{profile.n_used} averaged grooves.")

    # ── UI ───────────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QHBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        left = QWidget()
        left.setFixedWidth(280)
        column = QVBoxLayout(left)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(self._build_geometry_group())
        column.addWidget(self._build_mount_group())
        self.compute_btn = QPushButton("Compute")
        self.compute_btn.clicked.connect(self.compute)
        column.addWidget(self.compute_btn)
        self.status = QLabel("")
        self.status.setWordWrap(True)
        self.status.setStyleSheet("font-size: 11px; color: gray;")
        column.addWidget(self.status)
        column.addStretch(1)
        root.addWidget(left)

        self.figure = Figure(figsize=(6, 4.6))
        self.canvas = FigureCanvas(self.figure)
        root.addWidget(self.canvas, stretch=1)

    def _build_geometry_group(self):
        group = QGroupBox("Grating")
        form = QFormLayout(group)

        period = QLabel(f"{self._profile.period_nm:.2f} nm")
        period.setToolTip(
            "Measured from the detected groove spacing, and carried straight\n"
            "through. A .ggp export cannot hold this number.")
        form.addRow("Period (measured):", period)

        self.coating_combo = QComboBox()
        self.coating_combo.addItem(NO_COATING)
        self.coating_combo.addItems(materials.available())
        if self.coating_combo.count() > 1:
            self.coating_combo.setCurrentIndex(1)
        form.addRow("Coating:", self.coating_combo)

        # Off by default and off is the honest default: this panel has no facet
        # fit of its own, and inventing one would be worse than the fallback.
        self.blaze_check = QCheckBox("Fitted facet angle")
        self.blaze_check.setToolTip(
            "The blaze angle measured from this scan, if the Analysis tab\n"
            "produced one. Without it the reflection is evaluated on the mean\n"
            "surface, which for a sawtooth at grazing incidence is off by the\n"
            "whole blaze angle.")
        self.blaze_check.stateChanged.connect(
            lambda: self.blaze_spin.setEnabled(self.blaze_check.isChecked()))
        form.addRow(self.blaze_check)

        self.blaze_spin = QDoubleSpinBox()
        self.blaze_spin.setRange(0.1, 89.9)
        self.blaze_spin.setDecimals(2)
        self.blaze_spin.setValue(29.50)
        self.blaze_spin.setSuffix(" deg")
        self.blaze_spin.setEnabled(False)
        form.addRow("", self.blaze_spin)
        return group

    def _build_mount_group(self):
        group = QGroupBox("Mount and scan")
        form = QFormLayout(group)

        self.graze_spin = QDoubleSpinBox()
        self.graze_spin.setRange(0.05, 90.0)
        self.graze_spin.setDecimals(2)
        self.graze_spin.setValue(1.50)
        self.graze_spin.setSuffix(" deg")
        form.addRow("Graze (gamma):", self.graze_spin)

        self.azimuth_spin = QDoubleSpinBox()
        self.azimuth_spin.setRange(-89.9, 89.9)
        self.azimuth_spin.setDecimals(2)
        self.azimuth_spin.setValue(25.00)
        self.azimuth_spin.setSuffix(" deg")
        form.addRow("Azimuth (alpha):", self.azimuth_spin)

        self.start_spin = QDoubleSpinBox()
        self.start_spin.setRange(0.01, 10000.0)
        self.start_spin.setDecimals(2)
        self.start_spin.setValue(1.00)
        self.start_spin.setSuffix(" nm")
        form.addRow("Wavelength from:", self.start_spin)

        self.stop_spin = QDoubleSpinBox()
        self.stop_spin.setRange(0.02, 10000.0)
        self.stop_spin.setDecimals(2)
        self.stop_spin.setValue(5.00)
        self.stop_spin.setSuffix(" nm")
        form.addRow("to:", self.stop_spin)

        self.count_spin = QSpinBox()
        self.count_spin.setRange(2, 400)
        self.count_spin.setValue(120)
        form.addRow("Points:", self.count_spin)

        self.model_combo = QComboBox()
        self.model_combo.addItems(["local", "average", "facet"])
        self.model_combo.setToolTip(
            "How reflectivity is evaluated across the groove.\n"
            "'local' resolves it from the profile slope and is the default;\n"
            "only 'facet' consults the fitted facet angle above.")
        form.addRow("Reflectivity:", self.model_combo)
        return group

    # ── Solve ────────────────────────────────────────────────────────────────

    def _to_problem(self):
        coating = self.coating_combo.currentText()
        return self._profile.to_problem(
            coating=None if coating == NO_COATING else coating,
            blaze_angle=self.blaze_spin.value() if self.blaze_check.isChecked() else None,
        )

    def compute(self):
        """Solve and plot. Errors land in the status line, not a traceback."""
        start, stop = self.start_spin.value(), self.stop_spin.value()
        if stop <= start:
            self.status.setText("The wavelength range must increase.")
            return

        wavelengths = np.linspace(start, stop, self.count_spin.value())
        QGuiApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            self._scan = get_solver("scalar").solve(
                self._to_problem(),
                Illumination.offplane(graze=self.graze_spin.value(),
                                      azimuth=self.azimuth_spin.value()),
                wavelengths,
                reflectivity_model=self.model_combo.currentText(),
            )
        except (UnsupportedConfiguration, ValueError) as exc:
            # Both are the solver refusing rather than failing: an out-of-scope
            # configuration, or optical constants that do not cover this range.
            self._scan = None
            self.status.setText(str(exc))
            self._show_placeholder("No result")
            return
        finally:
            QGuiApplication.restoreOverrideCursor()

        self._report()
        self._draw()

    def _report(self):
        notes = self._scan.provenance.notes
        graze = notes.get("reflectivity_graze", "")
        self.status.setText(
            f"Scalar solver, not convergence-checked -- a preview, not a "
            f"defensible result.\n\n{graze}")

    # ── Draw ─────────────────────────────────────────────────────────────────

    def _show_placeholder(self, message):
        self.figure.clear()
        axes = self.figure.add_subplot(111)
        axes.text(0.5, 0.5, message, ha='center', va='center', fontsize=11,
                  color='gray', transform=axes.transAxes, wrap=True)
        axes.set_axis_off()
        self.canvas.draw()

    def _draw(self):
        scan = self._scan
        self.figure.clear()
        axes = self.figure.add_subplot(111)

        # Only orders that carry something. A grazing-incidence scan can span
        # dozens of orders, most of them evanescent over most of the range, and
        # a legend naming all of them is unreadable.
        carried = [(m, scan.order(m)) for m in scan.orders]
        carried = [(m, e) for m, e in carried if e.max() > 0.005]
        carried.sort(key=lambda pair: pair[1].max(), reverse=True)

        for m, efficiency in carried[:8]:
            axes.plot(scan.wavelengths, efficiency, lw=1.3, label=f"m = {m}")
        axes.plot(scan.wavelengths, scan.total, 'k--', lw=1.0, alpha=0.6,
                  label="total")

        axes.set_xlabel("Wavelength (nm)")
        axes.set_ylabel("Absolute efficiency")
        axes.set_title(f"Measured groove, period {self._profile.period_nm:.2f} nm")
        axes.grid(True, alpha=0.3)
        axes.set_ylim(bottom=0)
        if carried:
            axes.legend(fontsize=8, ncol=2)
        else:
            axes.text(0.5, 0.5, "No order carries more than 0.5%",
                      ha='center', va='center', color='gray',
                      transform=axes.transAxes)
        self.figure.tight_layout()
        self.canvas.draw()
