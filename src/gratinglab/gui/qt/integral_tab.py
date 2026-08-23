"""The integral solver's own tab: its knob, its result, its provenance.

Deliberately a parallel construction of ``scalar_tab.py`` rather than a
shared base class: the roadmap's own warning is that every solver seam was
designed against a sample of one, and the way to learn the actual common
shape is to build the second tab and *then* extract what proved common. The
duplication between these two files is that evidence, kept visible.

What genuinely differs from scalar: the options group (one accuracy knob,
``boundary_points``, no reflectivity or roughness -- the perfectly conducting
model has no material in it) and the state type behind ``build_options``.
"""

from __future__ import annotations

import csv
import dataclasses
from pathlib import Path
from typing import Sequence

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .. import orders as orders_module
from .. import provenance
from ..integral_options import IntegralOptionsState
from ..integral_options import build_options as build_integral_options
from .sizing import fit_width_to_contents

__all__ = ["IntegralTab"]

_ORDER_ROLE = Qt.ItemDataRole.UserRole


class IntegralTab(QWidget):
    """Integral's own controls, result, and provenance."""

    name = "integral"

    solve_requested = Signal()
    cancel_requested = Signal()
    convergence_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._fields: dict[str, QWidget] = {}
        self._scan = None
        self._energy = None
        self._lambda_over_period: float | None = None
        self._visible_orders: frozenset[int] = frozenset()
        self._previous_orders: frozenset[int] = frozenset()
        self._build()

    # -- construction ------------------------------------------------

    def _build(self) -> None:
        controls = QWidget()
        column = QVBoxLayout(controls)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(self._build_options_group())
        column.addWidget(self._build_orders_group())
        column.addLayout(self._build_action_buttons())
        column.addWidget(self._export_button)
        column.addWidget(self._progress)
        column.addStretch(1)

        scroller = QScrollArea()
        scroller.setWidget(controls)
        scroller.setWidgetResizable(True)
        fit_width_to_contents(scroller)

        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(self._build_efficiency_plot())
        right.addWidget(self._build_provenance())
        right.setStretchFactor(0, 1)
        right.setStretchFactor(1, 0)

        body = QSplitter(Qt.Orientation.Horizontal)
        body.addWidget(scroller)
        body.addWidget(right)
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setSizes([260, 700])

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(body)

        self._check_fields_match_integral_options()

    def _build_options_group(self) -> QWidget:
        defaults = IntegralOptionsState()
        group = QGroupBox("Integral solver")
        form = QFormLayout(group)
        points = QLineEdit(defaults.boundary_points)
        points.returnPressed.connect(self.solve_requested)
        points.setToolTip(
            "Equal-arc-length nodes on one groove period -- the accuracy\n"
            "knob the convergence sweep varies. Cost grows roughly cubically."
        )
        self._fields["boundary_points"] = points
        form.addRow("Boundary pts", points)

        note = QLabel(
            "Perfectly conducting boundary: efficiencies are relative to a\n"
            "perfect reflector and sum to 1. A coating named on the Setup\n"
            "tab is not consulted (the provenance panel says so)."
        )
        note.setWordWrap(True)
        form.addRow(note)
        return group

    def _build_orders_group(self) -> QWidget:
        group = QGroupBox("Orders")
        layout = QVBoxLayout(group)

        self._order_list = QListWidget()
        self._order_list.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self._order_list.itemChanged.connect(self._on_order_item_changed)
        layout.addWidget(self._order_list)

        buttons = QHBoxLayout()
        show_all = QPushButton("All")
        show_all.clicked.connect(self._show_all_orders)
        show_none = QPushButton("None")
        show_none.clicked.connect(self._show_no_orders)
        show_default = QPushButton("Default")
        show_default.clicked.connect(self._show_default_orders)
        for button in (show_all, show_none, show_default):
            buttons.addWidget(button)
        layout.addLayout(buttons)

        self._order_count_label = QLabel("")
        layout.addWidget(self._order_count_label)
        return group

    def _build_action_buttons(self) -> QVBoxLayout:
        self._solve_button = QPushButton("Solve")
        self._solve_button.clicked.connect(self.solve_requested)
        self._cancel_button = QPushButton("Cancel")
        self._cancel_button.clicked.connect(self.cancel_requested)
        self._cancel_button.setEnabled(False)
        self._export_button = QPushButton("Export CSV…")
        self._export_button.clicked.connect(self.export_csv)

        self._converge_button = QPushButton("Check convergence…")
        self._converge_button.setToolTip(
            "Sweep boundary_points until the answer stops moving, and report "
            "the coarsest setting that can be defended."
        )
        self._converge_button.clicked.connect(self.convergence_requested)

        self._progress = QProgressBar()
        self._progress.setRange(0, 0)
        self._progress.hide()

        row = QHBoxLayout()
        row.addWidget(self._solve_button)
        row.addWidget(self._cancel_button)

        column = QVBoxLayout()
        column.addLayout(row)
        column.addWidget(self._converge_button)
        return column

    def _build_efficiency_plot(self) -> QWidget:
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg,
            NavigationToolbar2QT,
        )
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(8, 3), layout="constrained")
        self._axes = self._figure.add_subplot(1, 1, 1)
        self._canvas = FigureCanvasQTAgg(self._figure)

        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(NavigationToolbar2QT(self._canvas, panel))
        layout.addWidget(self._canvas)
        return panel

    def _build_provenance(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 4, 8, 8)
        heading = QLabel("Provenance")
        heading.setStyleSheet("font-weight: bold")
        self._provenance = QTextEdit()
        self._provenance.setReadOnly(True)
        layout.addWidget(heading)
        layout.addWidget(self._provenance)
        return panel

    def _check_fields_match_integral_options(self) -> None:
        declared = {f.name for f in dataclasses.fields(IntegralOptionsState)}
        if self._fields.keys() != declared:
            raise AssertionError(
                "integral-option fields do not match IntegralOptionsState: "
                f"missing {sorted(declared - self._fields.keys())}, "
                f"extra {sorted(self._fields.keys() - declared)}"
            )

    # -- the MainWindow contract ---------------------------------------

    def build_options(self, problem, illumination, wavelengths) -> dict:
        options = IntegralOptionsState(
            boundary_points=self._fields["boundary_points"].text(),
        )
        return build_integral_options(problem, illumination, wavelengths, options)

    def set_running(self, running: bool) -> None:
        self._solve_button.setEnabled(not running)
        self._converge_button.setEnabled(not running)
        self._cancel_button.setEnabled(running)
        self._fields["boundary_points"].setEnabled(not running)
        if not running:
            self._progress.hide()

    def show_progress(self) -> None:
        self._progress.show()

    def show_progress_value(self, done: int, total: int) -> None:
        if self._progress.maximum() != total:
            self._progress.setRange(0, total)
        self._progress.setValue(done)

    def show_solving(self, wavelength_count: int) -> None:
        self._progress.setRange(0, 0)
        self._paint(provenance.solving_lines(self.name, wavelength_count))

    def show_result(self, scan, energy, lambda_over_period: float) -> None:
        self._scan, self._energy = scan, energy
        self._lambda_over_period = lambda_over_period

        summaries = orders_module.summarize(scan)
        self._visible_orders = orders_module.carry_over(
            self._visible_orders, self._previous_orders, summaries
        )
        self._previous_orders = frozenset(s.order for s in summaries)
        self._rebuild_order_list(summaries)

        self._draw_efficiency(scan)
        self._paint(provenance.provenance_lines(scan, energy, lambda_over_period))

    def show_convergence(self, scan, energy, lambda_over_period: float, report) -> None:
        self.show_result(scan, energy, lambda_over_period)
        if report.converged_at is not None:
            self._fields["boundary_points"].setText(str(report.converged_at))

    def show_cancelled(self) -> None:
        if self._scan is None:
            return
        self._paint(
            provenance.provenance_lines(
                self._scan, self._energy, self._lambda_over_period, cancelled=True
            )
        )

    def show_error(self, message: str) -> None:
        self._paint((provenance.Line(f"solve failed: {message}\n", "bad"),))

    def show_field_errors(self, errors) -> None:
        self._paint(provenance.error_lines(errors))

    # -- orders ----------------------------------------------------------

    def _rebuild_order_list(self, summaries: Sequence["orders_module.OrderSummary"]) -> None:
        self._order_list.blockSignals(True)
        self._order_list.clear()
        for summary in summaries:
            item = QListWidgetItem(orders_module.describe(summary))
            item.setData(_ORDER_ROLE, summary.order)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(
                Qt.CheckState.Checked
                if summary.order in self._visible_orders
                else Qt.CheckState.Unchecked
            )
            if not summary.ever_propagating:
                item.setForeground(self.palette().mid().color())
            self._order_list.addItem(item)
        self._order_list.blockSignals(False)
        self._update_order_count_label(len(summaries))

    def _update_order_count_label(self, total: int) -> None:
        self._order_count_label.setText(
            f"{len(self._visible_orders)} of {total} shown"
        )

    def _on_order_item_changed(self, item: QListWidgetItem) -> None:
        order = item.data(_ORDER_ROLE)
        if item.checkState() == Qt.CheckState.Checked:
            self._visible_orders = self._visible_orders | {order}
        else:
            self._visible_orders = self._visible_orders - {order}
        self._update_order_count_label(self._order_list.count())
        if self._scan is not None:
            self._draw_efficiency(self._scan)

    def _show_all_orders(self) -> None:
        self._visible_orders = frozenset(self._previous_orders)
        self._resync_order_checkboxes()

    def _show_no_orders(self) -> None:
        self._visible_orders = frozenset()
        self._resync_order_checkboxes()

    def _show_default_orders(self) -> None:
        if self._scan is None:
            return
        self._visible_orders = orders_module.default_visible(
            orders_module.summarize(self._scan)
        )
        self._resync_order_checkboxes()

    def _resync_order_checkboxes(self) -> None:
        self._order_list.blockSignals(True)
        for i in range(self._order_list.count()):
            item = self._order_list.item(i)
            item.setCheckState(
                Qt.CheckState.Checked
                if item.data(_ORDER_ROLE) in self._visible_orders
                else Qt.CheckState.Unchecked
            )
        self._order_list.blockSignals(False)
        self._update_order_count_label(self._order_list.count())
        if self._scan is not None:
            self._draw_efficiency(self._scan)

    # -- drawing -----------------------------------------------------

    def _paint(self, lines) -> None:
        self._provenance.setHtml(provenance.to_html(lines))

    def _draw_efficiency(self, scan) -> None:
        axes = self._axes
        axes.clear()
        for index, order in enumerate(scan.orders):
            if int(order) not in self._visible_orders:
                continue
            axes.plot(
                scan.wavelengths, scan.efficiency[:, index],
                lw=1.3, label=f"m={int(order):+d}",
            )
        axes.plot(scan.wavelengths, scan.total, "k--", lw=1.0, alpha=0.7, label="Σ")
        axes.axhline(1.0, color="#b00020", lw=0.8, ls=":", alpha=0.7)
        axes.set_xlabel("wavelength (nm)")
        axes.set_ylabel("efficiency")
        axes.set_ylim(bottom=0)
        axes.grid(alpha=0.25)
        handles, _ = axes.get_legend_handles_labels()
        if handles:
            axes.legend(fontsize=7, ncol=max(1, len(handles) // 8), loc="upper right")
        self._canvas.draw_idle()

    # -- export ------------------------------------------------------

    def export_csv(self) -> None:
        if self._scan is None:
            QMessageBox.information(self, "Nothing to export", "Solve first.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export efficiencies", "efficiency.csv", "CSV (*.csv)"
        )
        if not path:
            return
        rows = self._scan.to_records()
        with open(path, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        QMessageBox.information(
            self, "Exported", f"{len(rows)} rows to {Path(path).name}"
        )
