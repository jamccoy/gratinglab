"""The RCWA solver's own tab: its knob, its result, its provenance.

The third solver tab, built on the same skeleton as ``integral_tab.py`` and
``scalar_tab.py`` rather than on a shared base class, for the reason the
integral tab's docstring gives: the common structure is to be extracted once
the copies show what is common, and three copies now say more than two did.

What differs is the options group. RCWA has two settings: ``truncation``,
the harmonic window ``-M..M`` the convergence sweep varies, and ``slices``,
the staircase's slab count for a non-lamellar profile (blank follows the
truncation). Beneath them a note says what the solver will do with the
current profile, exact slab or staircase, because that is the one fact about
an RCWA result most worth seeing before trusting it.
"""

from __future__ import annotations

import csv
import dataclasses
from pathlib import Path
from typing import Sequence

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
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
from ..efficiency_plot import draw_efficiency
from ..theme import provenance_colors, style_figure
from . import theming
from ..rcwa_options import RcwaOptionsState, profile_note
from ..rcwa_options import build_options as build_rcwa_options
from .sizing import fit_width_to_contents

__all__ = ["RcwaTab"]

_ORDER_ROLE = Qt.ItemDataRole.UserRole
#: Set on an order that never propagates in the scan, drawn in the muted colour.
_DIM_ROLE = Qt.ItemDataRole.UserRole + 1


class RcwaTab(QWidget):
    """RCWA's own controls, result, and provenance."""

    name = "rcwa"

    solve_requested = Signal()
    cancel_requested = Signal()
    convergence_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._fields: dict[str, QWidget] = {}
        self._scan = None
        self._lines: tuple = ()
        self._order_colors: dict[int, str] = {}
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

        self._check_fields_match_rcwa_options()

    def _build_options_group(self) -> QWidget:
        defaults = RcwaOptionsState()
        group = QGroupBox("RCWA solver")
        form = QFormLayout(group)
        truncation = QLineEdit(defaults.truncation)
        truncation.returnPressed.connect(self.solve_requested)
        truncation.setToolTip(
            "M: the harmonics -M..M kept in the Fourier expansion -- the\n"
            "accuracy knob the convergence sweep varies. It must hold every\n"
            "order the scan reports, and every order a lossless substrate\n"
            "transmits. Cost grows roughly as M^3 per slab."
        )
        self._fields["truncation"] = truncation
        form.addRow("Truncation M", truncation)

        slices = QLineEdit(defaults.slices)
        slices.setPlaceholderText("auto")
        slices.returnPressed.connect(self.solve_requested)
        slices.setToolTip(
            "Slabs in the staircase of a profile that is not lamellar.\n"
            "Blank follows the truncation (max(8, 2M)), so one sweep refines\n"
            "both; a number pins it. A lamellar groove is one exact slab and\n"
            "ignores this."
        )
        self._fields["slices"] = slices
        form.addRow("Slices", slices)

        material = QLabel(
            "Needs a coating in the geometry panel: the grating is cut into "
            "it, and it fills the substrate. Efficiencies are absolute; a "
            "lossless substrate's transmitted orders are recorded too."
        )
        material.setWordWrap(True)
        form.addRow(material)

        # Filled in when the geometry is known (at each solve), because what
        # it says depends on the profile, not on anything in this group.
        self._profile_note = QLabel()
        self._profile_note.setWordWrap(True)
        form.addRow(self._profile_note)
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
            "Sweep truncation until the answer stops moving, and report "
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
        style_figure(self._figure, theming.current())
        theming.subscribe(self._on_theme)

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
        # Bold by font, not stylesheet: a stylesheet freezes the label's
        # palette at polish time, and it would keep the old theme's colour.
        bold = heading.font()
        bold.setBold(True)
        heading.setFont(bold)
        self._provenance = QTextEdit()
        self._provenance.setReadOnly(True)
        layout.addWidget(heading)
        layout.addWidget(self._provenance)
        return panel

    def _check_fields_match_rcwa_options(self) -> None:
        declared = {f.name for f in dataclasses.fields(RcwaOptionsState)}
        if self._fields.keys() != declared:
            raise AssertionError(
                "rcwa-option fields do not match RcwaOptionsState: "
                f"missing {sorted(declared - self._fields.keys())}, "
                f"extra {sorted(self._fields.keys() - declared)}"
            )

    # -- the MainWindow contract ---------------------------------------

    def build_options(self, problem, illumination, wavelengths) -> dict:
        options = RcwaOptionsState(
            truncation=_value(self._fields["truncation"]),
            slices=_value(self._fields["slices"]),
        )
        self._profile_note.setText(
            profile_note(problem, options.slices, options.truncation)
        )
        return build_rcwa_options(problem, illumination, wavelengths, options)

    def set_running(self, running: bool) -> None:
        self._solve_button.setEnabled(not running)
        self._converge_button.setEnabled(not running)
        self._cancel_button.setEnabled(running)
        # Only the knob: the sweep chooses it, so the field it would come from
        # is dead for the duration. `slices` is passed through untouched.
        self._fields["truncation"].setEnabled(not running)
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
            self._fields["truncation"].setText(str(report.converged_at))

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
            item.setData(_DIM_ROLE, not summary.ever_propagating)
            self._order_list.addItem(item)
        self._order_list.blockSignals(False)
        self._apply_swatches()
        self._update_order_count_label(len(summaries))

    def _apply_swatches(self) -> None:
        """Give each listed order its line's colour, and dim the dead ones.

        Both are redone on a theme change: the swatch because the colours
        move, the dimming because it was the old palette's muted colour.
        """
        # From the theme, not self.palette(): on a theme change this runs
        # before Qt has propagated the new application palette to the widget.
        from PySide6.QtGui import QColor

        theme = theming.current()
        muted, text = QColor(theme.muted), QColor(theme.text)
        # Restyling an item emits itemChanged, whose handler redraws, which
        # restyles again: blocked, as the checkbox resync already is.
        self._order_list.blockSignals(True)
        try:
            for i in range(self._order_list.count()):
                item = self._order_list.item(i)
                color = self._order_colors.get(item.data(_ORDER_ROLE))
                item.setIcon(theming.swatch(color) if color else theming.QIcon())
                item.setForeground(muted if item.data(_DIM_ROLE) else text)
        finally:
            self._order_list.blockSignals(False)

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
        self._lines = tuple(lines)
        self._provenance.setHtml(
            provenance.to_html(self._lines, provenance_colors(theming.current()))
        )

    def _on_theme(self, theme) -> None:
        style_figure(self._figure, theme)
        if self._scan is not None:
            self._draw_efficiency(self._scan)
        else:
            self._canvas.draw_idle()
        if self._lines:
            self._paint(self._lines)
        self._apply_swatches()

    def _draw_efficiency(self, scan) -> None:
        self._order_colors = draw_efficiency(
            self._axes, scan, self._visible_orders, theming.current()
        )
        self._apply_swatches()
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
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        QMessageBox.information(
            self, "Exported", f"{len(rows)} rows to {Path(path).name}"
        )


def _value(widget: QWidget) -> str:
    """The text of an input, whichever kind it is."""
    if isinstance(widget, QComboBox):
        return widget.currentText()
    return widget.text()
