"""The scalar solver's own tab: its options, its result, its provenance.

Not built against a `SolverTab` base class -- one concrete implementation
exists and a second doesn't, and formalizing "the shape every solver tab must
have" from a sample of one would be a guess dressed as an interface.
`MainWindow` calls a small, duck-typed set of methods on whatever lives in its
`_tabs` dict (`build_options`, `set_running`, `show_progress`, `show_result`,
`show_error`, `show_field_errors`, `show_cancelled`); see the comment at
`main_window.py`'s `_build_tabs` for the documented contract. Formalize a base
class once RCWA's tab exists and the actual common shape is known rather than
assumed.
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
from ..scalar_options import (
    REFLECTIVITY_MODELS,
    ROUGHNESS_MODELS,
    VISIBILITY_MODES,
    ScalarOptionsState,
)
from ..scalar_options import build_options as build_scalar_options
from .sizing import fit_width_to_contents

__all__ = ["ScalarTab"]

#: Qt's role for stashing arbitrary data on a QListWidgetItem. Used to carry
#: the order integer without parsing it back out of the label text.
_ORDER_ROLE = Qt.ItemDataRole.UserRole


class ScalarTab(QWidget):
    """Scalar's own controls, result, and provenance.

    Holds the last scan itself (`_scan`/`_energy`), not `MainWindow` -- each
    solver tab remembers its own last result, the way switching tabs must not
    lose what the previous one computed. Geometry and the groove-profile plot
    are *not* here: both are shared across every solver and live one level up.
    """

    #: Registry key this tab represents. A future RcwaTab sets its own.
    name = "scalar"

    solve_requested = Signal()
    cancel_requested = Signal()
    #: Sweep the accuracy knob instead of solving once. A separate signal
    #: rather than a flag on `solve_requested`: the window routes it to a
    #: different worker slot, and one signal meaning two things is how a tab
    #: starts deciding what a solve *is*.
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
        # The Orders All/None/Default row sets this column's floor and it is
        # wider than it looks: 3 buttons at Qt's push-button minimum, plus
        # spacing, groupbox margins, the scrollbar and the frame. Derived, not
        # a literal -- those metrics differ across the three CI platforms.
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

        self._check_fields_match_scalar_options()

    def _build_options_group(self) -> QWidget:
        defaults = ScalarOptionsState()
        group = QGroupBox("Scalar solver")
        form = QFormLayout(group)
        quadrature = QLineEdit(defaults.quadrature_points)
        quadrature.returnPressed.connect(self.solve_requested)
        self._fields["quadrature_points"] = quadrature
        form.addRow("Quadrature pts", quadrature)

        # Both only bite with a coating named on the Setup tab, which the
        # tooltips say rather than the form guessing and disabling itself --
        # a control that greys out for reasons the user cannot see is worse
        # than one that is honestly inert.
        reflectivity = QComboBox()
        reflectivity.addItems(REFLECTIVITY_MODELS)
        reflectivity.setToolTip(
            "How reflectivity is resolved across the groove cycle.\n"
            "local: Fresnel amplitude at every quadrature point, carried\n"
            "  inside the integral, so reflectivity is order-dependent.\n"
            "average: groove-cycle mean of the intensity, one factor per\n"
            "  wavelength.\n"
            "facet: one R at the active-facet angle (the pre-M16 model,\n"
            "  kept for reproducing earlier runs).\n"
            "Only has an effect when a coating is set."
        )
        self._fields["reflectivity_model"] = reflectivity
        form.addRow("Reflectivity", reflectivity)

        visibility = QComboBox()
        visibility.addItems(VISIBILITY_MODES)
        visibility.setToolTip(
            "Which shadows the visibility masks see.\n"
            "facet-normal: the local orientation test alone -- a point is\n"
            "  shadowed iff its own facet turns away from the ray.\n"
            "horizon: adds the shadows one part of the groove casts on\n"
            "  another, for the incident and each exit direction.\n"
            "With a coating, horizon needs a per-point reflectivity model\n"
            "(local or average) -- facet has no masks for it to narrow."
        )
        self._fields["visibility"] = visibility
        form.addRow("Visibility", visibility)

        roughness = QComboBox()
        roughness.addItems(ROUGHNESS_MODELS)
        roughness.setToolTip(
            "How the surface roughness damps reflectivity.\n"
            "Nevot-Croce carries the transmitted wave and is the right one\n"
            "near the critical angle; Debye-Waller is the common\n"
            "approximation. Only has an effect with a coating and a\n"
            "non-zero roughness."
        )
        self._fields["roughness_model"] = roughness
        form.addRow("Roughness", roughness)
        return group

    def _build_orders_group(self) -> QWidget:
        """Which orders the efficiency plot draws.

        Replaces the old `if values.max() < 1e-4: continue` -- see
        `gui/orders.py` for the rule this wires up. Toggling an item redraws
        the plot from the scan already held; it never re-solves.
        """
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

        # Its own row rather than a third button beside Solve and Cancel.
        # `TestLayoutFloors` pins this column at the 276 px its Orders
        # All/None/Default row needs, and a third text-width button here would
        # push past it -- which is that test doing its job, not an obstacle.
        self._converge_button = QPushButton("Check convergence…")
        self._converge_button.setToolTip(
            "Sweep quadrature_points until the answer stops moving, and report "
            "the coarsest setting that can be defended."
        )
        self._converge_button.clicked.connect(self.convergence_requested)

        self._progress = QProgressBar()
        # Starts indeterminate and becomes determinate on the first report, so
        # a backend that declares no `reports_progress` still looks alive
        # rather than stuck at 0%.
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

    def _check_fields_match_scalar_options(self) -> None:
        """Refuse to open rather than fail later on one code path.

        `build_options` builds `ScalarOptionsState` by keyword from this
        dict, so a renamed or forgotten field is a TypeError at solve time.
        """
        declared = {f.name for f in dataclasses.fields(ScalarOptionsState)}
        if self._fields.keys() != declared:
            raise AssertionError(
                "scalar-option fields do not match ScalarOptionsState: "
                f"missing {sorted(declared - self._fields.keys())}, "
                f"extra {sorted(self._fields.keys() - declared)}"
            )

    # -- the MainWindow contract ---------------------------------------

    def build_options(self, problem, illumination, wavelengths) -> dict:
        """Validate this tab's own fields against a shared, already-parsed
        geometry. Raises `FormErrors` -- caught by `MainWindow` before
        anything reaches the worker, same as a geometry error would be."""
        options = ScalarOptionsState(
            quadrature_points=_value(self._fields["quadrature_points"]),
            reflectivity_model=_value(self._fields["reflectivity_model"]),
            roughness_model=_value(self._fields["roughness_model"]),
            visibility=_value(self._fields["visibility"]),
        )
        return build_scalar_options(problem, illumination, wavelengths, options)

    def set_running(self, running: bool) -> None:
        self._solve_button.setEnabled(not running)
        self._converge_button.setEnabled(not running)
        self._cancel_button.setEnabled(running)
        # The sweep chooses the knob, so the field it would come from is dead
        # for the duration. Leaving an ignored input live is the M9 mistake --
        # a control that looks like it does something and does not.
        self._fields["quadrature_points"].setEnabled(not running)
        if not running:
            self._progress.hide()

    def show_progress(self) -> None:
        self._progress.show()

    def show_progress_value(self, done: int, total: int) -> None:
        """Advance the bar, switching it to determinate on the first report."""
        if self._progress.maximum() != total:
            self._progress.setRange(0, total)
        self._progress.setValue(done)

    def show_solving(self, wavelength_count: int) -> None:
        # Back to indeterminate for each new run: the previous run's range is
        # meaningless here, and leaving it would show a full bar for a solve
        # that has not started.
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
        """A completed sweep, through the ordinary result path.

        No new rendering: `provenance_lines` already reads the study off
        `scan.provenance.notes["convergence"]` -- that is what M3-C put there
        -- so the panel says "convergence: yes — quadrature_points=4096 is
        enough" with no help from here.
        """
        self.show_result(scan, energy, lambda_over_period)
        if report.converged_at is not None:
            # The actionable number, put where the next plain Solve will use
            # it. It is usually far cheaper than the value the sweep had to
            # reach to prove it.
            self._fields["quadrature_points"].setText(str(report.converged_at))

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


def _value(widget: QWidget) -> str:
    """The text of an input, whichever kind it is."""
    if isinstance(widget, QComboBox):
        return widget.currentText()
    return widget.text()
