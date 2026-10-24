"""The scalar solver's own tab: its options, its result, its provenance.

This tab is not built on a `SolverTab` base class. It was written when it was
the only solver tab, and defining the structure that every solver tab must have
from a sample of one would have been a guess presented as an interface. The
integral tab has since followed the same small set of methods without one. `MainWindow` calls a small, duck-typed set of methods on
whatever lives in its `_tabs` dict (`build_options`, `set_running`,
`show_progress`, `show_result`, `show_error`, `show_field_errors`,
`show_cancelled`); see the comment at `_build_tabs` in `main_window.py` for
the documented contract. A base class should be formalized once RCWA's tab
exists and the common structure is known from evidence.
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
#: Set on an order that never propagates in the scan, drawn in the muted colour.
_DIM_ROLE = Qt.ItemDataRole.UserRole + 1


class ScalarTab(QWidget):
    """Scalar's own controls, result, and provenance.

    The tab holds the last scan itself (`_scan`/`_energy`), instead of
    `MainWindow` holding it. Each solver tab remembers its own last result, so
    that switching tabs does not lose what the previous one computed. The
    geometry and the groove-profile plot are not held here, because both are
    shared by every solver and live one level up.
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

        This replaces the old `if values.max() < 1e-4: continue`; see
        `gui/orders.py` for the rule it connects to. Toggling an item redraws
        the plot from the scan already held and never re-solves.
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

    def _check_fields_match_scalar_options(self) -> None:
        """Refuse to open instead of failing later on one code path.

        `build_options` builds `ScalarOptionsState` by keyword from this
        dict, so a renamed or forgotten field would otherwise cause a
        TypeError at solve time.
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
        """Validate this tab's own fields against a shared, already-parsed geometry.

        Raises `FormErrors`, which `MainWindow` catches before anything
        reaches the worker, in the same way as a geometry error."""
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
        """Show a completed sweep through the ordinary result path.

        No new rendering is needed. `provenance_lines` already reads the study
        from `scan.provenance.notes["convergence"]`, where M3-C placed it, so
        the panel says "convergence: yes; quadrature_points=4096 is
        enough" without any help from this method.
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
