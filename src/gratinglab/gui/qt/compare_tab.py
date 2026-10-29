"""The Compare tab: one geometry, several solvers, one picture.

Every solver tab answers the problem its own way; this one asks several of
them at once and puts the answers side by side, which is the point of the
project (see :mod:`gratinglab.compare`). It holds no solver settings of its
own. Each ticked solver runs with whatever its own tab currently says, so
there is one place to set a truncation and the comparison cannot quietly use
a different one from the tab the user tuned.

It is a tab rather than a mode of the solver tabs, because a comparison has
a result of its own (the differences, and which solvers declined) that
belongs to no one solver, and because it needs a second panel that a single
solver's tab has no use for.

The layout follows the solver tabs: controls in a scrolling column on the
left, the plot above a text panel on the right. The text panel is where a
solver that refused the problem is named, with its reason; the plot simply
has no curves for it.

What is drawn, and what the panel says, is decided in
:mod:`gratinglab.gui.comparison`. Nothing here computes anything.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
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
from ..comparison import (
    ComparisonOutcome,
    comparison_lines,
    draw_comparison,
    pick_reference,
)
from ..theme import provenance_colors, style_figure
from . import theming
from .sizing import fit_width_to_contents

__all__ = ["CompareTab"]

_ORDER_ROLE = Qt.ItemDataRole.UserRole
_DIM_ROLE = Qt.ItemDataRole.UserRole + 1

#: The reference the lower panel starts from. The integral method is the
#: longest-standing rigorous solver here, and the only one that takes every
#: profile, vertical facets and undercuts included, so it is the one most
#: likely to have solved whatever else declined. Reconsidered when its
#: equal-arc default read 2.95 for R + A on the default blaze: that was the
#: mesh, not the method, and with the graded mesh the default (0.2.7) it
#: matches the C-method there to 2e-4. It is also the slowest, by far; the
#: reference only chooses what is subtracted, and pick_reference falls back
#: to another rigorous solver when it is not ticked.
_DEFAULT_REFERENCE = "integral"


class CompareTab(QWidget):
    """Tick solvers, run them on the shared geometry, see them together."""

    #: Not a registry name. The window uses it to route the worker's
    #: replies here, and `MainWindow.tabs` deliberately does not contain it.
    name = "compare"

    compare_requested = Signal()
    cancel_requested = Signal()

    def __init__(self, solvers, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        #: Every registered solver with a tab, in tab order. Line styles are
        #: fixed by position in this tuple, so unticking one never restyles
        #: the rest.
        self._solvers: tuple[str, ...] = tuple(solvers)
        self._outcome: ComparisonOutcome | None = None
        self._lines: tuple = ()
        self._order_colors: dict[int, str] = {}
        self._visible_orders: frozenset[int] = frozenset()
        self._previous_orders: frozenset[int] = frozenset()
        self._build()

    # -- construction ----------------------------------------------------

    def _build(self) -> None:
        controls = QWidget()
        column = QVBoxLayout(controls)
        column.setContentsMargins(0, 0, 0, 0)
        column.addWidget(self._build_solvers_group())
        column.addWidget(self._build_orders_group())
        column.addLayout(self._build_action_buttons())
        column.addWidget(self._progress)
        column.addStretch(1)

        scroller = QScrollArea()
        scroller.setWidget(controls)
        scroller.setWidgetResizable(True)
        fit_width_to_contents(scroller)

        right = QSplitter(Qt.Orientation.Vertical)
        right.addWidget(self._build_plot())
        right.addWidget(self._build_panel())
        right.setStretchFactor(0, 1)
        right.setStretchFactor(1, 0)
        right.setSizes([560, 160])

        body = QSplitter(Qt.Orientation.Horizontal)
        body.addWidget(scroller)
        body.addWidget(right)
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setSizes([260, 700])

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(body)

    def _build_solvers_group(self) -> QWidget:
        from ..docs import display_title, short_title

        group = QGroupBox("Solvers")
        layout = QVBoxLayout(group)
        note = QLabel("Each runs with the settings on its own tab.")
        note.setWordWrap(True)
        layout.addWidget(note)

        self._checks: dict[str, QCheckBox] = {}
        for name in self._solvers:
            box = QCheckBox(short_title(name))
            box.setToolTip(display_title(name))
            box.setChecked(True)
            box.toggled.connect(self._on_solvers_ticked)
            self._checks[name] = box
            layout.addWidget(box)

        form = QFormLayout()
        self._reference = QComboBox()
        for name in self._solvers:
            self._reference.addItem(short_title(name), name)
        start = self._reference.findData(_DEFAULT_REFERENCE)
        self._reference.setCurrentIndex(max(start, 0))
        self._reference.setToolTip(
            "The lower panel shows every other solver minus this one. "
            "Changing it redraws without solving again."
        )
        self._reference.currentIndexChanged.connect(self._redraw)
        form.addRow("Reference", self._reference)

        self._log_difference = QCheckBox("log |Δ|")
        self._log_difference.setToolTip(
            "Plot the size of the difference on a log scale, so an agreement\n"
            "to 1e-7 between rigorous solvers is visible beside a scalar\n"
            "result that is off by 0.1."
        )
        self._log_difference.toggled.connect(self._redraw)
        form.addRow("Difference", self._log_difference)
        layout.addLayout(form)
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

    def _build_action_buttons(self) -> QHBoxLayout:
        self._compare_button = QPushButton("Compare")
        self._compare_button.setToolTip(
            "Solve the geometry with every ticked solver, one after another."
        )
        self._compare_button.clicked.connect(self.compare_requested)
        self._cancel_button = QPushButton("Cancel")
        self._cancel_button.clicked.connect(self.cancel_requested)
        self._cancel_button.setEnabled(False)

        self._progress = QProgressBar()
        self._progress.setRange(0, 0)
        self._progress.hide()

        row = QHBoxLayout()
        row.addWidget(self._compare_button)
        row.addWidget(self._cancel_button)
        return row

    def _build_plot(self) -> QWidget:
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg,
            NavigationToolbar2QT,
        )
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(8, 5), layout="constrained")
        grid = self._figure.add_gridspec(2, 1, height_ratios=(3, 2))
        self._axes = self._figure.add_subplot(grid[0])
        self._difference_axes = self._figure.add_subplot(grid[1], sharex=self._axes)
        self._canvas = FigureCanvasQTAgg(self._figure)
        style_figure(self._figure, theming.current())
        theming.subscribe(self._on_theme)

        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(NavigationToolbar2QT(self._canvas, panel))
        layout.addWidget(self._canvas)
        return panel

    def _build_panel(self) -> QWidget:
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(8, 4, 8, 8)
        heading = QLabel("Comparison")
        # Bold by font, not stylesheet, as the solver tabs' heading is.
        bold = heading.font()
        bold.setBold(True)
        heading.setFont(bold)
        self._panel = QTextEdit()
        self._panel.setReadOnly(True)
        layout.addWidget(heading)
        layout.addWidget(self._panel)
        self._paint((provenance.Line(
            "Tick the solvers to compare and press Compare. The rigorous "
            "ones can take minutes on a long scan; Cancel stops them.\n", "dim",
        ),))
        return panel

    # -- what the window asks ----------------------------------------------

    def selected(self) -> tuple[str, ...]:
        """The ticked solvers, in tab order."""
        return tuple(n for n in self._solvers if self._checks[n].isChecked())

    def reference(self) -> str | None:
        return self._reference.currentData()

    def set_selected(self, names) -> None:
        wanted = set(names)
        for name, box in self._checks.items():
            box.setChecked(name in wanted)

    def set_reference(self, name: str) -> None:
        index = self._reference.findData(name)
        if index < 0:
            raise KeyError(name)
        self._reference.setCurrentIndex(index)

    def set_running(self, running: bool) -> None:
        self._compare_button.setEnabled(not running and bool(self.selected()))
        self._cancel_button.setEnabled(running)
        for box in self._checks.values():
            box.setEnabled(not running)
        if not running:
            self._progress.hide()

    def show_progress(self) -> None:
        self._progress.show()

    def show_progress_value(self, done: int, total: int) -> None:
        if self._progress.maximum() != total:
            self._progress.setRange(0, total)
        self._progress.setValue(done)

    def show_solving(self, names, wavelength_count: int) -> None:
        from ..docs import short_title

        self._progress.setRange(0, 0)
        self._paint((provenance.Line(
            f"comparing: {', '.join(short_title(n) for n in names)}, "
            f"{wavelength_count} wavelengths each…\n", "dim",
        ),))

    def show_outcome(self, outcome: ComparisonOutcome) -> None:
        self._outcome = outcome
        comparison = outcome.comparison
        if comparison is not None:
            summaries = _summaries(comparison)
            self._visible_orders = orders_module.carry_over(
                self._visible_orders, self._previous_orders, summaries
            )
            self._previous_orders = frozenset(s.order for s in summaries)
            self._rebuild_order_list(summaries)
        self._redraw()

    def show_cancelled(self) -> None:
        if self._outcome is None:
            self._paint((provenance.Line(
                "cancelled: the comparison stopped before any result\n", "dim"
            ),))
            return
        self._paint(comparison_lines(self._outcome, self.reference(), cancelled=True))

    def show_error(self, message: str) -> None:
        self._paint((provenance.Line(f"comparison failed: {message}\n", "bad"),))

    def show_field_errors(self, errors) -> None:
        self._paint(provenance.error_lines(errors))

    # -- orders ----------------------------------------------------------

    def _rebuild_order_list(self, summaries) -> None:
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
        self._update_order_count_label()

    def _apply_swatches(self) -> None:
        from PySide6.QtGui import QColor

        theme = theming.current()
        muted, text = QColor(theme.muted), QColor(theme.text)
        self._order_list.blockSignals(True)
        try:
            for i in range(self._order_list.count()):
                item = self._order_list.item(i)
                color = self._order_colors.get(item.data(_ORDER_ROLE))
                item.setIcon(theming.swatch(color) if color else theming.QIcon())
                item.setForeground(muted if item.data(_DIM_ROLE) else text)
        finally:
            self._order_list.blockSignals(False)

    def _update_order_count_label(self) -> None:
        self._order_count_label.setText(
            f"{len(self._visible_orders)} of {self._order_list.count()} shown"
        )

    def _on_order_item_changed(self, item: QListWidgetItem) -> None:
        order = item.data(_ORDER_ROLE)
        if item.checkState() == Qt.CheckState.Checked:
            self._visible_orders = self._visible_orders | {order}
        else:
            self._visible_orders = self._visible_orders - {order}
        self._update_order_count_label()
        self._redraw()

    def _show_all_orders(self) -> None:
        self._visible_orders = frozenset(self._previous_orders)
        self._resync_order_checkboxes()

    def _show_no_orders(self) -> None:
        self._visible_orders = frozenset()
        self._resync_order_checkboxes()

    def _show_default_orders(self) -> None:
        if self._outcome is None or self._outcome.comparison is None:
            return
        self._visible_orders = orders_module.default_visible(
            _summaries(self._outcome.comparison)
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
        self._update_order_count_label()
        self._redraw()

    # -- drawing ---------------------------------------------------------

    def _on_solvers_ticked(self, _checked: bool = False) -> None:
        # Nothing ticked means nothing to run; say so with the button rather
        # than with an error after it is pressed.
        if self._cancel_button.isEnabled():
            return  # running; set_running(False) re-evaluates
        self._compare_button.setEnabled(bool(self.selected()))

    def _paint(self, lines) -> None:
        self._lines = tuple(lines)
        self._panel.setHtml(
            provenance.to_html(self._lines, provenance_colors(theming.current()))
        )

    def _on_theme(self, theme) -> None:
        style_figure(self._figure, theme)
        self._redraw()
        if self._lines:
            self._paint(self._lines)

    def _redraw(self, *_args) -> None:
        outcome = self._outcome
        if outcome is None:
            self._canvas.draw_idle()
            return
        comparison = outcome.comparison
        if comparison is not None:
            reference = pick_reference(self.reference(), comparison.methods)
            self._order_colors = draw_comparison(
                self._axes, self._difference_axes, comparison, reference,
                self._visible_orders, theming.current(),
                solvers=self._solvers,
                log_difference=self._log_difference.isChecked(),
            )
        else:
            self._order_colors = {}
            for axes in (self._axes, self._difference_axes):
                axes.clear()
            style_figure(self._figure, theming.current())
        self._apply_swatches()
        self._paint(comparison_lines(outcome, self.reference()))
        self._canvas.draw_idle()


def _summaries(comparison) -> tuple["orders_module.OrderSummary", ...]:
    """One summary per order over every solver, for the order list.

    The peak is the largest any solver reached, so an order only one method
    puts power into is still offered by default.
    """
    import numpy as np

    peaks = comparison.efficiency.max(axis=(0, 1))
    propagating = np.zeros(len(comparison.orders), dtype=bool)
    for scan in comparison.scans:
        columns = np.searchsorted(comparison.orders, scan.orders)
        propagating[columns] |= scan.propagating.any(axis=0)
    return tuple(
        orders_module.OrderSummary(
            order=int(order),
            column=column,
            peak=float(peaks[column]),
            ever_propagating=bool(propagating[column]),
        )
        for column, order in enumerate(comparison.orders)
    )
