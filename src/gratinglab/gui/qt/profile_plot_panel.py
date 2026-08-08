"""The groove-profile plot: height vs. position.

Drawn from `Problem.height_nm(...)` alone -- no solver is involved in
computing a groove's shape, so this plot is shared across every solver a tab
might run rather than duplicated inside each one. It is redrawn once per
completed solve (from whichever solver just finished, since the geometry is
the same either way), not once per tab.
"""

from __future__ import annotations

from PySide6.QtWidgets import QVBoxLayout, QWidget

__all__ = ["ProfilePlotPanel"]


class ProfilePlotPanel(QWidget):
    """A single matplotlib axes, embedded with its own navigation toolbar."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        from matplotlib.backends.backend_qtagg import (
            FigureCanvasQTAgg,
            NavigationToolbar2QT,
        )
        from matplotlib.figure import Figure

        self._figure = Figure(figsize=(8, 3), layout="constrained")
        self._axes = self._figure.add_subplot(1, 1, 1)
        self._canvas = FigureCanvasQTAgg(self._figure)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(NavigationToolbar2QT(self._canvas, self))
        layout.addWidget(self._canvas)

    def draw(self, parsed) -> None:
        """Redraw from a solved geometry.

        Two periods shown, exactly as the Tk and first-light Qt windows drew
        it -- one period alone does not show that the profile repeats.
        """
        import numpy as np

        axes = self._axes
        axes.clear()
        t = np.linspace(0.0, 1.0, 600, endpoint=False)
        two_periods = np.concatenate([t, t + 1.0])
        height = parsed.problem.height_nm(two_periods)
        axes.plot(two_periods, height, color="#1f3b73", lw=1.8)
        axes.fill_between(two_periods, 0, height, color="#1f3b73", alpha=0.12)
        axes.set_xlabel("position / period  (two periods shown)")
        axes.set_ylabel("height (nm)")
        axes.set_title(
            f"{type(parsed.problem.profile).__name__} — "
            f"period {parsed.problem.period:g} nm, depth {parsed.problem.depth:.4g} nm",
            fontsize=10,
        )
        axes.grid(alpha=0.25)
        self._canvas.draw_idle()
