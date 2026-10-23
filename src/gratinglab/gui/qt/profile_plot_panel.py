"""The groove-profile plot: height vs. position.

The plot is drawn from `Problem.height_nm(...)` alone. No solver is involved
in computing the shape of a groove, so this plot is shared by every solver a
tab might run instead of being duplicated in each. It is redrawn once per
completed solve (from whichever solver has just finished, since the geometry
is the same in either case), not once per tab.
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

        Two periods are shown, as in the Tk window and the first Qt window,
        because a single period does not show that the profile repeats.
        """
        import numpy as np

        from .. import diagram

        axes = self._axes
        axes.clear()
        # Physical nm along d-hat, not the raw profile parameter. `t` runs
        # against the dispersion direction (conventions.md §3), so plotting
        # against `t` mirrors the groove -- and this panel sits beside a
        # diagram that draws the same object with rays on it. Two plots of one
        # grating with opposite handedness would be a small lie.
        t = np.linspace(-2.0, 0.0, 1200)
        x = diagram.x_nm(t, parsed.problem.period)
        height = parsed.problem.height_nm(t)
        axes.plot(x, height, color="#1f3b73", lw=1.8)
        axes.fill_between(x, 0, height, color="#1f3b73", alpha=0.12)
        axes.set_xlabel("position along d̂ (nm), two periods shown")
        axes.set_ylabel("height (nm)")
        axes.set_title(
            f"{type(parsed.problem.profile).__name__}: "
            f"period {parsed.problem.period:g} nm, depth {parsed.problem.depth:.4g} nm",
            fontsize=10,
        )
        axes.grid(alpha=0.25)
        self._canvas.draw_idle()
