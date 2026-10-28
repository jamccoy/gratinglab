"""
The matplotlib canvas and every plot the window can draw.

All of these run on the main thread. The worker hands back result dictionaries;
nothing here is ever called from another thread.
"""
from __future__ import annotations

import functools

import numpy as np

from . import *  # noqa: F401,F403 - sets QT_API before matplotlib's shim loads

import matplotlib
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from gratinglab.gui.qt import theming
from gratinglab.gui.theme import set_titles, style_figure

__all__ = ["PlotCanvas"]


def _remembered(method):
    """Record the call, so a theme change can draw the same plot again."""
    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        self._last = (method, args, kwargs)
        return method(self, *args, **kwargs)
    return wrapper


class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(8, 5), tight_layout=True)
        super().__init__(self.fig)
        self.setParent(parent)
        self._last = None
        theming.subscribe(self._on_theme)
        self.show_placeholder("Load an AFM file to see the raw profile")

    def _on_theme(self, theme):
        """Redraw the last plot in the new theme.

        Every plot clears the figure and builds fresh axes, which take their
        colours from matplotlib's settings at that moment, so drawing it again
        is the whole restyle. The figure itself predates the change and needs
        its own background set first.
        """
        style_figure(self.fig, theme)
        if self._last is not None:
            method, args, kwargs = self._last
            method(self, *args, **kwargs)
        else:
            self.draw_idle()

    @_remembered
    def show_placeholder(self, message):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=12,
                color=theming.current().muted, transform=ax.transAxes)
        ax.set_axis_off()
        self.draw()

    @_remembered
    def plot_raw(self, disp_um, profile_nm, filename, scan_size_um):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.plot(disp_um, profile_nm, color='steelblue', linewidth=1.0)
        ax.set_xlabel("Position (µm)", fontsize=11)
        ax.set_ylabel("Height (nm)", fontsize=11)
        set_titles(ax, "Raw profile", f"scan width {scan_size_um:.3f} µm",
                   theming.current())
        ax.grid(True, alpha=0.3)
        self.draw()

    @_remembered
    def plot_2d(self, data, scan_size_um, filename):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        im = ax.imshow(data * 1e9, aspect='auto', cmap='viridis',
                       extent=[0, scan_size_um, 0, data.shape[0]], origin='upper')
        self.fig.colorbar(im, ax=ax, label='Height (nm)')
        ax.set_xlabel("X (µm)", fontsize=11)
        ax.set_ylabel("Scan line", fontsize=11)
        set_titles(ax, "Topography", f"{data.shape[0]} scan lines",
                   theming.current())
        self.draw()

    @_remembered
    def plot_row_groups(self, disp_um, profiles_nm, group_info, filename, scan_size_um):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        n_groups = group_info['n_groups']
        cmap = matplotlib.colormaps['coolwarm']

        for i, profile in enumerate(profiles_nm):
            ax.plot(disp_um, profile, color=cmap(i / max(n_groups - 1, 1)),
                    linewidth=0.8, alpha=0.6)

        ax.plot(disp_um, np.mean(np.array(profiles_nm), axis=0), color=theming.current().ink,
                linewidth=1.8, label=f'Mean (n={n_groups})', zorder=5)
        ax.set_xlabel("Position (µm)", fontsize=11)
        ax.set_ylabel("Height (nm)", fontsize=11)
        set_titles(ax, "Row-group profiles",
                   f"{n_groups} groups of {group_info['rows_per_group']} rows, "
                   f"scan width {scan_size_um:.3f} µm", theming.current())
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        sm = matplotlib.cm.ScalarMappable(
            cmap=cmap,
            norm=matplotlib.colors.Normalize(vmin=0, vmax=group_info['n_rows']))
        sm.set_array([])
        self.fig.colorbar(sm, ax=ax, pad=0.02).set_label(
            "Image row (top → bottom)", fontsize=9)
        self.draw()

    @_remembered
    def plot_groove_detection(self, result, filename, settings):
        """Flattened profile with the grooves the analysis used"""
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        x, y = result['raw_x'], result['flat_y']
        centers = result['groove_centers']

        ax.plot(x, y, '-', color=theming.current().ink, linewidth=0.7, label='Flattened profile')
        ax.plot(x[centers], y[centers], 'ro', markersize=9,
                label=f'Detected grooves (N={len(centers)})')

        # Shade the excluded zone at each scan edge, so it is visible why a
        # groove at the boundary was dropped.
        edge = settings.edge_exclusion_periods
        if edge > 0 and len(centers) > 0:
            margin = edge * result['period_nm'] / 1000.0
            for lo, hi in ((x[0], x[0] + margin), (x[-1] - margin, x[-1])):
                ax.axvspan(lo, hi, color='red', alpha=0.10)
            ax.axvspan(x[0], x[0], color='red', alpha=0.10,
                       label=f'Edge exclusion ({edge}× period)')

        ax.set_xlabel("Position (µm)", fontsize=11)
        ax.set_ylabel("Height (nm)", fontsize=11)
        set_titles(ax, "Groove detection", f"{len(centers)} grooves",
                   theming.current())
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        self.draw()

    @_remembered
    def plot_blaze_angles(self, result, filename):
        """Angle distribution, fit quality, and angle against facet width"""
        self.fig.clear()
        ax1, ax2, ax3 = (self.fig.add_subplot(1, 3, i) for i in (1, 2, 3))

        angles = np.asarray(result['all_angles'])
        quality = result['quality']
        mean, sem = result['mean_angle'], result.get('sem', 0.0)

        ax1.hist(angles, bins=20, edgecolor=theming.current().ink, alpha=0.75, color='steelblue')
        ax1.axvline(mean, color='r', linestyle='--', linewidth=2,
                    label=f'Mean {mean:.2f}°')
        ax1.set_xlabel("Blaze angle (°)", fontsize=10)
        ax1.set_ylabel("Count", fontsize=10)
        set_titles(ax1, "Distribution",
                   f"N = {len(angles)}, σ = {result['std_angle']:.2f}°", theming.current())
        ax1.legend(fontsize=9)
        ax1.grid(True, alpha=0.3)

        r2 = np.array([q.get('blaze_r2', np.nan) for q in quality])
        ax2.hist(r2, bins=20, edgecolor=theming.current().ink, alpha=0.75, color='seagreen')
        ax2.set_xlabel("Fit R²", fontsize=10)
        ax2.set_ylabel("Count", fontsize=10)
        set_titles(ax2, "Fit quality", f"worst R² {np.nanmin(r2):.4f}",
                   theming.current())
        ax2.grid(True, alpha=0.3)

        widths = np.array([q.get('blaze_width_nm', np.nan) for q in quality])
        ax3.plot(widths, angles, 'o', markersize=4, alpha=0.45, color='darkorange')
        ax3.set_xlabel("Blaze facet width (nm)", fontsize=10)
        ax3.set_ylabel("Blaze angle (°)", fontsize=10)
        set_titles(ax3, "Angle against facet width", "narrow facets fit poorly",
                   theming.current())
        ax3.grid(True, alpha=0.3)

        self.fig.suptitle(f"Blaze angle {mean:.2f}° ± {sem:.3f}° (SEM)",
                          x=0.01, ha="left", fontsize=12, fontweight="semibold")
        self.draw()
