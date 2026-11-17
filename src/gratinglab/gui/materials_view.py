"""What the Materials tab resolves, draws and says, with no Qt in it.

The tab answers three things about one material before it is used: what name
to write for it, what its index looks like, and whether it covers the
wavelengths the dock is about to scan. :func:`material_name` builds the name
from the tab's fields, so nobody has to remember ``factors:SiO2@2.2``;
:func:`preview` resolves it and checks the scan; :func:`draw` and
:func:`lines` turn that into the plot and the panel text.

Called by :mod:`gratinglab.gui.qt.materials_tab`; tested in
``tests/test_materials_view.py``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..materials import OpticalConstants, lookup
from ..materials.optical import HC_EV_NM
from .provenance import Line
from .theme import order_colors, set_titles, style_axes

__all__ = ["KINDS", "MaterialPreview", "material_name", "preview", "draw", "lines"]

#: The ways to name a material, as (key, label, what the name field holds).
KINDS = (
    ("table", "Table in your library", "Name"),
    ("factors", "Scattering factors + density", "Formula"),
    ("typical", "Typical-material model", "Formula"),
)


def material_name(kind: str, name: str, density: str = "") -> str:
    """The coating name for the tab's fields; ``""`` when nothing is entered.

    A density is optional everywhere: a table is rescaled to it
    (``Au@17.5``), and a formula with none uses its default where one exists.
    """
    name, density = name.strip(), density.strip()
    if not name:
        return ""
    prefix = {"table": "", "factors": "factors:", "typical": "typical:"}[kind]
    return f"{prefix}{name}" + (f"@{density}" if density else "")


@dataclass(frozen=True)
class MaterialPreview:
    name: str
    constants: OpticalConstants
    #: The dock's scan, nm, or None when its fields do not parse.
    scan: tuple[float, float] | None
    #: X-ray-like over most of its range (index near 1): plot the decrement
    #: and absorption on a log axis rather than n and k.
    near_unity: bool


def _scan(start: str, stop: str) -> tuple[float, float] | None:
    try:
        low, high = float(start), float(stop)
    except ValueError:
        return None
    return (low, high) if 0 < low < high else None


def preview(name: str, scan_start: str = "", scan_stop: str = "") -> MaterialPreview:
    """Resolve ``name`` and note the scan. Raises what :func:`lookup` raises."""
    constants = lookup(name)
    # Most of the table, not all of it: an X-ray table's longest wavelengths
    # reach |n - 1| of 0.2, and one such corner should not flatten the plot.
    near = bool(constants.near_unity(constants.wavelength_nm).mean() > 0.5)
    return MaterialPreview(name, constants, _scan(scan_start, scan_stop), near)


def covers_scan(p: MaterialPreview) -> bool | None:
    """Whether the table covers the whole scan; None with no scan to check."""
    if p.scan is None:
        return None
    return bool(p.constants.covers(np.array(p.scan)).all())


def draw(axes, p: MaterialPreview, theme) -> None:
    """The index against wavelength, with the dock's scan shaded."""
    axes.clear()
    style_axes(axes, theme)
    c = p.constants
    colors = order_colors(range(2), theme)
    if p.near_unity:
        # Log axes drop a non-positive decrement, as near an edge.
        axes.plot(c.wavelength_nm, np.where(c.decrement > 0, c.decrement, np.nan),
                  color=colors[0], label="decrement δ")
        axes.plot(c.wavelength_nm, np.where(c.absorption > 0, c.absorption, np.nan),
                  color=colors[1], label="absorption β")
        axes.set_yscale("log")
        axes.set_ylabel("1 − n,  k")
    else:
        axes.plot(c.wavelength_nm, 1.0 - c.decrement, color=colors[0], label="n")
        axes.plot(c.wavelength_nm, c.absorption, color=colors[1], label="k")
        axes.set_yscale("linear")
        axes.set_ylabel("n,  k")
    axes.set_xscale("log")
    axes.set_xlabel("wavelength (nm)")
    if p.scan is not None:
        axes.axvspan(*p.scan, color=theme.muted, alpha=0.18, label="the dock's scan")
    low, high = c.range_nm
    set_titles(axes, p.name, f"{low:.4g} to {high:.4g} nm, {len(c.wavelength_nm)} points", theme)
    axes.legend(fontsize=8, loc="best")


def _short(source: str) -> str:
    """``source`` with the home directory written ``~``: the library a table
    came from still shows, without a user's whole path in front of it."""
    return source.replace(str(Path.home()), "~")


def lines(p: MaterialPreview) -> tuple[Line, ...]:
    """Where it came from, what it covers, and whether that includes the scan."""
    c = p.constants
    low, high = c.range_nm
    out = [
        Line("name: ", "dim"), Line(f"{p.name}\n"),
        Line("source: ", "dim"), Line(f"{_short(c.source)}\n"),
        Line("covers: ", "dim"),
        Line(f"{low:.4g} to {high:.4g} nm ({HC_EV_NM / high:.4g} to {HC_EV_NM / low:.4g} eV)\n"),
    ]
    covered = covers_scan(p)
    if covered is None:
        out.append(Line("The dock's wavelength scan is not set, so coverage is not checked.\n", "dim"))
    elif covered:
        out.append(Line(f"Covers the dock's scan, {p.scan[0]:g} to {p.scan[1]:g} nm.\n", "ok"))
    else:
        out.append(Line(
            f"Does not cover the dock's scan, {p.scan[0]:g} to {p.scan[1]:g} nm: a solve "
            "with it would be refused outside the table.\n", "warn"))
    if covered:
        grid = np.linspace(p.scan[0], p.scan[1], 200)
        for caveat in c.caveats(grid):
            out.append(Line(f"{caveat}\n", "warn"))
    return tuple(out)
