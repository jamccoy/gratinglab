"""What the Mirror tab computes, draws and says, with no Qt in it.

The stack comes from the geometry dock (Coating as the substrate, Layers,
Roughness, and the wavelength scan), so a coating is described in one place
whether it is about to be ruled or only reflected from. The tab's own
settings choose the scan: graze angle at one wavelength, or wavelength at one
graze. :func:`compute` returns a :class:`MirrorResult`; :func:`draw` and
:func:`lines` turn it into the two plots and the panel text.

Called by :mod:`gratinglab.gui.qt.mirror_tab`; tested in
``tests/test_mirror_view.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ..mirror import Mirror
from .provenance import Line
from .state import FieldError, FormErrors, FormState, _build_layers, _number
from .theme import order_colors, set_titles, style_axes

__all__ = ["MirrorSettings", "MirrorResult", "compute", "draw", "lines"]

#: Most points one scan evaluates; each one sweeps a thickness for the bulk test.
MAX_POINTS = 1000


@dataclass(frozen=True, slots=True)
class MirrorSettings:
    """The tab's own fields, as typed."""

    scan: str = "angle"  # "angle" or "wavelength"
    graze_start: str = "0.1"
    graze_stop: str = "6.0"
    graze_count: str = "300"
    graze: str = "2.0"  # deg, the fixed graze for a wavelength scan
    wavelength: str = "1.0"  # nm, the fixed wavelength for an angle scan
    polarization: str = "unpolarized"
    tolerance: str = "1e-3"


@dataclass(frozen=True)
class MirrorResult:
    mirror: Mirror
    scan: str
    x: NDArray[np.float64]  # graze in deg, or wavelength in nm
    fixed: float  # the wavelength (angle scan) or graze (wavelength scan)
    reflectivity: NDArray[np.float64]
    depths: dict[str, NDArray[np.float64]]
    #: Thickness the outermost layer needs to read as bulk, per point, or None.
    bulk: NDArray[np.float64] | None
    tolerance: float
    polarization: str


def compute(form: FormState, settings: MirrorSettings) -> MirrorResult:
    """Build the mirror and evaluate it, raising :class:`FormErrors` for
    every bad field at once."""
    errors: list[FieldError] = []
    substrate = form.coating.strip()
    if not substrate:
        errors.append(FieldError("coating", "a mirror needs a substrate: set Coating"))
    roughness = _number(form.roughness, "roughness", errors, minimum=0.0)
    layers = _build_layers(form.layers, errors)
    tolerance = _number(settings.tolerance, "tolerance", errors, minimum=0.0, exclusive_min=True, maximum=1.0)
    if settings.polarization not in ("s", "p", "unpolarized"):
        errors.append(FieldError("polarization", "must be s, p or unpolarized"))
    if settings.scan == "angle":
        start = _number(settings.graze_start, "graze_start", errors, minimum=0.0, exclusive_min=True, maximum=90.0)
        stop = _number(settings.graze_stop, "graze_stop", errors, minimum=0.0, exclusive_min=True, maximum=90.0)
        count = _number(settings.graze_count, "graze_count", errors, minimum=2, maximum=MAX_POINTS, integer=True)
        fixed = _number(settings.wavelength, "wavelength", errors, minimum=0.0, exclusive_min=True)
    elif settings.scan == "wavelength":
        start = _number(form.wavelength_start, "wavelength_start", errors, minimum=0.0, exclusive_min=True)
        stop = _number(form.wavelength_stop, "wavelength_stop", errors, minimum=0.0, exclusive_min=True)
        count = _number(form.wavelength_count, "wavelength_count", errors, minimum=2, maximum=MAX_POINTS, integer=True)
        fixed = _number(settings.graze, "graze", errors, minimum=0.0, exclusive_min=True, maximum=90.0)
    else:
        raise ValueError(f"scan must be 'angle' or 'wavelength', got {settings.scan!r}")
    if start is not None and stop is not None and stop <= start:
        errors.append(FieldError(f"{'graze' if settings.scan == 'angle' else 'wavelength'}_stop", "must exceed the start"))
    if errors:
        raise FormErrors(tuple(errors))

    mirror = Mirror(substrate, layers, roughness)
    x = np.linspace(start, stop, int(count))
    graze, wl = (x, np.full_like(x, fixed)) if settings.scan == "angle" else (np.full_like(x, fixed), x)
    reflectivity = mirror.reflectivity(graze, wl, settings.polarization)
    depths = mirror.penetration_depth(graze, wl)
    bulk = None
    if layers:
        bulk = np.array([
            mirror.bulk_thickness(g, w, tolerance, settings.polarization) for g, w in zip(graze, wl)
        ])
    return MirrorResult(mirror, settings.scan, x, fixed, reflectivity, depths, bulk, tolerance, settings.polarization)


def _x_label(result: MirrorResult) -> str:
    return "graze angle (deg)" if result.scan == "angle" else "wavelength (nm)"


def _fixed(result: MirrorResult) -> str:
    return f"at {result.fixed:g} nm" if result.scan == "angle" else f"at {result.fixed:g} deg graze"


def draw(r_axes, depth_axes, result: MirrorResult, theme, log_r: bool = True) -> None:
    """Reflectivity above; 1/e depth per material and the bulk thickness below."""
    for axes in (r_axes, depth_axes):
        axes.clear()
        style_axes(axes, theme)
    r_axes.plot(result.x, result.reflectivity, color=theme.ink)
    r_axes.set_yscale("log" if log_r else "linear")
    r_axes.set_ylabel("reflectivity")
    set_titles(r_axes, "Specular reflectivity", f"{result.polarization}, {_fixed(result)}", theme)

    colors = order_colors(range(len(result.depths)), theme)
    for i, (name, depth) in enumerate(result.depths.items()):
        depth_axes.plot(result.x, depth, color=colors[i], label=f"1/e depth in {name}")
    if result.bulk is not None:
        top = result.mirror.layers[-1]
        depth_axes.plot(result.x, result.bulk, color=theme.ink, linestyle="--",
                        label=f"{top.material} needed to read as bulk")
        depth_axes.axhline(top.thickness, color=theme.muted, linestyle=":",
                           label=f"{top.material} as deposited, {top.thickness:g} nm")
    depth_axes.set_yscale("log")
    depth_axes.set_ylabel("depth (nm)")
    depth_axes.set_xlabel(_x_label(result))
    set_titles(depth_axes, "How deep", f"bulk within {result.tolerance:g} of R", theme)
    depth_axes.legend(fontsize=8, loc="best")


def lines(result: MirrorResult) -> tuple[Line, ...]:
    """The panel: the stack, then the answer to "is the top layer bulk"."""
    m = result.mirror
    stack = " / ".join(
        ["vacuum", *(f"{layer.material} {layer.thickness:g} nm" for layer in reversed(m.layers)),
         f"{m.substrate} (substrate)"]
    )
    out = [Line("stack: ", "dim"), Line(stack + "\n"), Line(f"{_fixed(result)}, {result.polarization}", "dim")]
    if m.roughness:
        out.append(Line(f", {m.roughness:g} nm RMS at every interface", "dim"))
    out.append(Line("\n"))
    if result.bulk is None:
        out.append(Line("No layers: the substrate is the whole mirror, so there is no thickness to test.\n", "dim"))
        return tuple(out)
    top = m.layers[-1]
    enough = result.bulk <= top.thickness
    unit = "deg" if result.scan == "angle" else "nm"
    if enough.all():
        out += [Line(f"{top.material} at {top.thickness:g} nm reads as bulk over the whole scan", "ok"),
                Line(f": what is beneath it changes R by under {result.tolerance:g}, so a grating "
                     f"calculation may take {top.material} as the coating alone.\n")]
    elif not enough.any():
        finite = result.bulk[np.isfinite(result.bulk)]
        need = f"; it would need {finite.min():.3g} nm at best" if finite.size else ""
        out += [Line(f"{top.material} at {top.thickness:g} nm is not bulk anywhere in the scan", "warn"),
                Line(f"{need}. Model the stack, not {top.material} alone.\n")]
    else:
        x = result.x[enough]
        where = f"at {x[0]:.4g} {unit}" if x.size == 1 else f"between {x.min():.4g} and {x.max():.4g} {unit}"
        out += [Line(f"{top.material} at {top.thickness:g} nm reads as bulk at {enough.sum()} of "
                     f"{enough.size} points", "warn"),
                Line(f", {where}. Elsewhere the layers beneath change R by more than "
                     f"{result.tolerance:g}, so model the stack there.\n")]
    depth = result.depths[top.material]
    out.append(Line(f"1/e depth in {top.material}: {np.nanmin(depth):.3g} to {np.nanmax(depth):.3g} nm "
                    "over the scan.\n", "dim"))
    return tuple(out)
