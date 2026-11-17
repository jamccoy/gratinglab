"""What the Mirror tab computes, draws and says, with no Qt in it.

The stack comes from the geometry dock's Surface group (its Grating material
as the substrate, the layers and the roughness), so a coating is described in
one place whether it is about to be ruled or only reflected from. The tab's
own settings choose the scan: graze angle at one wavelength, or wavelength at
one graze, over the dock's wavelengths or the tab's own. :func:`compute`
returns a :class:`MirrorResult`; :func:`verdict` is the one-line answer,
:func:`draw` and :func:`lines` the two plots and the panel text. An angle
scan marks the graze angles named in the tab's own Mark angles field, so the
angles a grating will actually be used at can be read off the curve. The
efficiency window fills that field from its mount when it sends a coating
over (:func:`mount_markers`, :func:`marks_text`).

Called by :mod:`gratinglab.gui.qt.mirror_tab`; tested in
``tests/test_mirror_view.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from ..geometry import facet_graze
from ..mirror import Mirror
from ..profiles import Blazed
from .provenance import Line
from .state import FieldError, FormErrors, FormState, _build_layers, _number, build
from .theme import order_colors, set_titles, style_axes

__all__ = [
    "MirrorSettings", "MirrorResult", "compute", "draw", "lines", "verdict",
    "is_bulk_throughout", "mount_markers", "marks_text", "parse_marks",
]

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
    #: A wavelength scan uses the dock's wavelengths, or the three below.
    dock_wavelengths: bool = True
    wavelength_start: str = "1.0"
    wavelength_stop: str = "5.0"
    wavelength_count: str = "200"
    #: Graze angles to mark on an angle scan: ``"1.5"`` or
    #: ``"grating plane 1.36, blaze facet 1.5"``.
    marks: str = ""


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
    #: The dock's mount as ``(label, graze in deg)``, for an angle scan.
    markers: tuple[tuple[str, float], ...] = ()


def mount_markers(form: FormState) -> tuple[tuple[str, float], ...]:
    """The graze angles the dock's mount puts light at: onto the grating
    plane, and onto the blaze facet when the profile has one. Empty when the
    form does not describe a grating, which is no concern of a mirror's."""
    try:
        parsed = build(form)
    except FormErrors:
        return ()
    light, profile = parsed.illumination, parsed.problem.profile
    out = [("grating plane", float(np.degrees(facet_graze(light.gamma, 0.0, light.alpha))))]
    if isinstance(profile, Blazed):
        blaze = np.radians(profile.blaze_angle)
        out.append(("blaze facet", float(np.degrees(facet_graze(light.gamma, blaze, light.alpha)))))
    return tuple(out)


def marks_text(markers) -> str:
    """Markers as the Mark angles field holds them."""
    return ", ".join(f"{label} {graze:.4g}" for label, graze in markers)


def parse_marks(text: str, errors: list[FieldError]) -> tuple[tuple[str, float], ...]:
    """``"grating plane 1.36, 2"`` as ``(("grating plane", 1.36), ("marked", 2.0))``.

    Each entry is one angle in degrees, with a label before or after it if
    you want one. A bad entry is an error on the ``marks`` field.
    """
    out = []
    for entry in (part.strip() for part in text.split(",")):
        if not entry:
            continue
        words = entry.split()
        numbers = []
        for word in words:
            try:
                numbers.append((word, float(word)))
            except ValueError:
                pass
        if len(numbers) != 1 or not 0 < numbers[0][1] <= 90:
            errors.append(FieldError("marks", f"{entry!r} should be one graze angle in "
                                              "degrees, with a label if you like"))
            continue
        label = " ".join(w for w in words if w != numbers[0][0]) or "marked"
        out.append((label, numbers[0][1]))
    return tuple(out)


def compute(form: FormState, settings: MirrorSettings) -> MirrorResult:
    """Build the mirror and evaluate it, raising :class:`FormErrors` for
    every bad field at once."""
    errors: list[FieldError] = []
    substrate = form.coating.strip()
    if not substrate:
        errors.append(FieldError("coating", "a mirror needs a substrate"))
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
        source = form if settings.dock_wavelengths else settings
        start = _number(source.wavelength_start, "wavelength_start", errors, minimum=0.0, exclusive_min=True)
        stop = _number(source.wavelength_stop, "wavelength_stop", errors, minimum=0.0, exclusive_min=True)
        count = _number(source.wavelength_count, "wavelength_count", errors, minimum=2, maximum=MAX_POINTS, integer=True)
        fixed = _number(settings.graze, "graze", errors, minimum=0.0, exclusive_min=True, maximum=90.0)
    else:
        raise ValueError(f"scan must be 'angle' or 'wavelength', got {settings.scan!r}")
    markers = parse_marks(settings.marks, errors) if settings.scan == "angle" else ()
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
    return MirrorResult(mirror, settings.scan, x, fixed, reflectivity, depths, bulk, tolerance,
                        settings.polarization, markers)


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
    for i, (label, graze) in enumerate(result.markers):
        if result.x[0] <= graze <= result.x[-1]:
            for axes in (r_axes, depth_axes):
                axes.axvline(graze, color=theme.muted, linewidth=1)
            # At the foot, clear of the curve near grazing, and on alternate
            # sides of their lines: the two are often a fraction of a degree apart.
            left = i % 2 == 0
            r_axes.annotate(f"{label} {graze:.3g}°", (graze, 0.0), xycoords=("data", "axes fraction"),
                            xytext=(-3 if left else 3, 4), textcoords="offset points",
                            ha="right" if left else "left", va="bottom", rotation=90,
                            fontsize=8, color=theme.muted)
    depth_axes.set_yscale("log")
    depth_axes.set_ylabel("depth (nm)")
    depth_axes.set_xlabel(_x_label(result))
    set_titles(depth_axes, "How deep", f"bulk within {result.tolerance:g} of R", theme)
    depth_axes.legend(fontsize=8, loc="best")


def is_bulk_throughout(result: MirrorResult) -> bool:
    """The outermost layer reads as bulk at every point of the scan, so a
    grating calculation may take it as the grating material alone."""
    return result.bulk is not None and bool((result.bulk <= result.mirror.layers[-1].thickness).all())


def verdict(result: MirrorResult) -> Line:
    """The answer in one line: is the top layer thick enough to be bulk?"""
    m = result.mirror
    if result.bulk is None:
        return Line(f"Bare {m.substrate}: there is no layer to test.", "dim")
    top = m.layers[-1]
    enough = result.bulk <= top.thickness
    what = f"{top.material} at {top.thickness:g} nm"
    if enough.all():
        return Line(f"{what} reads as bulk across the scan.", "ok")
    if not enough.any():
        return Line(f"{what} is not bulk anywhere in the scan.", "warn")
    return Line(f"{what} is bulk at {enough.sum()} of {enough.size} points.", "warn")


def lines(result: MirrorResult) -> tuple[Line, ...]:
    """The panel: the stack, then what the verdict means for modelling."""
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
    out.append(verdict(result))
    if enough.all():
        out.append(Line(f" What is beneath it changes R by under {result.tolerance:g}, so a grating "
                        f"calculation may take {top.material} as the grating material alone.\n"))
    elif not enough.any():
        finite = result.bulk[np.isfinite(result.bulk)]
        need = f" It would need {finite.min():.3g} nm at best." if finite.size else ""
        out.append(Line(f"{need} Model the stack, not {top.material} alone.\n"))
    else:
        x = result.x[enough]
        where = f"at {x[0]:.4g} {unit}" if x.size == 1 else f"between {x.min():.4g} and {x.max():.4g} {unit}"
        out.append(Line(f" That is {where}. Elsewhere the layers beneath change R by more than "
                        f"{result.tolerance:g}, so model the stack there.\n"))
    depth = result.depths[top.material]
    out.append(Line(f"1/e depth in {top.material}: {np.nanmin(depth):.3g} to {np.nanmax(depth):.3g} nm "
                    "over the scan.\n", "dim"))
    for label, graze in result.markers:
        out.append(Line(f"Marked: {label} at {graze:.3g} deg graze.\n", "dim"))
    return tuple(out)
