"""Write what a window is showing as a Python script.

"Copy as Python" in the efficiency and coatings windows puts one of these on
the clipboard: a short script that rebuilds the problem on screen and solves
it the way the window would, for a user who wants to take it from there in
code of their own, such as a sweep the window has no control for.

The script is generated from the *built* objects (the ``Problem``, the
``Illumination``, the wavelengths and the tab's options), not from the form's
text, so it is what the solver was actually given. Nothing here touches Qt;
``tests/test_gui_script.py`` runs the scripts it writes and compares them with
the window's own result.

Called by :mod:`gratinglab.gui.qt.main_window` and
:mod:`gratinglab.gui.qt.coatings_window`.
"""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np
from pydantic import BaseModel

from ..profiles import FromProfileData

__all__ = ["efficiency_script", "mirror_script"]

_HEADER = '"""{title}\n\nWritten by gratinglab {version}, "Copy as Python". Lengths are nm, angles degrees.\n"""\n'


def _version() -> str:
    from .. import __version__

    return __version__


def _value(value: Any) -> str:
    """A value as source: a model as its constructor call, the rest by repr."""
    if isinstance(value, BaseModel):
        return _call(value)
    if isinstance(value, tuple):
        inner = ", ".join(_value(v) for v in value)
        return f"({inner},)" if len(value) == 1 else f"({inner})"
    if isinstance(value, (np.floating, np.integer)):
        return repr(value.item())
    return repr(value)


def _call(model: BaseModel, *, spread: bool = False, **replaced: str) -> str:
    """``Class(field=value, ...)`` for a pydantic model, leaving out fields
    still at their default. ``replaced`` gives source for named fields, and
    ``spread`` puts one argument on a line, as a person would write a long one."""
    parts = []
    for name, field in type(model).model_fields.items():
        value = getattr(model, name)
        if name in replaced:
            parts.append(f"{name}={replaced[name]}")
        elif field.is_required() or value != field.default:
            parts.append(f"{name}={_value(value)}")
    if spread and len(parts) > 2:
        return f"{type(model).__name__}(\n" + "".join(f"    {part},\n" for part in parts) + ")"
    return f"{type(model).__name__}({', '.join(parts)})"


def _wavelengths(values) -> str:
    """The scan as ``np.linspace`` when it is one, else the numbers."""
    values = np.asarray(values, dtype=float)
    if len(values) == 1:
        return f"np.array([{values[0]!r}])"
    even = np.linspace(values[0], values[-1], len(values))
    if np.allclose(values, even, rtol=1e-12, atol=0.0):
        return f"np.linspace({float(values[0])!r}, {float(values[-1])!r}, {len(values)})"
    return f"np.array({[float(v) for v in values]!r})"


def efficiency_script(
    parsed,
    methods: Sequence[str] = (),
    options: Mapping[str, Mapping[str, Any]] | None = None,
    profile_path: str = "",
) -> str:
    """The script for the efficiency window.

    ``parsed`` is the built geometry (``state.Parsed``). With one name in
    ``methods`` the script solves with it; with several it compares them, as
    the Compare tab does; with none it stops at the problem, which is what
    the Setup and Grating Geometry tabs show. ``options`` holds each method's
    settings by name. ``profile_path`` is the file a "From file" profile was
    read from: a measured groove is loaded again rather than written out.
    """
    options = options or {}
    problem, light = parsed.problem, parsed.illumination
    names = {"Illumination", "Problem"}
    lines: list[str] = []

    if isinstance(problem.profile, FromProfileData):
        lines += [
            "from gratinglab.io.boundary import read_profile_file",
            "",
            f"profile = read_profile_file(r{profile_path!r})",
        ]
        problem_source = _call(problem, spread=True, profile="profile")
    else:
        names.add(type(problem.profile).__name__)
        problem_source = _call(problem, spread=True)
        lines.append("")
    if problem.layers:
        names.add("ConformalLayer")

    body = [
        f"problem = {problem_source}",
        f"illumination = {_call(light)}",
        f"wavelengths = {_wavelengths(parsed.wavelengths)}",
        "",
    ]
    if len(methods) == 1:
        (method,) = methods
        names.add("solve")
        settings = options.get(method, {})
        if settings:
            call = (f"scan = solve(\n    problem, illumination, wavelengths, method={method!r},\n"
                    + "".join(f"    {key}={_value(value)},\n" for key, value in settings.items()) + ")")
        else:
            call = f"scan = solve(problem, illumination, wavelengths, method={method!r})"
        body += [
            call,
            "",
            "# One row per wavelength, one column per order.",
            "for m in scan.orders:",
            "    print(f\"order {m:+d}: peak {scan.order(int(m)).max():.4f}\")",
        ]
        title = f"A grating solved with the {method} solver."
    elif methods:
        names.update({"sweep", "align"})
        chosen = {m: dict(options[m]) for m in methods if options.get(m)}
        body += [
            f"methods = {list(methods)!r}",
            f"options = {_value_dict(chosen)}",
            "comparison = align(sweep(problem, illumination, wavelengths, methods, options))",
            "",
            "reference = methods[0]",
            "for other in methods[1:]:",
            "    print(other, \"against\", reference, comparison.summary(other, reference))",
        ]
        title = f"A grating solved with {', '.join(methods)}, and compared."
    else:
        body += ["print(problem)", "print(illumination)"]
        title = "A grating and its illumination."

    imports = ["import numpy as np", "", f"from gratinglab import {', '.join(sorted(names))}"]
    return "\n".join([_HEADER.format(title=title, version=_version()), *imports, *lines, *body]) + "\n"


def _value_dict(mapping: Mapping[str, Mapping[str, Any]]) -> str:
    inner = ", ".join(
        f"{name!r}: {{{', '.join(f'{k!r}: {_value(v)}' for k, v in opts.items())}}}"
        for name, opts in mapping.items()
    )
    return "{" + inner + "}"


def mirror_script(result) -> str:
    """The script for the coatings window's Mirror tab, from its last result
    (``mirror_view.MirrorResult``): the stack, the scan, and the three
    things the tab reports."""
    mirror = result.mirror
    names = {"Mirror"}
    if mirror.layers:
        names.add("ConformalLayer")
    args = [f"substrate={mirror.substrate!r}"]
    if mirror.layers:
        args.append(f"layers={_value(tuple(mirror.layers))}")
    if mirror.roughness:
        args.append(f"roughness={mirror.roughness!r}")

    scan = _wavelengths(result.x)
    if result.scan == "angle":
        axes = [f"graze = {scan}  # degrees from the surface", f"wavelength = {result.fixed!r}  # nm"]
    else:
        axes = [f"graze = {result.fixed!r}  # degrees from the surface", f"wavelength = {scan}  # nm"]
    body = [
        f"mirror = Mirror({', '.join(args)})",
        *axes,
        "",
        f"reflectivity = mirror.reflectivity(graze, wavelength, polarization={result.polarization!r})",
        "depth = mirror.penetration_depth(graze, wavelength)  # 1/e intensity depth, nm, by material",
        "",
        "print(f\"reflectivity from {reflectivity.min():.4g} to {reflectivity.max():.4g}\")",
        "for material, d in depth.items():",
        "    print(f\"1/e depth in {material}: {d.min():.3g} to {d.max():.3g} nm\")",
    ]
    if mirror.layers:
        one = "graze[0], wavelength" if result.scan == "angle" else "graze, wavelength[0]"
        body += [
            "",
            "# The thinnest top layer that reads as bulk, at the first point of the scan.",
            f"needed = mirror.bulk_thickness({one}, tolerance={result.tolerance!r}, "
            f"polarization={result.polarization!r})",
            f"print(f\"{mirror.layers[-1].material} needed to read as bulk: {{needed:.3g}} nm\")",
        ]
    title = "A coating stack as a flat mirror."
    imports = ["import numpy as np", "", f"from gratinglab import {', '.join(sorted(names))}", ""]
    return "\n".join([_HEADER.format(title=title, version=_version()), *imports, *body]) + "\n"
