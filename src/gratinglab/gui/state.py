"""Turn form fields into a validated problem.

This module contains the GUI's only substantial logic, and it is kept
separate from the widgets so that it can be tested headlessly. Nothing here
imports a GUI toolkit, and nothing here computes physics: the module parses
strings, validates ranges and passes the result to the tested core.

The separation is important because widget code is the least testable part of
any project, while the value of this project lies in defensible numbers.
Keeping parsing and validation here means that a bug in the untested layer
can produce a wrong layout but not a wrong number.

Errors are collected instead of being raised at the first failure, so a form
with three bad fields reports three problems and the user need not fix them
one at a time.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from endstation.forms import FieldError as _FieldError
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ..illumination import Illumination
from ..problem import ConformalLayer, Problem
from ..profiles import Blazed, Lamellar, Sinusoidal

__all__ = [
    "FieldError",
    "FormErrors",
    "FormState",
    "Parsed",
    "PROFILE_KINDS",
    "PROFILE_FIELDS",
    "MOUNTS",
    "ANGLE_LABELS",
    "build",
    "validate",
]

#: Selectable profile kinds, in the order they should appear in the UI.
PROFILE_KINDS = ("Blazed", "Lamellar", "Sinusoidal", "From file")

#: Which form fields each profile kind uses. The UI hides the rest,
#: because a duty cycle means nothing for a sinusoid and offering it invites
#: the reasonable assumption that it does something.
#:
#: Here rather than in the window for the same reason as ANGLE_LABELS: what a
#: control *means* is form logic, and "the dropdown offers a kind nothing maps"
#: should be a test rather than an empty panel.
PROFILE_FIELDS: dict[str, frozenset[str]] = {
    "Blazed": frozenset({"blaze_angle", "antiblaze_angle"}),
    "Lamellar": frozenset({"depth_fraction", "duty_cycle"}),
    "Sinusoidal": frozenset({"depth_fraction"}),
    "From file": frozenset({"profile_path"}),
}

#: Selectable mounts.
MOUNTS = ("Classical", "Conical", "Off-plane")

#: What the two angle fields *mean* in each mount. The form always carries two
#: angles; only their interpretation changes, so the UI relabels rather than
#: swapping widgets.
ANGLE_LABELS: dict[str, tuple[str, str | None]] = {
    "Classical": ("α incidence (deg)", None),
    "Conical": ("θ polar (deg)", "φ azimuth (deg)"),
    "Off-plane": ("α azimuth (deg)", "γ graze (deg)"),
}


#: One thing wrong with one field. Shared with `metrology.gui.state`, which
#: had its own identical copy, and with the other projects here -- see
#: :mod:`endstation.forms`. Re-exported so every call site keeps importing it
#: from the module it is used in.
FieldError = _FieldError


class FormErrors(ValueError):
    """Raised by :func:`build` when a form cannot be turned into a problem."""

    def __init__(self, errors: tuple[FieldError, ...]) -> None:
        self.errors = errors
        super().__init__("; ".join(str(e) for e in errors))


@dataclass(frozen=True, slots=True)
class FormState:
    """The grating geometry, exactly as typed, held as strings as in a form.

    It holds the geometry only (period, profile, mount, wavelengths), which is
    shared by every solver a tab might run. A solver's own knobs, such as
    scalar's ``quadrature_points``, live in that solver's own state, e.g.
    :class:`gratinglab.gui.scalar_options.ScalarOptionsState`, because this
    dataclass cannot know which solver is about to use it.

    The defaults describe a soft-X-ray off-plane case, which is the project's
    primary application and exercises the relevant physics immediately.
    """

    period: str = "315.15"

    profile_kind: str = "Blazed"
    blaze_angle: str = "29.5"
    antiblaze_angle: str = "70.5"
    depth_fraction: str = "0.30"
    duty_cycle: str = "0.50"
    profile_path: str = ""

    mount: str = "Off-plane"
    alpha: str = "25.0"
    gamma: str = "1.5"

    wavelength_start: str = "1.0"
    wavelength_stop: str = "5.0"
    wavelength_count: str = "200"

    #: Material name, or empty for none. A `Problem` field like any other --
    #: it belongs beside the period rather than on a tab, because one
    #: `Problem` feeds every solver and two places to set it would be two
    #: sources of truth. Empty means relative efficiency, which is a correct
    #: answer to a different question rather than a missing input.
    coating: str = ""
    #: RMS surface roughness in nm. Damps the reflectivity, so it does nothing
    #: without a coating -- and the solver does not pretend otherwise.
    roughness: str = "0.0"
    #: Conformal layers on the grating, innermost first, written
    #: ``"Au 20, Cr 5"``: a material name and a thickness in nm per layer
    #: (:class:`~gratinglab.problem.ConformalLayer`). Empty means a bare grating.
    #: A `Problem` field, so it lives beside the coating; only a solver that
    #: declares `handles_layers` will take a problem that has any.
    layers: str = ""

    def with_field(self, name: str, value: Any) -> "FormState":
        """Return a copy with one field changed."""
        return replace(self, **{name: value})


@dataclass(frozen=True, slots=True)
class Parsed:
    """A geometry that made sense, ready to hand to any solver.

    There is no ``options`` field. The geometry is shared by every solver a
    tab might run, but each solver's options, and their validation (e.g.
    scalar's Nyquist guard in :mod:`gratinglab.gui.scalar_options`), belong to
    that solver and are built separately by whichever tab is about to solve.
    """

    problem: Problem
    illumination: Illumination
    wavelengths: NDArray[np.float64]

    @property
    def lambda_over_period(self) -> float:
        """Ratio at the long-wavelength end, used as the scalar validity indicator."""
        return float(self.wavelengths.max() / self.problem.period)


def _number(
    raw: str,
    field: str,
    errors: list[FieldError],
    *,
    minimum: float | None = None,
    maximum: float | None = None,
    integer: bool = False,
    exclusive_min: bool = False,
    exclusive_max: bool = False,
) -> float | None:
    """Parse one numeric field, appending an error rather than raising."""
    text = raw.strip()
    if not text:
        errors.append(FieldError(field, "required"))
        return None
    try:
        value = float(text)
    except ValueError:
        errors.append(FieldError(field, f"{raw!r} is not a number"))
        return None
    if not np.isfinite(value):
        errors.append(FieldError(field, "must be finite"))
        return None
    if integer and value != int(value):
        errors.append(FieldError(field, "must be a whole number"))
        return None
    if minimum is not None:
        too_small = value <= minimum if exclusive_min else value < minimum
        if too_small:
            relation = "greater than" if exclusive_min else "at least"
            errors.append(FieldError(field, f"must be {relation} {minimum:g}"))
            return None
    if maximum is not None:
        too_big = value >= maximum if exclusive_max else value > maximum
        if too_big:
            relation = "less than" if exclusive_max else "at most"
            errors.append(FieldError(field, f"must be {relation} {maximum:g}"))
            return None
    return value


def _build_profile(form: FormState, errors: list[FieldError]):
    """Construct the profile named by ``profile_kind``."""
    kind = form.profile_kind
    if kind not in PROFILE_KINDS:
        errors.append(
            FieldError("profile_kind", f"unknown; choose one of {', '.join(PROFILE_KINDS)}")
        )
        return None

    if kind == "Blazed":
        blaze = _number(
            form.blaze_angle, "blaze_angle", errors,
            minimum=0.0, maximum=90.0, exclusive_min=True, exclusive_max=True,
        )
        antiblaze = _number(
            form.antiblaze_angle, "antiblaze_angle", errors,
            minimum=0.0, maximum=90.0, exclusive_min=True,
        )
        if blaze is None or antiblaze is None:
            return None
        return Blazed(blaze_angle=blaze, antiblaze_angle=antiblaze)

    if kind == "Lamellar":
        depth = _number(
            form.depth_fraction, "depth_fraction", errors, minimum=0.0, exclusive_min=True
        )
        duty = _number(
            form.duty_cycle, "duty_cycle", errors,
            minimum=0.0, maximum=1.0, exclusive_min=True, exclusive_max=True,
        )
        if depth is None or duty is None:
            return None
        return Lamellar(depth_fraction=depth, duty_cycle=duty)

    if kind == "Sinusoidal":
        depth = _number(
            form.depth_fraction, "depth_fraction", errors, minimum=0.0, exclusive_min=True
        )
        return None if depth is None else Sinusoidal(depth_fraction=depth)

    # From file
    if not form.profile_path.strip():
        errors.append(FieldError("profile_path", "choose a profile file"))
        return None
    path = Path(form.profile_path).expanduser()
    if not path.is_file():
        errors.append(FieldError("profile_path", f"no such file: {path}"))
        return None
    try:
        # A gratinglab .json or a PCGrate .ggp, by suffix. The JSON's period is
        # applied to the form when the file is loaded (GeometryPanel), not here:
        # the period field stays the one place the period is read from.
        from ..io.boundary import read_profile_file

        return read_profile_file(path)
    except Exception as exc:  # the reader raises ValueError with a useful message
        errors.append(FieldError("profile_path", str(exc).split("\n")[0]))
        return None


def _build_illumination(form: FormState, errors: list[FieldError]) -> Illumination | None:
    """Construct the illumination for the selected mount.

    The form always carries two angles; the mount decides what they mean.
    See :data:`ANGLE_LABELS`, which the UI uses to relabel the same widgets.
    """
    mount = form.mount
    if mount not in MOUNTS:
        errors.append(FieldError("mount", f"unknown; choose one of {', '.join(MOUNTS)}"))
        return None

    if mount == "Classical":
        alpha = _number(
            form.alpha, "alpha", errors,
            minimum=-90.0, maximum=90.0, exclusive_min=True, exclusive_max=True,
        )
        return None if alpha is None else Illumination.classical(
            alpha=alpha, polarization="unpolarized"
        )

    if mount == "Conical":
        theta = _number(
            form.alpha, "alpha", errors,
            minimum=0.0, maximum=90.0, exclusive_max=True,
        )
        phi = _number(form.gamma, "gamma", errors, minimum=-180.0, maximum=180.0)
        if theta is None or phi is None:
            return None
        try:
            return Illumination.conical(theta=theta, phi=phi, polarization="unpolarized")
        except Exception as exc:
            errors.append(FieldError("mount", str(exc).split("\n")[0]))
            return None

    graze = _number(
        form.gamma, "gamma", errors,
        minimum=0.0, maximum=90.0, exclusive_min=True,
    )
    azimuth = _number(
        form.alpha, "alpha", errors,
        minimum=-90.0, maximum=90.0, exclusive_min=True, exclusive_max=True,
    )
    if graze is None or azimuth is None:
        return None
    return Illumination.offplane(
        graze=graze, azimuth=azimuth, polarization="unpolarized"
    )


def build(form: FormState) -> Parsed:
    """Validate a form and construct the objects a solver needs.

    Raises
    ------
    FormErrors
        Carrying every problem found, so the UI can mark all bad fields at
        once instead of revealing them one refresh at a time.
    """
    errors: list[FieldError] = []

    period = _number(form.period, "period", errors, minimum=0.0, exclusive_min=True)
    profile = _build_profile(form, errors)
    illumination = _build_illumination(form, errors)

    start = _number(
        form.wavelength_start, "wavelength_start", errors, minimum=0.0, exclusive_min=True
    )
    stop = _number(
        form.wavelength_stop, "wavelength_stop", errors, minimum=0.0, exclusive_min=True
    )
    count = _number(
        form.wavelength_count, "wavelength_count", errors, minimum=1, integer=True
    )
    if start is not None and stop is not None and stop <= start:
        errors.append(
            FieldError("wavelength_stop", f"must exceed the start ({start:g} nm)")
        )

    # With the rest, not after the raise below: `build` exists to carry *every*
    # problem found so the UI can mark all the bad fields at once.
    roughness = _number(form.roughness, "roughness", errors, minimum=0.0)
    coating = form.coating.strip() or None
    if coating is not None and coating not in _material_names():
        known = ", ".join(_material_names()) or "none installed"
        errors.append(
            FieldError("coating", f"unknown material {coating!r}; available: {known}")
        )

    layers = _build_layers(form.layers, errors)

    if errors:
        raise FormErrors(tuple(errors))

    assert period is not None and profile is not None and illumination is not None
    assert start is not None and stop is not None and count is not None
    assert roughness is not None

    problem = Problem(
        period=period,
        profile=profile,
        coating=coating,
        roughness=roughness,
        layers=layers,
    )
    wavelengths = np.linspace(start, stop, int(count))

    return Parsed(problem=problem, illumination=illumination, wavelengths=wavelengths)


def _build_layers(text: str, errors: list[FieldError]) -> tuple[ConformalLayer, ...]:
    """Parse ``"Au 20, Cr 5"`` into conformal layers, innermost first.

    Every bad entry is reported, in one error on the ``layers`` field that
    names each, so a long stack is fixed in one pass. A material must be in
    the library, as the coating must.
    """
    entries = [part.strip() for part in text.split(",") if part.strip()]
    layers: list[ConformalLayer] = []
    problems: list[str] = []
    known = None
    for position, entry in enumerate(entries, start=1):
        words = entry.split()
        if len(words) != 2:
            problems.append(
                f"layer {position} ({entry!r}) should be a material and a "
                "thickness in nm, e.g. 'Au 20'"
            )
            continue
        name, thickness_text = words
        try:
            thickness = float(thickness_text)
        except ValueError:
            problems.append(f"layer {position}: {thickness_text!r} is not a number")
            continue
        if not np.isfinite(thickness) or thickness <= 0.0:
            problems.append(f"layer {position}: the thickness must be positive")
            continue
        if known is None:
            known = _material_names()
        if name not in known:
            problems.append(
                f"layer {position}: unknown material {name!r}; available: "
                + (", ".join(known) or "none installed")
            )
            continue
        layers.append(ConformalLayer(material=name, thickness=thickness))
    if problems:
        errors.append(FieldError("layers", "; ".join(problems)))
    return tuple(layers)


def _material_names() -> tuple[str, ...]:
    """Vendored material names, for validating the coating field.

    The import is made inside the function so that importing `state`, which
    the headless tests do constantly, does not touch the data directory, and
    so that a trimmed `materials/data/` produces "none installed" instead of
    an import error.
    """
    from ..materials import available

    return available()


def validate(form: FormState) -> tuple[FieldError, ...]:
    """Errors in a form, without raising, for live feedback as fields change."""
    try:
        build(form)
    except FormErrors as exc:
        return exc.errors
    return ()
