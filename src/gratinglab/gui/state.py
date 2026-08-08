"""Turning form fields into a validated problem.

This module is the GUI's only interesting logic, and it is deliberately
separated from the widgets so it can be tested headlessly. **Nothing here
imports tkinter**, and nothing here computes physics -- it parses strings,
validates ranges, and hands the result to the tested core.

The split matters: widget code is the least testable part of any project, and
this project's value is defensible numbers. Keeping parsing and validation here
means a bug in the untested layer produces a wrong *layout*, never a wrong
*number*.

Errors are **collected, not raised on the first failure**, so a form with three
bad fields reports three problems rather than making the user fix them one at a
time.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from ..illumination import Illumination
from ..problem import Problem
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

#: Which form fields each profile kind actually uses. The UI hides the rest,
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


@dataclass(frozen=True, slots=True)
class FieldError:
    """One thing wrong with one field."""

    field: str
    message: str

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return f"{self.field}: {self.message}"


class FormErrors(ValueError):
    """Raised by :func:`build` when a form cannot be turned into a problem."""

    def __init__(self, errors: tuple[FieldError, ...]) -> None:
        self.errors = errors
        super().__init__("; ".join(str(e) for e in errors))


@dataclass(frozen=True, slots=True)
class FormState:
    """The grating geometry, exactly as typed. Strings, because that is what a
    form holds.

    Geometry only -- period, profile, mount, wavelengths -- shared across
    every solver a tab might run. A solver's own knobs (scalar's
    ``quadrature_points``, say) live in that solver's own state, e.g.
    :class:`gratinglab.gui.scalar_options.ScalarOptionsState`, not here; this
    dataclass has no way to know which solver is about to use it.

    Defaults describe a soft-X-ray off-plane case, which is the project's
    primary application and exercises the interesting physics immediately.
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

    def with_field(self, name: str, value: Any) -> "FormState":
        """Return a copy with one field changed."""
        return replace(self, **{name: value})


@dataclass(frozen=True, slots=True)
class Parsed:
    """A geometry that made sense, ready to hand to any solver.

    No ``options`` field: geometry is shared across every solver a tab might
    run, but each solver's own options (and their own validation, e.g.
    scalar's Nyquist guard in :mod:`gratinglab.gui.scalar_options`) are that
    solver's business, built separately from whichever tab is about to solve.
    """

    problem: Problem
    illumination: Illumination
    wavelengths: NDArray[np.float64]

    @property
    def lambda_over_period(self) -> float:
        """Ratio at the long-wavelength end -- the scalar validity indicator."""
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
        from ..io.ggp import read_profile

        return read_profile(path)
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
        Carrying **every** problem found, so the UI can mark all bad fields at
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

    if errors:
        raise FormErrors(tuple(errors))

    assert period is not None and profile is not None and illumination is not None
    assert start is not None and stop is not None and count is not None

    problem = Problem(period=period, profile=profile)
    wavelengths = np.linspace(start, stop, int(count))

    return Parsed(problem=problem, illumination=illumination, wavelengths=wavelengths)


def validate(form: FormState) -> tuple[FieldError, ...]:
    """Errors in a form, without raising. For live feedback as fields change."""
    try:
        build(form)
    except FormErrors as exc:
        return exc.errors
    return ()
