"""Decide what the provenance panel says, without a GUI toolkit.

The provenance panel is what makes the GUI more than a demonstration, because
it is where a result states how far it can be trusted. The choice of colour
for each line is therefore a question of correctness, and questions of
correctness do not belong in a widget.

This module makes that decision and returns :class:`Line` tuples, and the
widget layer only paints them. The rule the tests protect is that output which
is not wrong must not be styled as wrong. ``converged is None`` means "not yet
checked", which was true for every solver until the convergence harness was
added, and a coating-free run giving relative efficiency is the correct
default. Both were once rendered identically to a real problem, which made a
fully successful run appear broken. The convergence harness now exists, so
``None`` means "this solve did not sweep" instead of "nothing can", and it is
still not a failure.
"""

from __future__ import annotations

import html
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Sequence

if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ..checks import EnergyReport
    from ..result import EfficiencyScan
    from .state import FieldError

__all__ = [
    "Line",
    "TAG_COLORS",
    "provenance_lines",
    "error_lines",
    "to_html",
    "about_text",
]


@dataclass(frozen=True, slots=True)
class Line:
    """One run of panel text and the tag that colours it.

    A ``tag`` of ``None`` is ordinary body text. Lines are emitted in reading
    order and a single visual line may be several ``Line`` objects, so that a
    neutral label can introduce a coloured value, for example
    ``"energy balance: "`` dimmed, followed by the range in green or red.
    """

    text: str
    tag: str | None = None


#: Only four, and each has one meaning. ``dim`` is for facts about
#: the run, ``ok``/``bad`` for a verdict that was reached, and ``warn``
#: for a validity guard the solver itself raised.
TAG_COLORS: dict[str, str] = {
    "warn": "#a5370d",
    "bad": "#b00020",
    "ok": "#0a7d28",
    "dim": "#666666",
}

#: Tags that read as "something is wrong here". The invariant this module
#: exists to protect is about which text may carry one.
ALARM_TAGS = frozenset({"warn", "bad"})


def provenance_lines(
    scan: "EfficiencyScan",
    energy: "EnergyReport",
    lambda_over_period: float,
    *,
    cancelled: bool = False,
) -> tuple[Line, ...]:
    """Everything the panel says about one completed solve.

    ``lambda_over_period`` is passed in instead of a parsed form, so this
    function needs nothing from the GUI and is easy to call from a test.

    ``cancelled`` marks the case where the user asked to stop and the result on
    screen is the previous one. The solve now does stop, because the solver
    checks for cancellation at every wavelength, so the text no longer needs
    to allow for work continuing in the background.
    """
    provenance = scan.provenance
    lines: list[Line] = []

    truncation = provenance.truncation
    knob = (
        f"{truncation} {_knob_label(provenance.method)}  ·  "
        if truncation is not None
        else ""
    )
    lines.append(
        Line(
            f"{provenance.method} {provenance.version}  ·  "
            f"{knob}"
            f"{_milliseconds(provenance.wall_time_s)}  ·  "
            f"λ/period = {lambda_over_period:.4g}\n"
        )
    )

    if cancelled:
        # It used to say the calculation was "still finishing", which was true
        # then and would be a lie now: the solver checks for cancellation at
        # every wavelength and raises out of its own loop.
        lines.append(
            Line(
                "cancelled: the calculation stopped; showing the previous "
                "result\n",
                "dim",
            )
        )

    lines.append(Line("convergence: ", "dim"))
    study = provenance.notes.get("convergence")
    if provenance.converged is None:
        # Not a failure. `gratinglab.convergence` exists, but a plain solve
        # does not run it -- a sweep costs several solves and is the caller's
        # decision, not a hidden tax on pressing Solve. So this stays the
        # ordinary state and must not be styled as alarming.
        lines.append(Line("not checked; see gratinglab.convergence\n", "dim"))
    elif provenance.converged:
        lines.append(Line("yes", "ok"))
        lines.append(Line(_evidence(study), "dim"))
    else:
        lines.append(Line("NO", "bad"))
        lines.append(Line(_evidence(study), "dim"))

    normalization = provenance.notes.get("normalization")
    if normalization:
        # The detail comes from whether a coating was *resolved*, which is a
        # different question from whether its reflectivity was applied -- so
        # "relative (Au ...)" is a real and honest state, meaning the material
        # is known and not yet used. Deriving the wording from `normalization`
        # alone would report "no coating" for it.
        # Two keys because two solvers mean two things by it: scalar's
        # `coating` is the material whose Fresnel amplitude was applied, the
        # integral solver's `material` is the one its boundary condition read.
        # Either is the honest detail; neither present means none was used.
        detail = (
            provenance.notes.get("coating")
            or provenance.notes.get("material")
            or "no coating"
        )
        lines.append(Line("normalization: ", "dim"))
        lines.append(Line(f"{normalization} ({detail})\n", "dim"))

    if scan.absorption is not None:
        # A physical result, not a problem: `dim`, never an alarm tag. It also
        # explains the deficit in Σ below, which would otherwise read as a
        # method losing power it never lost.
        lines.append(Line("absorption: ", "dim"))
        lines.append(
            Line(
                f"A ∈ [{scan.absorption.min():.4f}, "
                f"{scan.absorption.max():.4f}]\n",
                "dim",
            )
        )

    # `check_energy_balance` folds absorption into the total when a scan
    # carries one, so the conserved quantity is R + A there. Naming it Σ
    # anyway would collide with the Σ curve on the plot, which is R alone.
    conserved = "Σ + A" if scan.absorption is not None else "Σ"
    lines.append(Line("energy balance: ", "dim"))
    lines.append(
        Line(
            f"{conserved} ∈ [{energy.total.min():.4f}, {energy.total.max():.4f}]"
            f"{'' if energy.passed else '  (EXCEEDS UNITY)'}\n",
            "ok" if energy.passed else "bad",
        )
    )

    for warning in provenance.warnings:
        lines.append(Line(f"  ⚠ {warning}\n", "warn"))

    return tuple(lines)


def _knob_label(method: str) -> str:
    """What this method's truncation count counts.

    The label is read from the registry instead of from a table kept here, so
    a solver added later names its own knob without any edit to this module.
    The header previously said "quadrature pts" for every method, which was
    wrong for the integral solver's boundary points from the time that tab
    was added. A fallback is returned for imported reference data, whose
    `method` need not be a registered solver at all.
    """
    from ..solvers.base import get_solver

    try:
        knob = get_solver(method).capabilities.accuracy_knob
    except (KeyError, ImportError):
        return "truncation"
    return knob.replace("_", " ") if knob else "truncation"


def _evidence(study: Any) -> str:
    """What a convergence verdict rests on, in one dim trailing clause.

    A bare "yes" is only a claim; the knob value and the swept range allow a
    reader to check it. `converged_at` is the more useful of the two for
    acting on the result, since it is usually cheaper than the value needed
    to establish it. For a verdict that arrived without its study, the
    function returns a bare newline and does not invent detail.
    """
    if not isinstance(study, dict):
        return "\n"
    knob = study.get("knob", "accuracy")
    at, values = study.get("converged_at"), study.get("values") or []
    swept = f", swept to {values[-1]}" if values else ""
    if at is None:
        return f"; {knob} still moving{swept}\n"
    return f"; {knob}={at} is enough{swept}\n"


def error_lines(errors: Sequence["FieldError"]) -> tuple[Line, ...]:
    """Report which fields stopped the solve, the panel's other function."""
    lines = [Line(f"{len(errors)} field(s) need attention\n", "bad")]
    lines.extend(Line(f"  • {e.field}: {e.message}\n", "bad") for e in errors)
    return tuple(lines)


def solving_lines(method: str, wavelength_count: int) -> tuple[Line, ...]:
    """Shown while a solve is in flight, in place of a stale verdict."""
    return (
        Line(f"solving: {method}, {wavelength_count} wavelengths…\n", "dim"),
    )


def to_html(lines: Sequence[Line]) -> str:
    """Render lines as the small HTML subset Qt's rich text engine accepts.

    The text is escaped, so a solver warning containing ``<`` or ``&`` (the
    docs do contain ampersands) cannot silently become markup or disappear.
    """
    parts: list[str] = []
    for line in lines:
        body = html.escape(line.text).replace("\n", "<br>")
        color = TAG_COLORS.get(line.tag) if line.tag else None
        parts.append(f'<span style="color:{color}">{body}</span>' if color else body)
    return f'<div style="white-space:pre-wrap">{"".join(parts)}</div>'


def alarming_text(lines: Sequence[Line]) -> str:
    """Every character styled as a problem, concatenated.

    It exists for the tests that assert that something is not flagged. Such
    an assertion is otherwise awkward to express, and it was previously made
    by walking Tk tag index pairs.
    """
    return "".join(line.text for line in lines if line.tag in ALARM_TAGS)


def about_text(version: str) -> str:
    """The About box, as a string so it can be asserted without a dialog."""
    return (
        "GratingLab\n"
        f"Version {version}\n\n"
        "An open comparison platform for grating efficiency analysis: "
        "scalar, RCWA, C-method and integral-method solvers driven from "
        "one problem spec.\n\n"
        "BSD-3-Clause."
    )


def _milliseconds(wall_time_s: Any) -> str:
    """Wall time for the header, tolerant of a solver that recorded none."""
    if wall_time_s is None:
        return "time not recorded"
    return f"{wall_time_s * 1e3:.0f} ms"
