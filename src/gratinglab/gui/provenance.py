"""What the provenance panel says, decided without a toolkit.

The panel is the reason the GUI is worth having rather than a toy: it is where
a result declares how far it can be trusted. That makes *which line gets which
colour* a correctness question, not a styling one, and correctness questions do
not belong in a widget.

So this module owns the decision and returns :class:`Line` tuples; the widget
layer only paints them. The rule the tests exist to protect is that **things
which are not wrong must not look wrong**: ``converged is None`` means "not yet
checked", true for every solver until the convergence harness lands, and a
coating-free run giving relative efficiency is the correct default. Both once
rendered identically to a real problem, which is how a fully successful run
could read as broken.
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
    neutral label can introduce a coloured value -- ``"energy balance: "``
    dimmed, then the range in green or red.
    """

    text: str
    tag: str | None = None


#: Only four, and their meanings are load-bearing. ``dim`` is for facts about
#: the run, ``ok``/``bad`` for a verdict that was actually reached, and ``warn``
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

    ``lambda_over_period`` is passed in rather than a parsed form, so this
    function needs nothing from the GUI and is trivial to call from a test.

    ``cancelled`` marks the case where the user asked to stop and the result on
    screen is the *previous* one. The solve itself cannot be interrupted -- see
    ``gui/qt/worker.py`` -- so saying "cancelled" without saying what is still
    running would overstate what happened.
    """
    provenance = scan.provenance
    lines: list[Line] = []

    truncation = provenance.truncation
    lines.append(
        Line(
            f"{provenance.method} {provenance.version}  ·  "
            f"{truncation} quadrature pts  ·  "
            f"{_milliseconds(provenance.wall_time_s)}  ·  "
            f"λ/period = {lambda_over_period:.4g}\n"
        )
    )

    if cancelled:
        lines.append(
            Line(
                "cancelled — showing the previous result; the calculation you "
                "stopped is still finishing\n",
                "dim",
            )
        )

    lines.append(Line("convergence: ", "dim"))
    if provenance.converged is None:
        # Not a failure. No solver has a convergence harness yet, so this is
        # the honest, universal state -- and must not be styled as alarming.
        lines.append(Line("not yet checked\n", "dim"))
    elif provenance.converged:
        lines.append(Line("yes\n", "ok"))
    else:
        lines.append(Line("NO\n", "bad"))

    normalization = provenance.notes.get("normalization")
    if normalization:
        detail = "no coating" if normalization == "relative" else "coating set"
        lines.append(Line("normalization: ", "dim"))
        lines.append(Line(f"{normalization} ({detail})\n", "dim"))

    lines.append(Line("energy balance: ", "dim"))
    lines.append(
        Line(
            f"Σ ∈ [{energy.total.min():.4f}, {energy.total.max():.4f}]"
            f"{'' if energy.passed else '  — EXCEEDS UNITY'}\n",
            "ok" if energy.passed else "bad",
        )
    )

    for warning in provenance.warnings:
        lines.append(Line(f"  ⚠ {warning}\n", "warn"))

    return tuple(lines)


def error_lines(errors: Sequence["FieldError"]) -> tuple[Line, ...]:
    """The panel's other job: say which fields stopped the solve."""
    lines = [Line(f"{len(errors)} field(s) need attention\n", "bad")]
    lines.extend(Line(f"  • {e.field}: {e.message}\n", "bad") for e in errors)
    return tuple(lines)


def solving_lines(method: str, wavelength_count: int) -> tuple[Line, ...]:
    """Shown while a solve is in flight, in place of a stale verdict."""
    return (
        Line(f"solving — {method}, {wavelength_count} wavelengths…\n", "dim"),
    )


def to_html(lines: Sequence[Line]) -> str:
    """Render lines as the small HTML subset Qt's rich text engine accepts.

    Escaped, so a solver warning containing ``<`` or ``&`` -- and the docs do
    carry ampersands -- cannot silently become markup or vanish.
    """
    parts: list[str] = []
    for line in lines:
        body = html.escape(line.text).replace("\n", "<br>")
        color = TAG_COLORS.get(line.tag) if line.tag else None
        parts.append(f'<span style="color:{color}">{body}</span>' if color else body)
    return f'<div style="white-space:pre-wrap">{"".join(parts)}</div>'


def alarming_text(lines: Sequence[Line]) -> str:
    """Every character styled as a problem, concatenated.

    Exists for the tests that assert something is *not* flagged, which is
    otherwise awkward to phrase and was previously done by walking Tk tag
    index pairs.
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
