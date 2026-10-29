"""What the Compare tab draws and says, without a GUI toolkit.

The Compare tab runs one geometry through several solvers and puts the
answers on one pair of axes: efficiency above, the difference from a chosen
reference solver below. Like :mod:`gratinglab.gui.provenance` and
:mod:`gratinglab.gui.efficiency_plot`, the decisions live here, so they can be
tested without a window and the widget only paints them.

Colour means the order and line style means the solver, never the other way
round. The order colours are the ones every solver tab already uses
(:func:`~gratinglab.gui.theme.order_colors`), so m = +1 is the same orange on
the Compare tab as on the Scalar tab. A solver keeps its line style whichever
is the reference and whichever others are ticked, because the style is fixed
by the solver's name rather than by its position in the run.

A solver that refused the problem, or failed on it, is not dropped from the
picture without a word. :class:`ComparisonOutcome` carries every solver that
was asked for, and :func:`comparison_lines` names each one that produced
nothing, with the reason it gave.
"""

from __future__ import annotations

from collections.abc import Collection, Sequence
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from .provenance import Line
from .theme import Theme, order_colors, set_titles, style_axes

if TYPE_CHECKING:  # pragma: no cover
    from ..checks import EnergyReport
    from ..compare import Comparison

__all__ = [
    "VERIFICATION_NOTE",
    "ComparisonOutcome",
    "line_style",
    "pick_reference",
    "draw_comparison",
    "comparison_lines",
]

#: Line styles in the order solvers are registered. Four solvers, four
#: styles: solid, dashed, dash-dot and dotted all stay distinct at 1.3 pt,
#: which no fifth matplotlib style does reliably.
_STYLES = ("-", "--", "-.", ":")

#: Closes every comparison the panel reports. Two solvers agreeing is evidence
#: that neither has a bug the other lacks, not that either is right: they can
#: share an assumption, a material table or a geometry convention. So the
#: panel says so where agreement is shown, dimmed because it is a fact about
#: the methods, not a fault in this run.
VERIFICATION_NOTE = (
    "\nagreement between solvers is not verification: they share the problem, "
    "the optical constants and the conventions, so check against measurement "
    "or a published benchmark\n"
)

#: More legend entries than this and the order legend takes a second column,
#: as the solver tabs' legend does.
_LEGEND_ROWS = 20


@dataclass(frozen=True)
class ComparisonOutcome:
    """Everything one Compare run produced, across the thread boundary.

    ``requested`` is every solver that was ticked, in run order. Each of them
    ends up in exactly one of ``comparison.methods``, ``refused`` or
    ``failed``. ``comparison`` is ``None`` when no solver produced a scan.
    ``refused`` holds solvers whose capabilities declined the problem
    (an ``UnsupportedConfiguration``, or a form error on that solver's own
    tab), ``failed`` those that raised anything else; both map a solver name
    to the reason it gave.
    """

    requested: tuple[str, ...]
    comparison: "Comparison | None"
    energies: dict[str, "EnergyReport"] = field(default_factory=dict)
    refused: dict[str, str] = field(default_factory=dict)
    failed: dict[str, str] = field(default_factory=dict)

    @property
    def solved(self) -> tuple[str, ...]:
        return () if self.comparison is None else self.comparison.methods


def line_style(name: str, solvers: Sequence[str]) -> str:
    """The line style of solver ``name``, fixed by its place in ``solvers``.

    Pass every registered solver, not only the ticked ones, so that unticking
    one never restyles the rest.
    """
    try:
        index = list(solvers).index(name)
    except ValueError:
        index = len(solvers)
    return _STYLES[index % len(_STYLES)]


def pick_reference(wanted: str | None, solved: Sequence[str]) -> str | None:
    """The reference actually used: ``wanted`` if it solved, else the first
    rigorous solver that did, else the first that did.

    A reference that refused cannot be subtracted from anything, and falling
    back silently would make the lower panel mean something other than its
    title says, so the caller reports the substitution (see
    :func:`comparison_lines`).
    """
    if not solved:
        return None
    if wanted in solved:
        return wanted
    from ..solvers.base import get_solver

    for name in solved:
        try:
            if get_solver(name).capabilities.rigorous:
                return name
        except KeyError:
            continue
    return solved[0]


def _title(name: str) -> str:
    from .docs import short_title

    return short_title(name)


def draw_comparison(
    efficiency_axes,
    difference_axes,
    comparison: "Comparison",
    reference: str,
    visible: Collection[int],
    theme: Theme,
    *,
    solvers: Sequence[str],
    log_difference: bool = False,
) -> dict[int, str]:
    """Draw ``comparison`` into the two axes; return the colour of each order.

    The upper axes get each solver's efficiency in every visible order. The
    lower get ``solver - reference`` for every other solver, or ``|solver -
    reference|`` on a log scale when ``log_difference`` is set, which is what
    makes a 1e-7 agreement between rigorous methods visible beside a 0.1
    disagreement with the scalar one. Exact zeros are left off a log axis
    rather than drawn at an invented floor.
    """
    colors = order_colors((int(m) for m in comparison.orders), theme)
    wavelengths = comparison.wavelengths
    shown = [
        (column, int(order))
        for column, order in enumerate(comparison.orders)
        if int(order) in visible
    ]
    ref_index = comparison.methods.index(reference)

    for axes in (efficiency_axes, difference_axes):
        axes.clear()
        style_axes(axes, theme)
        axes.grid(alpha=0.35, color=theme.grid)

    for index, name in enumerate(comparison.methods):
        style = line_style(name, solvers)
        for column, order in shown:
            efficiency_axes.plot(
                wavelengths, comparison.efficiency[index, :, column],
                color=colors[order], ls=style, lw=1.3,
            )
        if index == ref_index:
            continue
        delta = comparison.efficiency[index] - comparison.efficiency[ref_index]
        for column, order in shown:
            values = delta[:, column]
            if log_difference:
                values = np.abs(values)
                values = np.where(values > 0.0, values, np.nan)
            difference_axes.plot(
                wavelengths, values, color=colors[order], ls=style, lw=1.3,
            )

    efficiency_axes.axhline(1.0, color=theme.bad, lw=0.8, ls=":", alpha=0.8)
    efficiency_axes.set_ylabel("efficiency")
    efficiency_axes.set_ylim(bottom=0)
    count = len(comparison.methods)
    set_titles(
        efficiency_axes, "Efficiency by solver",
        f"{count} solver{'s' if count != 1 else ''}, "
        f"{len(shown)} of {len(comparison.orders)} orders shown",
        theme,
    )

    if log_difference:
        difference_axes.set_yscale("log")
        difference_axes.set_ylabel("|Δ efficiency|")
    else:
        difference_axes.set_yscale("linear")
        difference_axes.axhline(0.0, color=theme.ink, lw=0.8, alpha=0.6)
        difference_axes.set_ylabel("Δ efficiency")
    difference_axes.set_xlabel("wavelength (nm)")
    set_titles(
        difference_axes, f"Difference from {_title(reference)}",
        _worst_detail(comparison, reference, shown), theme,
    )

    _legends(
        efficiency_axes, difference_axes, comparison, shown, colors, theme,
        solvers, reference,
    )
    return colors


def _worst_detail(comparison, reference: str, shown) -> str:
    """The largest disagreement over the drawn orders, for the title's right."""
    others = [m for m in comparison.methods if m != reference]
    if not others:
        return "nothing to compare with"
    if not shown:
        return "no orders shown"
    columns = [c for c, _ in shown]
    ref = comparison.efficiency[comparison.methods.index(reference)][:, columns]
    worst_name, worst = None, -1.0
    for name in others:
        delta = comparison.efficiency[comparison.methods.index(name)][:, columns] - ref
        value = float(np.abs(delta).max())
        if value > worst:
            worst_name, worst = name, value
    return f"largest |Δ| {worst:.2g} ({_title(worst_name)})"


def _legends(axes, difference_axes, comparison, shown, colors, theme, solvers, reference) -> None:
    """Two keys outside the axes: orders by colour above, solvers by style below.

    Two keys rather than one entry per curve, because the curves are the
    product of the two and a 4 x 10 legend would be forty entries saying the
    same thing forty times. One key per axes, rather than both on one, so the
    figure's constrained layout makes room for each: a second legend added
    to the same axes is invisible to the layout and was clipped.
    """
    from matplotlib.lines import Line2D

    keys = []
    if shown:
        keys.append(axes.legend(
            [Line2D([], [], color=colors[m], lw=2.0) for _, m in shown],
            [f"m={m:+d}" for _, m in shown],
            fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0, frameon=False,
            ncol=1 if len(shown) <= _LEGEND_ROWS else 2,
        ))
    keys.append(difference_axes.legend(
        [Line2D([], [], color=theme.ink, ls=line_style(name, solvers), lw=1.3)
         for name in comparison.methods],
        [_title(name) + (" (ref)" if name == reference else "")
         for name in comparison.methods],
        fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0, frameon=False, handlelength=3.0,
    ))
    for legend in keys:
        for text in legend.get_texts():
            text.set_color(theme.ink)


def comparison_lines(
    outcome: ComparisonOutcome,
    wanted_reference: str | None,
    *,
    cancelled: bool = False,
) -> tuple[Line, ...]:
    """What the Compare tab's panel says about a run, solver by solver.

    One line per requested solver, in run order: the time it took and its
    energy verdict if it solved, its reason if it refused or failed. Then the
    headline numbers against the reference, over every order, so the panel
    states the disagreement the plot shows. A refusal is reported in the
    ``warn`` colour, because it is the solver declining a problem it cannot
    model, which is correct behaviour and not a fault; a failure is ``bad``.
    """
    lines: list[Line] = []
    if cancelled:
        lines.append(Line(
            "cancelled: the comparison stopped; showing the previous one\n", "dim",
        ))
    for name in outcome.requested:
        lines.append(Line(f"{_title(name)}: ", None))
        if name in outcome.refused:
            lines.append(Line(f"refused: {outcome.refused[name]}\n", "warn"))
        elif name in outcome.failed:
            lines.append(Line(f"failed: {outcome.failed[name]}\n", "bad"))
        elif name in outcome.solved:
            lines.extend(_solved_line(name, outcome))
        else:
            lines.append(Line("not run\n", "dim"))

    reference = pick_reference(wanted_reference, outcome.solved)
    if reference is None:
        lines.append(Line("\nno solver produced a result to compare\n", "bad"))
        return tuple(lines)
    if wanted_reference is not None and reference != wanted_reference:
        lines.append(Line(
            f"\nreference: {_title(wanted_reference)} did not solve, so the "
            f"differences are from {_title(reference)}\n", "warn",
        ))
    else:
        lines.append(Line(f"\nreference: {_title(reference)}\n", "dim"))

    comparison = outcome.comparison
    for name in comparison.methods:
        if name == reference:
            continue
        summary = comparison.summary(name, reference)
        lines.append(Line(
            f"  {_title(name)} − {_title(reference)}: max |Δ| "
            f"{summary['max_abs_difference']:.3g} at "
            f"{summary['at_wavelength']:.4g} nm, m={summary['at_order']:+d}; "
            f"rms {summary['rms_difference']:.3g}\n", None,
        ))
    lines.append(Line(VERIFICATION_NOTE, "dim"))
    return tuple(lines)


def _solved_line(name: str, outcome: ComparisonOutcome) -> list[Line]:
    scan = outcome.comparison.scans[outcome.comparison.methods.index(name)]
    provenance = scan.provenance
    wall = provenance.wall_time_s
    timing = "time not recorded" if wall is None else f"{wall * 1e3:.0f} ms"
    parts = [Line(f"solved, {timing}", "dim")]
    if not _rigorous(name):
        parts.append(Line(", approximate method", "dim"))
    energy = outcome.energies.get(name)
    if energy is not None:
        parts.append(Line(", energy balance ", "dim"))
        if energy.passed:
            parts.append(Line("ok", "ok"))
        elif energy.lossless:
            deviation = float(np.abs(energy.total - 1.0).max())
            parts.append(Line(f"off unity by up to {deviation:.2g}", "bad"))
        else:
            parts.append(Line(f"EXCEEDS UNITY by {energy.max_excess:.2g}", "bad"))
    parts.append(Line("\n"))
    # The solver's own warnings, word for word and in the solver tabs' form.
    # A verdict alone says that something is wrong; the warning is where the
    # solver says why and what to change. Without it, the integral solver on
    # an equal-arc mesh read "EXCEEDS UNITY by 2.3" here while its own tab
    # said "pass corner_grading=4".
    parts.extend(Line(f"    ⚠ {warning}\n", "warn") for warning in provenance.warnings)
    return parts


def _rigorous(name: str) -> bool:
    from ..solvers.base import get_solver

    try:
        return get_solver(name).capabilities.rigorous
    except KeyError:
        return True
