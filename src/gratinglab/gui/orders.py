"""Which diffraction orders the efficiency plot draws.

Previously one line in the widget layer::

    if values.max() < 1e-4:
        continue

That rule was untested, could not be overridden and was not visible: an
order below the threshold disappeared from the plot with no indication that
it had. This is unacceptable in this project for the same reason that
renormalising away an energy defect is. The threshold is retained as a
default, because it is a sensible one, but each order is now listed, labelled
with the peak it was judged on, and can be shown or hidden by the user.

The logic lives in this module instead of in the panel, so that the default
can be tested without a window and the rule is documented where a reader will
find it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Sequence


if TYPE_CHECKING:  # pragma: no cover - imports for type checkers only
    from ..result import EfficiencyScan

__all__ = [
    "OrderSummary",
    "DEFAULT_THRESHOLD",
    "DEFAULT_LIMIT",
    "summarize",
    "default_visible",
    "carry_over",
    "describe",
]

#: Peak efficiency below which an order is hidden by default. Inherited from
#: the rule it replaces, so the plot looks the same on first open.
DEFAULT_THRESHOLD = 1e-4

#: At most this many curves by default, largest peaks first. New, and the part
#: that protects the plot: an X-ray geometry can propagate hundreds of
#: orders, and the old code's only defence was a legend that grew columns.
DEFAULT_LIMIT = 20


@dataclass(frozen=True, slots=True)
class OrderSummary:
    """One diffraction order, as the visibility panel needs to describe it.

    Attributes
    ----------
    order
        The diffraction order `m`.
    column
        Its column in ``scan.efficiency``. It is carried explicitly so that a
        caller never has to search ``scan.orders`` again and cannot get the
        mapping wrong.
    peak
        Largest efficiency over the scanned wavelengths. The default rule is
        applied to this value, and the panel shows it next to each entry.
    ever_propagating
        Whether the order propagates at any wavelength in the scan. An order
        that is evanescent throughout has nothing to plot. It is listed anyway,
        greyed out, because the fact that the order exists and carries no
        power is itself informative.
    """

    order: int
    column: int
    peak: float
    ever_propagating: bool


def summarize(scan: "EfficiencyScan") -> tuple[OrderSummary, ...]:
    """Describe every order in a scan, in the scan's own order."""
    peaks = scan.efficiency.max(axis=0)
    propagating = scan.propagating.any(axis=0)
    return tuple(
        OrderSummary(
            order=int(order),
            column=column,
            peak=float(peaks[column]),
            ever_propagating=bool(propagating[column]),
        )
        for column, order in enumerate(scan.orders)
    )


def default_visible(
    summaries: Sequence[OrderSummary],
    *,
    threshold: float = DEFAULT_THRESHOLD,
    limit: int = DEFAULT_LIMIT,
) -> frozenset[int]:
    """The orders worth drawing before the user says otherwise.

    Two rules are applied in sequence: the peak must be at or above
    ``threshold``, and then the ``limit`` largest of those are kept. Ties are
    broken by order index so that the result is deterministic, since a plot
    that redrew differently for the same scan would be a bug.

    The function returns a set of orders instead of columns, so the choice
    survives a re-solve that changes how many orders propagate;
    :func:`carry_over` relies on this.
    """
    candidates = [s for s in summaries if s.peak >= threshold]
    candidates.sort(key=lambda s: (-s.peak, s.order))
    return frozenset(s.order for s in candidates[:limit])


def carry_over(
    previous_visible: Iterable[int],
    previous_orders: Iterable[int],
    summaries: Sequence[OrderSummary],
    *,
    threshold: float = DEFAULT_THRESHOLD,
    limit: int = DEFAULT_LIMIT,
) -> frozenset[int]:
    """Reapply a user's choices to a freshly solved scan.

    Adjusting the blaze angle and re-solving must not silently re-check the
    orders that a user deliberately unchecked. An order that existed before
    therefore keeps its previous state, shown or hidden. An order that has
    only just appeared receives the default rule. An order that no longer
    exists is dropped instead of remembered, since a stale selection would
    reappear on the next unrelated re-solve.

    ``previous_orders`` contains every order the last scan had, including the
    hidden ones, and it is what distinguishes an order hidden on purpose from
    one not seen before. With only the visible set the two cases are
    indistinguishable, and an unchecked order would silently re-check itself
    on the next solve.

    Pass empty iterables for the first solve; everything is then new and the
    result is exactly :func:`default_visible`.
    """
    seen_before = frozenset(previous_orders)
    still_exists = {s.order for s in summaries}
    kept = frozenset(previous_visible) & still_exists
    fresh = [s for s in summaries if s.order not in seen_before]
    return kept | default_visible(fresh, threshold=threshold, limit=limit)


def describe(summary: OrderSummary) -> str:
    """One entry's label, e.g. ``m=+3   peak 0.4127``.

    The order is always signed, because its sign is physically meaningful
    (`+3` and `-3` diffract to opposite sides), and an unsigned `3` in a list
    next to `-3` would read as ambiguous instead of positive.
    """
    if not summary.ever_propagating:
        return f"m={summary.order:+d}   evanescent"
    return f"m={summary.order:+d}   peak {summary.peak:.4f}"
