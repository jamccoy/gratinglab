"""Which diffraction orders the efficiency plot draws.

Previously one line in the widget layer::

    if values.max() < 1e-4:
        continue

An untested, unoverridable, invisible rule: an order below the threshold
vanished from the plot with nothing to say it had. That is the wrong shape for
this project -- the same objection as renormalising away an energy defect. The
threshold survives as a *default*, because it is a sensible one, but it is now
listed, labelled with the peak it was judged on, and overridable.

The logic lives here rather than in the panel so the default can be tested
without a window, and so the rule is written down somewhere a reader will find
it.
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
#: that actually protects the plot: an X-ray geometry can propagate hundreds of
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
        Its column in ``scan.efficiency``. Carried explicitly so a caller never
        has to search ``scan.orders`` again, and cannot get the mapping wrong.
    peak
        Largest efficiency over the scanned wavelengths. What the default rule
        judges, and what the panel shows next to each entry.
    ever_propagating
        Whether the order propagates at any wavelength in the scan. An order
        that is evanescent throughout has nothing to plot; it is listed anyway,
        greyed, because "this order exists and carries nothing" is information.
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

    Two rules, in this sequence: peak at or above ``threshold``, then the
    ``limit`` largest of those. Ties are broken by order index so the result is
    deterministic -- a plot that redrew differently for the same scan would be
    its own small bug.

    Returning a set of *orders* rather than columns means the choice survives a
    re-solve that changes how many orders propagate, which is what
    :func:`carry_over` relies on.
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

    Nudging the blaze angle and re-solving must not silently re-check the
    orders someone deliberately unchecked. So: an order that existed before
    keeps whatever it was, shown or hidden; an order that has only just
    appeared gets the default rule; an order that no longer exists is dropped
    rather than remembered, since a stale selection would resurrect on the next
    unrelated re-solve.

    ``previous_orders`` -- every order the last scan had, not just the visible
    ones -- is what separates "hidden on purpose" from "not seen before". With
    only the visible set the two are indistinguishable, and an unchecked order
    would silently re-check itself on the next solve.

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

    Signed, always, because the sign of the order is the physics -- `+3` and
    `-3` diffract to opposite sides -- and an unsigned `3` in a list next to
    `-3` reads as ambiguous rather than positive.
    """
    if not summary.ever_propagating:
        return f"m={summary.order:+d}   evanescent"
    return f"m={summary.order:+d}   peak {summary.peak:.4f}"
