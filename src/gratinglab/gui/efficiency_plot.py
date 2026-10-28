"""The efficiency-against-wavelength plot both solver tabs draw.

The scalar and integral tabs carried the same drawing code twice, and a
theme would have had to be threaded through both copies. It lives here once,
without Qt, so it can be drawn into any matplotlib axes, including one in a
test.

The legend sits outside the axes, on the right. Inside, it covered the
curves it was naming: at the reference geometry the top-right corner is
where the low orders peak.
"""

from __future__ import annotations

from collections.abc import Collection
from typing import TYPE_CHECKING

from .theme import Theme, order_colors, style_axes

if TYPE_CHECKING:  # pragma: no cover
    from ..result import EfficiencyScan

__all__ = ["draw_efficiency"]

#: More entries than this and the legend takes a second column.
_LEGEND_ROWS = 20


def draw_efficiency(
    axes, scan: "EfficiencyScan", visible: Collection[int], theme: Theme
) -> dict[int, str]:
    """Draw ``scan`` into ``axes`` and return the colour used for each order.

    Orders not in ``visible`` are not drawn, but they still take part in the
    colour assignment, so toggling one never recolours the rest. The returned
    mapping is what the order list's swatches show.
    """
    colors = order_colors((int(m) for m in scan.orders), theme)
    axes.clear()
    style_axes(axes, theme)
    for index, order in enumerate(scan.orders):
        order = int(order)
        if order not in visible:
            continue
        axes.plot(
            scan.wavelengths, scan.efficiency[:, index],
            color=colors[order], lw=1.3, label=f"m={order:+d}",
        )
    axes.plot(
        scan.wavelengths, scan.total,
        color=theme.ink, ls="--", lw=1.0, alpha=0.7, label="Σ",
    )
    transmitted = (
        scan.transmission.total if scan.transmission is not None else 0.0
    )
    if scan.transmission is not None:
        axes.plot(
            scan.wavelengths, scan.transmission.total,
            color=theme.ink, lw=1.0, ls=":", alpha=0.9, label="T",
        )
    if scan.absorption is not None:
        # Under a finite-conductivity solve Σ alone is the reflected fraction
        # and sits below the unity line by exactly the absorbed one. Drawing
        # A and Σ + A puts the R + A = 1 theorem on the same axes as its own
        # discretisation error, rather than leaving a correct result looking
        # like a deficit.
        axes.plot(
            scan.wavelengths, scan.absorption,
            color=theme.warn, lw=1.0, ls="-.", alpha=0.9, label="A",
        )
        axes.plot(
            scan.wavelengths, scan.total + transmitted + scan.absorption,
            color=theme.ink, lw=1.0, ls="-", alpha=0.5,
            label="Σ + T + A" if scan.transmission is not None else "Σ + A",
        )
    axes.axhline(1.0, color=theme.bad, lw=0.8, ls=":", alpha=0.8)
    axes.set_xlabel("wavelength (nm)")
    axes.set_ylabel("efficiency")
    axes.set_ylim(bottom=0)
    axes.grid(alpha=0.35, color=theme.grid)
    handles, _ = axes.get_legend_handles_labels()
    if handles:
        legend = axes.legend(
            fontsize=7,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.0),
            borderaxespad=0.0,
            ncol=1 if len(handles) <= _LEGEND_ROWS else 2,
            frameon=False,
        )
        for text in legend.get_texts():
            text.set_color(theme.ink)
    return colors
