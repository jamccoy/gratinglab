"""The app icon, drawn rather than stored.

Rendered from matplotlib because this machine has no SVG rasteriser
(``rsvg-convert`` / ``inkscape``), and generating the artwork keeps it
reproducible and editable in the same language as the rest of the project.

The motif is a blazed sawtooth with diffracted orders fanning out of it -- the
thing this software computes. Deliberately few, thick elements so it stays
legible at 16 px, which is where most icons fail.

``tools/make_icon.py`` turns this into a macOS ``.icns``.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["render", "BACKGROUND", "FACET", "ORDER_COLORS"]

#: Deep blue-black, so the spectrum reads against it at any size.
BACKGROUND = "#101a2e"
FACET = "#e8ecf5"

#: Long to short wavelength, matching the fan direction. Four rather than six:
#: at 16 px, six strokes blur into an unreadable smear, four stay distinct.
ORDER_COLORS = ("#ff4d4d", "#ffc233", "#4ade80", "#4dc4ff")


def render(path: str | Path, size: int = 1024) -> Path:
    """Draw the icon to a square PNG.

    Parameters
    ----------
    size
        Pixels per side. 1024 is the largest macOS asks for; every smaller
        variant is downsampled from it.
    """
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np
    from matplotlib.figure import Figure
    from matplotlib.patches import FancyBboxPatch, Polygon

    figure = Figure(figsize=(size / 100, size / 100), dpi=100)
    axes = figure.add_axes((0, 0, 1, 1))
    axes.set_xlim(0, 1)
    axes.set_ylim(0, 1)
    axes.axis("off")

    # Rounded-square plate, inset like a macOS app icon.
    axes.add_patch(
        FancyBboxPatch(
            (0.06, 0.06), 0.88, 0.88,
            boxstyle="round,pad=0,rounding_size=0.20",
            facecolor=BACKGROUND, edgecolor="none",
        )
    )

    # Diffracted orders fanning from the blaze point. Drawn before the grating
    # so the facets sit cleanly on top.
    #
    # Everything here is sized for the 16 px case: thick strokes, wide angular
    # separation, few elements. A finer design looks better at 512 and turns to
    # mush in the Finder list view, which is where an icon earns its keep.
    origin = (0.26, 0.42)
    for index, color in enumerate(ORDER_COLORS):
        angle = np.radians(26 + index * 19.0)
        axes.plot(
            [origin[0], origin[0] + 0.62 * np.cos(angle)],
            [origin[1], origin[1] + 0.62 * np.sin(angle)],
            color=color, lw=size * 0.030, solid_capstyle="round",
        )

    # Incident ray, arriving at a shallow graze from the left.
    axes.plot(
        [0.10, origin[0]], [0.60, origin[1]],
        color="#8fa3c8", lw=size * 0.022, solid_capstyle="round",
    )

    # Two big sawtooth teeth. Three read better at 512 px but disappear at 16;
    # two survive the downsample and still say "blazed grating".
    baseline, height = 0.16, 0.26
    left, width = 0.14, 0.36
    for tooth in range(2):
        start = left + tooth * width
        axes.add_patch(
            Polygon(
                [
                    (start, baseline),
                    (start + width * 0.82, baseline + height),
                    (start + width * 0.82, baseline),
                ],
                closed=True, facecolor=FACET, edgecolor="none",
            )
        )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=100, facecolor="none", transparent=True)
    return path
