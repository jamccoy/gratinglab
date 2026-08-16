"""The app icon, drawn rather than stored.

Rendered from matplotlib because this machine has no SVG rasteriser
(``rsvg-convert`` / ``inkscape``), and generating the artwork keeps it
reproducible and editable in the same language as the rest of the project.

The motif is a blazed sawtooth with an AFM cantilever and tip riding across it --
the measurement this software performs. Deliberately few, thick elements so it
stays legible at 16 px, which is where most icons fail.

Sibling project GratingLab uses a sawtooth with diffracted orders fanning out of
it. Keeping the tip here rather than the orders is what distinguishes the two in
a Dock at small sizes: this one measures the groove, that one models what the
groove does to light.

matplotlib is imported inside ``render`` rather than at module level: importing
``afm_analysis.gui`` must stay cheap, since that is what raises the friendly
"PySide6 is an optional extra" message. ``tests/test_gui_boundary.py`` checks it.

``tools/make_icon.py`` turns this into a macOS ``.icns``.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["render", "BACKGROUND", "FACET", "TIP", "SUBSTRATE"]

#: Warm near-black, so the pale facets and the tip read against it at any size.
BACKGROUND = "#1b1712"

#: The grating itself: bright, high contrast against the plate.
FACET = "#f2e9dc"

#: Cantilever and tip. Amber against the cool facets, so the two shapes stay
#: separable when the downsample blurs their edges together.
TIP = "#ffb03a"

#: Shadowed under-face of each tooth, giving the sawtooth its asymmetry at a
#: glance. A flat silhouette reads as a triangle wave rather than a *blazed*
#: grating, which is the whole point of the motif.
SUBSTRATE = "#8d7a63"


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

    # Two big blazed teeth. Three read better at 512 px but disappear at 16; two
    # survive the downsample and still say "blazed grating". The long face rises
    # gently and the short face drops vertically -- that asymmetry is the blaze.
    baseline, height = 0.17, 0.30
    left, width = 0.13, 0.37
    for tooth in range(2):
        start = left + tooth * width
        peak_x = start + width * 0.86
        axes.add_patch(
            Polygon(
                [(start, baseline), (peak_x, baseline + height), (peak_x, baseline)],
                closed=True, facecolor=FACET, edgecolor="none",
            )
        )
        # Thin shadow along the steep face, so the drop is visible even when the
        # facet and the plate blur into each other.
        axes.add_patch(
            Polygon(
                [(peak_x, baseline), (peak_x, baseline + height),
                 (peak_x - width * 0.07, baseline)],
                closed=True, facecolor=SUBSTRATE, edgecolor="none",
            )
        )

    # The tip rides on the *blaze* facet -- the long rising face of the first
    # tooth, which is the face this software measures. Placing the apex on that
    # line rather than near it is what makes the icon read as a measurement
    # rather than as two unrelated shapes.
    facet_start = (left, baseline)
    facet_peak = (left + width * 0.86, baseline + height)
    along = 0.58
    tip_apex = (
        facet_start[0] + along * (facet_peak[0] - facet_start[0]),
        facet_start[1] + along * (facet_peak[1] - facet_start[1]),
    )

    tip_height, tip_half_width = 0.20, 0.075
    tip_top_y = tip_apex[1] + tip_height

    # Cantilever: one straight beam from the upper right down to the tip. A
    # drawn-out cantilever outline vanishes at 16 px; a single thick stroke does
    # not. Drawn before the tip so the wedge sits cleanly on top of it.
    axes.plot(
        [0.90, tip_apex[0]], [0.78, tip_top_y],
        color=TIP, lw=size * 0.030, solid_capstyle="round",
    )

    # The tip: a symmetric downward wedge, apex on the facet. An asymmetric
    # wedge looked like a crease in the beam rather than a probe.
    axes.add_patch(
        Polygon(
            [(tip_apex[0], tip_apex[1]),
             (tip_apex[0] - tip_half_width, tip_top_y),
             (tip_apex[0] + tip_half_width, tip_top_y)],
            closed=True, facecolor=TIP, edgecolor="none",
        )
    )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=100, facecolor="none", transparent=True)
    return path
