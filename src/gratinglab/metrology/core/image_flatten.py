"""
Image flattening: removing background from the 2-D scan, before rows are averaged.

This is the *first* of two flattening stages, and the two are easy to confuse:

- **Image flattening** (here) acts on the 2-D array as it came off the
  instrument, correcting scan lines relative to one another.
- **Profile flattening** (``core.processing.flatten_profile``) acts on the 1-D
  height trace produced after each row group has been averaged.

Both are "flattening" in the usual AFM sense. The distinction is what they see.

**Affine image flattening cannot change a blaze angle.** Measured, not assumed:
on ``20250820_280C_00004.txt`` - the scan with the worst row-offset spread in
this project's dataset, 2.74 nm - ``none``, ``plane`` and ``align_rows`` all give
30.8524 degrees - identical to four decimal places, and to within 5e-15 in
full precision - and ``align_rows`` gives
0.0000 degrees of change across all eight compare-mode samples.

The reason is structural. ``raw_data_multi_group`` averages rows into a profile,
and ``flatten_profile`` then fits and removes a background from that profile.
Subtracting a per-row constant shifts a band's average by a constant; subtracting
a plane adds a constant and a linear ramp. Every profile-flattening method fits at
least a first-order polynomial, so it removes exactly those terms again. What goes
out here comes out there regardless.

That is worth knowing before benchmarking these and concluding they are broken.
What they *are* for:

- **Seeing the data.** A raw AFM image with sample tilt renders as a gradient
  with the grating barely visible. Levelled, it looks like a grating.
- **Methods that are not affine.** A per-row quadratic changes the shape within
  each row rather than merely its offset or tilt, and does survive the averaging -
  measured at 0.0003 degrees, small but not zero. The registry below is the place
  to add such a method.
- **Being explicit.** The step now happens in a named, recorded place rather than
  in Gwyddion before export, where nothing captured that it had happened at all.
"""
from __future__ import annotations

import numpy as np

__all__ = ["flatten_image", "IMAGE_FLATTEN_METHODS", "VALID_IMAGE_FLATTEN_METHODS",
           "DEFAULT_IMAGE_FLATTEN_METHOD", "row_offset_spread"]

DEFAULT_IMAGE_FLATTEN_METHOD = "align_rows"


def _none(data: np.ndarray) -> np.ndarray:
    """Leave the image alone. The pre-2026-08 behaviour."""
    return data


def _plane(data: np.ndarray) -> np.ndarray:
    """
    Subtract a single tilted plane fitted to the whole image.

    Corrects the sample sitting at an angle under the scanner - the most common
    background there is, and the one that makes a raw image unreadable.

    Fitted by least squares over every pixel, so a deep grating biases the fit
    slightly toward its own mean. That does not matter here (see the module
    docstring), but it would if this were ever the only flattening applied.
    """
    rows, cols = data.shape
    y, x = np.mgrid[0:rows, 0:cols]
    basis = np.column_stack([x.ravel(), y.ravel(), np.ones(data.size)])
    coefficients, *_ = np.linalg.lstsq(basis, data.ravel(), rcond=None)
    return data - (basis @ coefficients).reshape(rows, cols)


def _align_rows(data: np.ndarray) -> np.ndarray:
    """
    Put every scan line on a common level by subtracting its median.

    AFM scan lines drift in Z between one another - feedback settling, thermal
    creep - leaving rows offset by a few nanometres. In this project's data that
    spread runs 0.44 to 2.74 nm against groove depths near 100 nm.

    The median rather than the mean: a grating profile is not symmetric about its
    own mean, and a row that happens to catch more groove bottom than land would
    be pulled by a mean where the median barely moves.
    """
    return data - np.median(data, axis=1, keepdims=True)


#: Methods by name. Add one here and it appears in the GUI, the settings
#: validator and the wiki table without further wiring - which is the point.
IMAGE_FLATTEN_METHODS = {
    "none": _none,
    "plane": _plane,
    "align_rows": _align_rows,
}

VALID_IMAGE_FLATTEN_METHODS = tuple(IMAGE_FLATTEN_METHODS)


def flatten_image(data, method: str = DEFAULT_IMAGE_FLATTEN_METHOD) -> np.ndarray:
    """
    Apply one image-flattening method.

    Parameters:
        data: 2-D height array, in metres
        method: a key of IMAGE_FLATTEN_METHODS

    Returns:
        A new array of the same shape. The input is never modified in place -
        callers hold on to the raw image to show a before/after comparison.
    """
    if method not in IMAGE_FLATTEN_METHODS:
        raise ValueError(
            f"unknown image flattening method {method!r}. "
            f"Available: {', '.join(VALID_IMAGE_FLATTEN_METHODS)}")

    data = np.asarray(data, dtype=float)
    if data.ndim != 2:
        raise ValueError(f"expected a 2-D image, got shape {data.shape}")

    return IMAGE_FLATTEN_METHODS[method](data)


def row_offset_spread(data) -> float:
    """
    Standard deviation of the per-row medians - how badly scan lines disagree.

    A one-number summary of what `align_rows` corrects, for the Import tab to
    display so the choice can be made on evidence rather than habit. Zero after
    `align_rows` by construction.
    """
    data = np.asarray(data, dtype=float)
    return float(np.std(np.median(data, axis=1)))
