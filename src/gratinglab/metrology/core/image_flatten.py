"""
Image flattening: removing background from the 2-D scan, before rows are averaged.

This is the first of two flattening stages, and the two are easily confused:

- **Image flattening** (here) acts on the 2-D array as it came off the
  instrument, correcting scan lines relative to one another.
- **Profile flattening** (``core.processing.flatten_profile``) acts on the 1-D
  height trace produced after each row group has been averaged.

Both are "flattening" in the usual AFM sense; they differ in the data they
operate on.

Affine image flattening cannot change a blaze angle, and this was measured
directly. On ``20250820_280C_00004.txt``, the scan with the worst row-offset
spread in this project's dataset (2.74 nm), ``none``, ``plane`` and
``align_rows`` all give 30.8524 degrees, identical to four decimal places and
to within 5e-15 in full precision, and ``align_rows`` gives 0.0000 degrees of
change across all eight compare-mode samples.

The reason lies in the structure of the pipeline. ``raw_data_multi_group``
averages rows into a profile, and ``flatten_profile`` then fits and removes a
background from that profile. Subtracting a per-row constant shifts a band's
average by a constant; subtracting a plane adds a constant and a linear ramp.
Every profile-flattening method fits at least a first-order polynomial, so it
removes exactly those terms again, whether or not they were removed here.

This should be kept in mind when benchmarking these methods, since they would
otherwise appear to be broken. Their uses are as follows:

- **Seeing the data.** A raw AFM image with sample tilt renders as a gradient
  with the grating barely visible. Once levelled, it looks like a grating.
- **Methods that are not affine.** A per-row quadratic changes the shape within
  each row in addition to its offset or tilt, and does survive the averaging,
  measured at 0.0003 degrees, which is small but not zero. The registry below
  is the place to add such a method.
- **Being explicit.** The step now happens in a named, recorded place instead
  of in Gwyddion before export, where nothing recorded that it had happened.
"""
from __future__ import annotations

import warnings

import numpy as np

__all__ = ["flatten_image", "IMAGE_FLATTEN_METHODS", "VALID_IMAGE_FLATTEN_METHODS",
           "DEFAULT_IMAGE_FLATTEN_METHOD", "row_offset_spread"]

DEFAULT_IMAGE_FLATTEN_METHOD = "align_rows"


def _none(data: np.ndarray) -> np.ndarray:
    """Leave the image alone, which was the behaviour before 2026-08."""
    return data


def _plane(data: np.ndarray) -> np.ndarray:
    """
    Subtract a single tilted plane fitted to the whole image.

    This corrects for the sample sitting at an angle under the scanner, which
    is the most common background and the one that makes a raw image
    unreadable.

    The plane is fitted by least squares over every pixel, so a deep grating
    biases the fit slightly toward its own mean. That does not matter here (see
    the module docstring), but it would if this were ever the only flattening
    applied.
    """
    rows, cols = data.shape
    y, x = np.mgrid[0:rows, 0:cols]
    basis = np.column_stack([x.ravel(), y.ravel(), np.ones(data.size)])
    # Masked pixels (NaN, see core/region_mask.py) take no part in the fit.
    known = np.isfinite(data.ravel())
    coefficients, *_ = np.linalg.lstsq(basis[known], data.ravel()[known],
                                       rcond=None)
    return data - (basis @ coefficients).reshape(rows, cols)


def _align_rows(data: np.ndarray) -> np.ndarray:
    """
    Put every scan line on a common level by subtracting its median.

    AFM scan lines drift in Z relative to one another (through feedback
    settling and thermal creep), leaving rows offset by a few nanometres. In
    this project's data that spread runs 0.44 to 2.74 nm against groove depths
    near 100 nm.

    The median is used in place of the mean because a grating profile is not
    symmetric about its own mean. A row that happens to catch more groove
    bottom than land would shift a mean, whereas the median barely moves.
    """
    return data - _row_medians(data)


def _row_medians(data):
    """
    Each line's median, over the pixels that are not masked (NaN).

    A masked particle would otherwise move the median of every line it sits
    on. Plain ``np.median`` when nothing is masked, so the result is unchanged
    for every scan analysed before masks existed.
    """
    if np.isfinite(data).all():
        return np.median(data, axis=1, keepdims=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)  # a fully masked line
        return np.nanmedian(data, axis=1, keepdims=True)


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
        A new array of the same shape. The input is never modified in place,
        because callers hold on to the raw image to show a before/after
        comparison.
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
    Standard deviation of the per-row medians, a measure of how far scan lines disagree.

    This is a one-number summary of what `align_rows` corrects, which the
    Import tab displays so the choice of method can be based on evidence. It
    is zero after `align_rows` by construction.
    """
    data = np.asarray(data, dtype=float)
    medians = _row_medians(data).ravel()
    return float(np.std(medians[np.isfinite(medians)]))
