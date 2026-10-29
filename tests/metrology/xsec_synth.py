"""
Synthetic SEM cross-sections of a groove whose shape is known exactly.

An image is rendered from a :class:`~gratinglab.profiles.Trapezoidal` the way
a secondary-electron detector sees a cleave face: vacuum dark, solid grey, a
bright rim just inside every surface where electrons escape, a Gaussian beam,
shot noise, then quantised to the 8-bit levels SmartSEM stores. Pixel coverage
is computed on a supersampled grid, so a pixel straddling the edge takes the
fraction of its area that is solid and the true boundary is known to far better
than a pixel.
"""
import numpy as np
from scipy import ndimage

from gratinglab.profiles import Trapezoidal


def render(
    profile: Trapezoidal,
    *,
    period_nm: float = 163.3,
    pixel_nm: float = 2.233,
    width_px: int = 640,
    height_px: int = 160,
    surface_row: float = 50.0,
    phase_nm: float = 37.0,
    vacuum: float = 40.0,
    material: float = 110.0,
    rim: float = 0.0,
    rim_nm: float = 4.0,
    beam_px: float = 0.7,
    noise: float = 0.0,
    slope: float = 0.0,
    height_scale: float = 1.0,
    seed: int = 0,
    supersample: int = 8,
) -> np.ndarray:
    """An 8-bit cross-section, material below, crest at ``surface_row``.

    ``slope`` tilts the cleave edge (rise per unit length), ``height_scale``
    compresses every height as a tilted stage would.
    """
    ss = supersample
    sub = pixel_nm / ss
    xs = (np.arange(width_px * ss) + 0.5) * sub
    zs = (np.arange(height_px * ss) + 0.5) * sub  # downward from the top edge
    t = np.mod((xs - phase_nm) / period_nm, 1.0)
    crest_depth = profile.depth * period_nm
    # Height above the floor, then the depth below the top edge of the frame.
    surface = (
        surface_row * pixel_nm
        + (crest_depth - profile.height(t) * period_nm) * height_scale
        - slope * xs
    )
    solid = zs[:, None] >= surface[None, :]

    value = np.where(solid, material, vacuum).astype(np.float64)
    if rim:
        inside = ndimage.distance_transform_edt(solid) * sub
        value += np.where(solid, rim * np.exp(-inside / rim_nm), 0.0)

    pixels = value.reshape(height_px, ss, width_px, ss).mean(axis=(1, 3))
    if beam_px:
        pixels = ndimage.gaussian_filter(pixels, beam_px, mode="nearest")
    if noise:
        pixels = pixels + np.random.default_rng(seed).normal(0.0, noise, pixels.shape)
    return np.clip(np.round(pixels), 0, 254).astype(np.uint8)
