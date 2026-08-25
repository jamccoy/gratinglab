"""
Tip convolution: simulate it, and undo what can be undone.

An AFM image is not the surface. It is the surface **dilated** by the tip
shape -- at every pixel the tip descends until it touches, and the recorded
height is where the *apex* stopped, which is above the true surface anywhere
the contact happened on the tip's flank instead of its apex. Villarrubia
(1997, J. Res. NIST 102, 425) gives the exact morphological statement and the
two operations this module implements:

- ``dilate`` -- image formation. What a tip of the given shape would report
  for a known surface. Exact, and used to *make* test data: dilating the
  synthetic scan produces an image whose true surface is known, which no real
  scan allows.
- ``erode`` -- reconstruction. The **least upper bound** on the true surface
  consistent with the image and the tip: everywhere the apex made contact the
  bound *is* the surface, and everywhere else the tip physically could not
  reach, so no algorithm can recover what is there. The certainty map marks
  the difference, which is the entire point: a reconstruction without one
  manufactures exactly the false confidence it exists to remove.

The tip is a cone with a spherical apex cap -- ``radius_nm`` and
``half_angle_deg`` (measured from the tip axis, so the flank rises at
``90 - half_angle`` degrees from the surface). The consequence worth knowing
before trusting any corrected groove: a facet steeper than the flank is
touched only at its top corner and is **unrecoverable**, and a facet close to
the flank angle is the marginal case -- an 18 deg tip has a 72 deg flank,
almost parallel to a 70.5 deg anti-blaze facet. That near-parallelism, not
the nanometre apex, is what rounds the troughs of blazed gratings.

**One dimension, deliberately.** The operations run along the fast-scan axis
(rows), not over the full 2-D image. For a grating the surface is invariant
along the grooves, and dilating an invariant surface with a cone reduces
exactly to dilating each cross-section with the cone's cross-section -- the
groove-parallel extent of the tip touches nothing the transverse extent does
not. What a 2-D pass would add is sensitivity to along-groove features, and
it would need the slow-axis pixel pitch, which the text-export format does
not record.

All heights are in metres (the unit the loaders promise), tip parameters in
nm, and the pixel pitch comes from the scan width exactly as
``core.processing.raw_data`` computes it.
"""
from dataclasses import dataclass

import numpy as np

#: Values ``AnalysisSettings.tip_correction`` may take. 'none' is the default:
#: the correction changes measured numbers and must be asked for, not applied
#: silently.
VALID_TIP_CORRECTIONS = ('none', 'erosion')

#: Height differences below this (metres) are numerically meaningless: 1e-12 m
#: is a thousandth of a nanometre, far under any AFM noise floor. Used to call
#: two candidate contacts a tie, which marks the pixel uncertain.
_TIE_TOL_M = 1e-12


def tip_cross_section(dx_nm, radius_nm, half_angle_deg, max_height_nm):
    """The tip's transverse profile, sampled at the pixel pitch.

    Returns heights in nm above the apex at lateral offsets
    ``k * dx_nm`` for ``k in [-K, K]`` (odd length, apex at the centre,
    ``t[K] == 0``). ``K`` is chosen so the profile just exceeds
    ``max_height_nm`` -- parts of the tip higher than the tallest feature can
    never make contact and would only widen every window.

    Sphere of radius ``R`` capping a cone of half-angle ``theta``: the two
    meet where the sphere's tangent matches the cone flank, at lateral offset
    ``R cos(theta)`` and height ``R (1 - sin(theta))``.
    """
    if radius_nm <= 0:
        raise ValueError(f"tip radius must be positive, got {radius_nm} nm")
    if not 0.0 < half_angle_deg < 90.0:
        raise ValueError(
            f"tip half angle must lie in (0, 90) degrees, got {half_angle_deg}")

    theta = np.radians(half_angle_deg)
    r_tangent = radius_nm * np.cos(theta)
    t_tangent = radius_nm * (1.0 - np.sin(theta))

    # Lateral extent at which the tip profile reaches max_height_nm.
    if max_height_nm <= t_tangent:
        r_max = float(np.sqrt(max(radius_nm**2 - (radius_nm - max_height_nm)**2,
                                  0.0)))
    else:
        r_max = float(r_tangent + (max_height_nm - t_tangent) * np.tan(theta))

    # At least the nearest neighbours, so a sub-pixel tip degrades to a
    # near-no-op rather than an empty window.
    k = max(int(np.ceil(r_max / dx_nm)), 1)
    r = np.abs(np.arange(-k, k + 1)) * dx_nm

    sphere = radius_nm - np.sqrt(np.maximum(radius_nm**2 - r**2, 0.0))
    cone = t_tangent + (r - r_tangent) / np.tan(theta)
    return np.where(r <= r_tangent, sphere, cone)


def _pixel_pitch_nm(data, scan_x_size):
    """nm per column, matching ``raw_data``'s displacement axis exactly."""
    return scan_x_size * 1000.0 / (data.shape[1] - 1)


def _shifted(data, k, fill):
    """``out[:, j] = data[:, j + k]``, with ``fill`` where that runs off."""
    out = np.full_like(data, fill)
    if k == 0:
        return data.copy()
    if k > 0:
        out[:, :-k] = data[:, k:]
    else:
        out[:, -k:] = data[:, :k]
    return out


def dilate(data, scan_x_size, *, radius_nm, half_angle_deg):
    """Image formation: what this tip reports for the surface ``data``.

    ``img(x) = max_u [ surface(x + u) - tip(u) ]`` -- the apex height when the
    tip, centred at ``x``, rests on whichever point it touches first. Always
    ``>= data`` pointwise: a tip can hide a trough but never dig one.

    Exists for synthesis and testing; the pipeline itself only erodes.
    """
    dx_nm = _pixel_pitch_nm(data, scan_x_size)
    relief_nm = float(data.max() - data.min()) * 1e9
    tip_m = tip_cross_section(dx_nm, radius_nm, half_angle_deg, relief_nm) * 1e-9
    k_max = len(tip_m) // 2

    img = np.full_like(data, -np.inf)
    for k in range(-k_max, k_max + 1):
        np.maximum(img, _shifted(data, k, -np.inf) - tip_m[k + k_max], out=img)
    return img


@dataclass(frozen=True)
class TipCorrection:
    """An eroded image, and how much of it is actually the surface.

    ``data`` is the least upper bound on the true surface: equal to it at
    every certain pixel, above it (by an unknowable amount) elsewhere. The
    parameters are carried so downstream records can state what was assumed
    rather than that something was.
    """

    data: np.ndarray       #: corrected heights, metres, same shape as input
    certain: np.ndarray    #: bool; True where the apex provably made contact
    radius_nm: float
    half_angle_deg: float
    footprint_px: int      #: window width used; 3 means the tip is sub-pixel

    @property
    def certain_fraction(self) -> float:
        """Fraction of pixels where the correction recovered the surface."""
        return float(np.mean(self.certain))

    @property
    def summary(self) -> str:
        """One line for a log or a metrics sidecar."""
        note = (" -- tip is sub-pixel at this pitch, correction is a near-no-op"
                if self.footprint_px <= 3 else "")
        return (f"erosion, R = {self.radius_nm:g} nm, "
                f"half angle = {self.half_angle_deg:g} deg: "
                f"{100.0 * self.certain_fraction:.1f}% of pixels certain, "
                f"the rest upper bounds{note}")


def erode(data, scan_x_size, *, radius_nm, half_angle_deg) -> TipCorrection:
    """Reconstruct the least upper bound on the surface behind an image.

    ``rec(x) = min_u [ img(x + u) + tip(u) ]`` -- the lowest the surface can
    be at ``x`` given that the tip, wherever it stood, stopped where the image
    says. Sits between the truth and the image: ``surface <= rec <= img``.

    The certainty map re-runs the contact argument: re-dilating ``rec``
    reproduces the image exactly, and for each image pixel the surface point
    achieving that maximum is where the tip touched. A pixel is certain when
    it is such a contact point for some image pixel **and** the contact is
    unique -- a tie means the contact cannot be localised (the flank lying
    flush along a facet, a wedge gripped at both walls), and every tied point
    stays an upper bound. Window-truncated pixels at the scan edges are never
    marked certain.
    """
    dx_nm = _pixel_pitch_nm(data, scan_x_size)
    relief_nm = float(data.max() - data.min()) * 1e9
    tip_m = tip_cross_section(dx_nm, radius_nm, half_angle_deg, relief_nm) * 1e-9
    k_max = len(tip_m) // 2

    rec = np.full_like(data, np.inf)
    for k in range(-k_max, k_max + 1):
        np.minimum(rec, _shifted(data, k, np.inf) + tip_m[k + k_max], out=rec)

    # Contact scan: for each image pixel, the best (highest) candidate
    # rec(x+u) - tip(u) is the redilated image; its argmax u is the contact.
    best = np.full_like(data, -np.inf)
    second = np.full_like(data, -np.inf)
    best_k = np.zeros(data.shape, dtype=np.int64)
    for k in range(-k_max, k_max + 1):
        candidate = _shifted(rec, k, -np.inf) - tip_m[k + k_max]
        improves = candidate > best
        second = np.where(improves, best, np.maximum(second, candidate))
        best = np.where(improves, candidate, best)
        best_k = np.where(improves, k, best_k)

    unique = np.isfinite(best) & (best - second > _TIE_TOL_M)
    certain = np.zeros(data.shape, dtype=bool)
    rows, cols = np.nonzero(unique)
    certain[rows, cols + best_k[rows, cols]] = True

    return TipCorrection(
        data=rec,
        certain=certain,
        radius_nm=radius_nm,
        half_angle_deg=half_angle_deg,
        footprint_px=2 * k_max + 1,
    )


def apply_tip_correction(data, scan_x_size, settings):
    """The pipeline entry point: dispatch on ``settings.tip_correction``.

    Returns ``(data, correction)`` -- the array to analyse and the
    :class:`TipCorrection` behind it, or ``(data, None)`` untouched when the
    setting is ``'none'``. Callers record ``correction.summary`` wherever
    their output lands, because a corrected depth and an uncorrected one are
    different measurements and must not be filed under the same description.
    """
    if settings.tip_correction == 'none':
        return data, None
    correction = erode(data, scan_x_size,
                       radius_nm=settings.tip_radius_nm,
                       half_angle_deg=settings.tip_half_angle_deg)
    return correction.data, correction
