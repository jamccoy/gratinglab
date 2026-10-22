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

The tip is a spherical apex cap of ``radius_nm`` on two straight flanks, one
either side of the apex in the scan direction. Each flank has its own half
angle, **measured from the tip axis**, so a flank of half angle ``a`` rises at
``90 - a`` degrees from the surface. ``half_angle_deg`` is the flank on the
-x side (towards column 0) and ``half_angle_right_deg`` the one on the +x
side; leaving the second out gives a symmetric cone, which is what every
earlier version of this module assumed.

The asymmetry is not a refinement for its own sake. A real tip is a pyramid,
and which of its faces meet the surface depends on how the cantilever lies
relative to the fast-scan axis. Bruker's ScanAsyst-Air, for instance, has a
15 deg front face, a 25 deg back face and 17.5 deg sides (each +/- 2.5 deg).
Scanned across the cantilever the profile sees the two sides, a symmetric 72.5
deg flank either way. Scanned along it, one flank rises at 75 deg and the other
at only **65**, which cannot follow a 70 deg anti-blaze facet at all. See
:data:`TIP_PRESETS`.

The consequence worth knowing before trusting any corrected groove: a facet
steeper than the flank facing it is touched only at its top corner and is
**unrecoverable**, and a facet close to the flank angle is the marginal case.
That near-parallelism, not the nanometre apex, is what rounds the troughs of
blazed gratings.

One dimension only. The operations run along the fast-scan axis (rows), not
over the full 2-D image. For a grating the surface is invariant along the
grooves, and dilating an invariant surface with the tip reduces exactly to
dilating each cross-section with the tip's cross-section in the scan plane --
the groove-parallel extent of the tip touches nothing the transverse extent
does not. What a 2-D pass would add is sensitivity to along-groove features,
and it would need the slow-axis pixel pitch, which the text-export format does
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


@dataclass(frozen=True)
class TipPreset:
    """A named tip, as the three numbers the model takes."""

    label: str
    radius_nm: float
    half_angle_deg: float          #: the -x flank
    half_angle_right_deg: float    #: the +x flank
    note: str


#: Probes worth having by name. Bruker's figures for the ScanAsyst-Air
#: (brukerafmprobes.com, p-3726): tip radius 2 nm nominal and 12 nm maximum,
#: front angle 15, back angle 25 and side angle 17.5 deg, each +/- 2.5 deg and
#: each a half angle from the tip axis. The radius is a radius as quoted, not a
#: width. A worn tip is blunter than nominal, so a correction made with the
#: nominal radius is the optimistic one.
#:
#: Which face looks which way along the scan is a fact about the instrument and
#: the scan angle, not about the tip, so the along-cantilever case comes in both
#: orientations and choosing between them is the user's.
TIP_PRESETS = {
    'scanasyst-air-across': TipPreset(
        'ScanAsyst-Air, scanning across the cantilever', 2.0, 17.5, 17.5,
        "the two 17.5 deg side faces: a 72.5 deg flank either way"),
    'scanasyst-air-along-front-minus': TipPreset(
        'ScanAsyst-Air, along the cantilever, front face to -x', 2.0, 15.0, 25.0,
        "front face (15 deg, a 75 deg flank) to -x, back face (25 deg, a 65 "
        "deg flank) to +x"),
    'scanasyst-air-along-front-plus': TipPreset(
        'ScanAsyst-Air, along the cantilever, front face to +x', 2.0, 25.0, 15.0,
        "back face (25 deg, a 65 deg flank) to -x, front face (15 deg, a 75 "
        "deg flank) to +x"),
}


def format_angles(half_angle_deg, half_angle_right_deg=None) -> str:
    """The flank angles the way a reader wants them: one figure when they agree."""
    if half_angle_right_deg is None or half_angle_right_deg == half_angle_deg:
        return f"half angle = {half_angle_deg:g} deg"
    return (f"half angles = {half_angle_deg:g} deg (-x) / "
            f"{half_angle_right_deg:g} deg (+x)")


def _check_angle(name, value):
    if not 0.0 < value < 90.0:
        raise ValueError(f"tip {name} must lie in (0, 90) degrees, got {value}")


def _flank_reach(radius_nm, theta, height_nm):
    """Lateral offset at which one side of the tip reaches ``height_nm``."""
    t_tangent = radius_nm * (1.0 - np.sin(theta))
    if height_nm <= t_tangent:
        return float(np.sqrt(max(radius_nm**2 - (radius_nm - height_nm)**2, 0.0)))
    return float(radius_nm * np.cos(theta) + (height_nm - t_tangent) * np.tan(theta))


def _flank_profile(r, radius_nm, theta):
    """Height above the apex at lateral distance ``r >= 0`` on one side."""
    r_tangent = radius_nm * np.cos(theta)
    t_tangent = radius_nm * (1.0 - np.sin(theta))
    sphere = radius_nm - np.sqrt(np.maximum(radius_nm**2 - r**2, 0.0))
    cone = t_tangent + (r - r_tangent) / np.tan(theta)
    return np.where(r <= r_tangent, sphere, cone)


def tip_cross_section(dx_nm, radius_nm, half_angle_deg, max_height_nm,
                      half_angle_right_deg=None):
    """The tip's profile in the scan plane, sampled at the pixel pitch.

    Returns heights in nm above the apex at lateral offsets ``k * dx_nm`` for
    ``k in [-K, K]`` (odd length, apex at the centre, ``t[K] == 0``), negative
    offsets on the ``half_angle_deg`` flank and positive ones on
    ``half_angle_right_deg``, which defaults to the same. ``K`` is chosen so
    the profile just exceeds ``max_height_nm`` on its wider side -- parts of
    the tip higher than the tallest feature can never make contact and would
    only widen every window.

    The sphere of radius ``R`` meets each flank where its tangent matches it,
    at lateral offset ``R cos(theta)`` and height ``R (1 - sin(theta))`` for
    that side's ``theta``, so the slope is continuous on both sides however
    different they are.
    """
    if half_angle_right_deg is None:
        half_angle_right_deg = half_angle_deg
    if radius_nm <= 0:
        raise ValueError(f"tip radius must be positive, got {radius_nm} nm")
    _check_angle("half angle", half_angle_deg)
    _check_angle("right half angle", half_angle_right_deg)

    left = np.radians(half_angle_deg)
    right = np.radians(half_angle_right_deg)
    r_max = max(_flank_reach(radius_nm, left, max_height_nm),
                _flank_reach(radius_nm, right, max_height_nm))

    # At least the nearest neighbours, so a sub-pixel tip degrades to a
    # near-no-op rather than an empty window.
    k = max(int(np.ceil(r_max / dx_nm)), 1)
    offsets = np.arange(-k, k + 1) * dx_nm
    r = np.abs(offsets)
    return np.where(offsets < 0,
                    _flank_profile(r, radius_nm, left),
                    _flank_profile(r, radius_nm, right))


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


def dilate(data, scan_x_size, *, radius_nm, half_angle_deg,
           half_angle_right_deg=None):
    """Image formation: what this tip reports for the surface ``data``.

    ``img(x) = max_u [ surface(x + u) - tip(u) ]`` -- the apex height when the
    tip, centred at ``x``, rests on whichever point it touches first; ``u`` is
    the offset from the apex, positive to +x. Always ``>= data`` pointwise: a
    tip can hide a trough but never dig one.

    Exists for synthesis and testing; the pipeline itself only erodes.
    """
    dx_nm = _pixel_pitch_nm(data, scan_x_size)
    relief_nm = float(data.max() - data.min()) * 1e9
    tip_m = tip_cross_section(dx_nm, radius_nm, half_angle_deg, relief_nm,
                              half_angle_right_deg) * 1e-9
    k_max = len(tip_m) // 2

    img = np.full_like(data, -np.inf)
    for k in range(-k_max, k_max + 1):
        np.maximum(img, _shifted(data, k, -np.inf) - tip_m[k + k_max], out=img)
    return img


@dataclass(frozen=True)
class TipCorrection:
    """An eroded image, and how much of it is the surface.

    ``data`` is the least upper bound on the true surface: equal to it at
    every certain pixel, above it (by an unknowable amount) elsewhere. The
    parameters are carried so downstream records can state what was assumed
    rather than that something was.
    """

    data: np.ndarray       #: corrected heights, metres, same shape as input
    certain: np.ndarray    #: bool; True where the apex provably made contact
    radius_nm: float
    half_angle_deg: float  #: the -x flank
    footprint_px: int      #: window width used; 3 means the tip is sub-pixel
    half_angle_right_deg: float = None  #: the +x flank; None means symmetric

    @property
    def symmetric(self) -> bool:
        return (self.half_angle_right_deg is None
                or self.half_angle_right_deg == self.half_angle_deg)

    @property
    def angles(self) -> str:
        """The flank angles, the way a reader wants them written."""
        return format_angles(self.half_angle_deg, self.half_angle_right_deg)

    @property
    def certain_fraction(self) -> float:
        """Fraction of pixels where the correction recovered the surface."""
        return float(np.mean(self.certain))

    @property
    def summary(self) -> str:
        """One line for a log or a metrics sidecar."""
        note = (" -- tip is sub-pixel at this pitch, correction is a near-no-op"
                if self.footprint_px <= 3 else "")
        return (f"erosion, R = {self.radius_nm:g} nm, {self.angles}: "
                f"{100.0 * self.certain_fraction:.1f}% of pixels certain, "
                f"the rest upper bounds{note}")


def erode(data, scan_x_size, *, radius_nm, half_angle_deg,
          half_angle_right_deg=None) -> TipCorrection:
    """Reconstruct the least upper bound on the surface behind an image.

    ``rec(y) = min_u [ img(y - u) + tip(u) ]`` -- the lowest the surface can
    be at ``y`` given that the tip, wherever it stood, stopped where the image
    says: with the apex at ``y - u`` the point ``y`` lies under the tip at
    offset ``u``. Sits between the truth and the image:
    ``surface <= rec <= img``.

    Note the ``y - u``. Written as a loop over shifts ``k = -u`` it pairs
    ``img(y + k)`` with ``tip(-k)``, the tip **reflected**. For a symmetric
    tip the reflection is invisible, which is why it matters to say so: with
    an asymmetric one, pairing ``img(y + k)`` with ``tip(k)`` erodes by the
    mirror image of the tip. On the asymmetric test wall that version is no
    longer the least upper bound (4.09 nm above the surface on average against
    3.53) and re-imaging it misses the image by up to 2.7 nm, where the right
    one reproduces it to 1e-14. ``test_tip.py`` catches it three ways.

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
    tip_m = tip_cross_section(dx_nm, radius_nm, half_angle_deg, relief_nm,
                              half_angle_right_deg) * 1e-9
    k_max = len(tip_m) // 2

    rec = np.full_like(data, np.inf)
    for k in range(-k_max, k_max + 1):
        np.minimum(rec, _shifted(data, k, np.inf) + tip_m[k_max - k], out=rec)

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
        half_angle_right_deg=half_angle_right_deg,
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
                       half_angle_deg=settings.tip_half_angle_deg,
                       half_angle_right_deg=settings.tip_half_angle_right_deg)
    return correction.data, correction
