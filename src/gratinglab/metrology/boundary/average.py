"""
Groove averaging and normalisation for PCGrate boundary profiles

Ported from afm_scan_avg_profile_rev3.py. The arithmetic is unchanged, and the
port is verified by reproducing that script's output byte-for-byte, but the
profile front-end now comes from gratinglab.metrology.core instead of a private
copy.
"""
import numpy as np
from scipy.interpolate import interp1d, make_interp_spline
from scipy.ndimage import uniform_filter1d

from ..core.processing import extract_single_groove


def flatten_endpoints(x, y):
    """
    Flatten a curve so it starts and ends at zero.

    This is not the same as core.processing.flatten_profile. That
    one removes a fitted background to measure facet angles accurately; this one
    forces the endpoints to zero so a single groove tiles seamlessly into a
    periodic boundary. The two serve different purposes with different
    mathematics and are not interchangeable.
    """
    x0, x1 = x[0], x[-1]
    y0, y1 = y[0], y[-1]
    slope = (y1 - y0) / (x1 - x0)
    line = y0 + slope * (x - x0)
    return y - line


def select_grooves(x, y, groove_centers, period_nm, margin=0.0,
                   min_half_width=10, require_full_period=True):
    """
    The grooves `average_grooves` averages, each centred on its own axis in µm.

    Separate so that `groove_residuals` measures exactly the set that went
    into the average, chosen by the same rules rather than a second copy of
    them. The parameters mean what they mean in `average_grooves`.

    Returns:
        a list of (x_g, y_g), possibly empty
    """
    grooves = []
    dx_nm = (x[1] - x[0]) * 1000
    half_width = int(round((period_nm / 2) * (1 + margin) / dx_nm))

    for center in groove_centers:
        L = min(half_width, center, len(x) - 1 - center)
        if require_full_period and L < half_width:
            continue  # the scan edge cuts into this groove's period
        if L > min_half_width:
            x_g, y_g = extract_single_groove(x, y, center, period_nm, margin,
                                             allow_asymmetric=False)
            # Re-create the exactly symmetric axis the original used, so that
            # floating point matches to the last bit.
            x_g = np.linspace(-L * dx_nm, L * dx_nm, len(y_g)) / 1000
            grooves.append((x_g, y_g))
    return grooves


def average_grooves(x, y, groove_centers, period_nm, margin=0.0, n_points=2000,
                    min_half_width=10, require_full_period=True):
    """
    Extract each groove symmetrically and average them onto a common axis.

    Parameters:
        min_half_width: skip grooves whose symmetric half-extent is at or below
            this many samples, a floor on how coarsely a groove may be sampled.
        require_full_period: skip any groove whose full one-period window does
            not fit inside the scan. The common axis below is the narrowest
            window in the set, and normalize_profile stretches that window to a
            full period, so a single groove the scan edge cuts into would
            truncate every other groove to its window and then stretch the
            result. On the synthetic sawtooth one partial groove narrowed the
            window to 0.66 of a period, turned a 70 degree anti-blaze facet
            into 58 and took 13% off the depth. EDGE_EXCLUSION_PERIODS upstream
            usually keeps such grooves out, but it is a user setting, so this
            check does not rely on it. False reproduces the retired rev3
            script, which tests/metrology/test_ggp_equivalence.py pins.

    Returns:
        x_common, y_avg, y_std, n_used
    """
    grooves = select_grooves(x, y, groove_centers, period_nm, margin=margin,
                             min_half_width=min_half_width,
                             require_full_period=require_full_period)

    if len(grooves) == 0:
        raise ValueError(
            "No grooves could be extracted: none has a full period inside the "
            "scan, or the sampling is coarser than min_half_width allows.")

    # Common symmetric extent = most restrictive of the set
    max_left = min(-g[0][0] for g in grooves)
    max_right = min(g[0][-1] for g in grooves)
    half_range = min(max_left, max_right)

    x_common = np.linspace(-half_range, half_range, n_points)

    y_aligned = []
    for g in grooves:
        interp_func = interp1d(g[0], g[1], kind='cubic', fill_value='extrapolate')
        y_aligned.append(interp_func(x_common))
    y_aligned = np.array(y_aligned)

    return x_common, np.mean(y_aligned, axis=0), np.std(y_aligned, axis=0), len(grooves)


def groove_residuals(grooves, x_common, y_avg, n_passes=3, n_steps=6):
    """
    How far each groove departs from the average, once registered to it.

    This is the RMS residual of each averaged groove about the average: the
    part of the measured profile that does not repeat from one groove to the
    next. It is reported, not applied. `BoundaryProfile.to_problem` and
    docs/theory/metrology.md section 4 say why it is not a Nevot-Croce sigma.

    A plain difference against `y_avg` does not measure that. Each groove is
    centred on the sample nearest its minimum, so the set is misregistered by
    up to half a pixel, and on a facet steeper than 45 degrees half a pixel of
    shift is several nanometres of height. On the synthetic sawtooth with no
    noise at all, the plain difference reads 2.1 nm. So each groove is first
    registered to the average by a sub-pixel shift and a height offset, fitted
    by Gauss-Newton on a cubic spline of its own samples, and the average is
    recomputed from the registered grooves, `n_passes` times over. The shift is
    bounded to one pixel either way, and the window is the aligned window the
    average uses less one pixel at each end, so a shifted groove is never
    extrapolated.

    The residual is then taken at each groove's own samples rather than on the
    dense common axis, because interpolating onto that axis smooths pixel-scale
    structure: it read white noise 10 to 15% low. It is scaled by
    sqrt(N / (N - 1)), since the average contains the groove it is compared
    with, and is undefined for a single groove.

    `y_avg` itself is not changed. It is still the unregistered mean, and is
    what is exported.

    Returns:
        dict with
          rms_nm: over every groove, or None for fewer than two grooves
          per_groove_nm: one RMS per groove, in the order of `grooves`
          shift_px: each groove's registration shift, in pixels
    """
    n = len(grooves)
    if n < 2:
        return {'rms_nm': None, 'per_groove_nm': np.zeros(0),
                'shift_px': np.zeros(0)}

    dx = grooves[0][0][1] - grooves[0][0][0]
    inner = (x_common > x_common[0] + dx) & (x_common < x_common[-1] - dx)
    xi = x_common[inner]
    reference = np.asarray(y_avg, dtype=float)[inner]
    splines = [make_interp_spline(g[0], g[1], k=3) for g in grooves]
    slopes = [spline.derivative() for spline in splines]
    shifts = np.zeros(n)

    for _ in range(n_passes):
        registered = np.empty((n, len(xi)))
        for k, (spline, slope) in enumerate(zip(splines, slopes)):
            d = shifts[k]
            for _ in range(n_steps):
                r = spline(xi + d) - reference
                # Least squares in (shift, offset). The offset's normal
                # equation removes the means, which leaves one unknown.
                g = slope(xi + d)
                g = g - g.mean()
                denom = np.dot(g, g)
                if denom == 0:
                    break
                d = float(np.clip(d - np.dot(g, r - r.mean()) / denom, -dx, dx))
            shifts[k] = d
            curve = spline(xi + d)
            registered[k] = curve - np.mean(curve - reference)
        reference = registered.mean(axis=0)

    at_reference = make_interp_spline(xi, reference, k=3)
    mean_square = np.empty(n)
    for k, (x_g, y_g) in enumerate(grooves):
        xs = x_g - shifts[k]
        keep = (xs >= xi[0]) & (xs <= xi[-1])
        r = y_g[keep] - at_reference(xs[keep])
        mean_square[k] = np.mean((r - r.mean()) ** 2)
    mean_square *= n / (n - 1)

    return {'rms_nm': float(np.sqrt(np.mean(mean_square))),
            'per_groove_nm': np.sqrt(mean_square),
            'shift_px': shifts / dx}


def normalize_profile(x_avg, y_avg, period_nm, apply_smoothing=True,
                      smoothing_window=5):
    """
    Turn an averaged groove into a period-normalised, seamlessly tiling profile.

    The four steps below are performed in the order the original used:
      1. flatten endpoints, then shift so the groove edges sit at zero
      2. scale x to [0, 1] over one period; express y as a fraction of the period
      3. roll the array so the minimum lands on the boundary, and force the
         endpoints to exactly zero so successive periods join without a step
      4. optionally smooth (wrapping at the boundary), then re-force the endpoints

    Returns:
        x_normalized, y_normalized, edge_height_nm
    """
    y_avg_flat = flatten_endpoints(x_avg, y_avg)

    # Shift so groove edges are at zero
    mid_idx = len(y_avg_flat) // 2
    left_min = np.min(y_avg_flat[:mid_idx])
    right_min = np.min(y_avg_flat[mid_idx:])
    edge_height = (left_min + right_min) / 2
    y_avg_shifted = y_avg_flat - edge_height

    # Normalise: x to [0, 1] over one period, y as a fraction of the period
    x_normalized = (x_avg - x_avg[0]) / (x_avg[-1] - x_avg[0])
    y_normalized = y_avg_shifted / period_nm

    # Minimum to zero
    y_min_at_zero = y_normalized - np.min(y_normalized)

    # Phase shift so the minimum sits on the period boundary
    imin = np.argmin(y_min_at_zero)
    y_phase_shifted = np.roll(y_min_at_zero, -imin)
    y_phase_shifted[0] = 0.0
    y_phase_shifted[-1] = 0.0

    if apply_smoothing:
        y_phase_shifted = uniform_filter1d(y_phase_shifted,
                                           size=smoothing_window, mode='wrap')
        y_phase_shifted[0] = 0.0
        y_phase_shifted[-1] = 0.0

    return x_normalized, y_phase_shifted, edge_height


def profile_metrics(x_normalized, y_normalized, period_nm, n_grooves):
    """Slope and curvature statistics of the normalised profile"""
    dx_norm = x_normalized[1] - x_normalized[0]
    dy_dx = np.gradient(y_normalized, dx_norm)
    d2y_dx2 = np.gradient(dy_dx, dx_norm)

    max_slope = np.max(np.abs(dy_dx))

    return {
        'period_nm': period_nm,
        'n_grooves': n_grooves,
        'groove_depth': np.max(y_normalized),
        'peak_to_valley': np.max(y_normalized) - np.min(y_normalized),
        'rms_slope': np.sqrt(np.mean(dy_dx ** 2)),
        'max_slope': max_slope,
        'max_angle_deg': np.arctan(max_slope) * 180 / np.pi,
        'max_curvature': np.max(np.abs(d2y_dx2)),
    }
