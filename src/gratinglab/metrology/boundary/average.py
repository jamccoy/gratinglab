"""
Groove averaging and normalisation for PCGrate boundary profiles

Ported from afm_scan_avg_profile_rev3.py. The arithmetic is unchanged - the
port is verified by reproducing that script's output byte-for-byte - but the
profile front-end now comes from afm_analysis.core instead of a private copy.
"""
import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import uniform_filter1d

from ..core.processing import extract_single_groove


def flatten_endpoints(x, y):
    """
    Flatten a curve so it starts and ends at zero.

    This is deliberately not the same as core.processing.flatten_profile. That
    one removes a fitted background to measure facet angles accurately; this one
    forces the endpoints to zero so a single groove tiles seamlessly into a
    periodic boundary. Different jobs, different maths - they are not
    interchangeable.
    """
    x0, x1 = x[0], x[-1]
    y0, y1 = y[0], y[-1]
    slope = (y1 - y0) / (x1 - x0)
    line = y0 + slope * (x - x0)
    return y - line


def average_grooves(x, y, groove_centers, period_nm, margin=0.0, n_points=2000,
                    min_half_width=10):
    """
    Extract each groove symmetrically and average them onto a common axis.

    Parameters:
        min_half_width: skip grooves whose symmetric half-extent is at or below
            this many samples. A groove near a scan edge cannot be extracted
            symmetrically - the scan stops part-way through it - so including it
            would distort the average. This is the boundary-profile equivalent of
            EDGE_EXCLUSION_PERIODS in the blaze-angle path.

    Returns:
        x_common, y_avg, y_std, n_used
    """
    grooves = []
    dx_nm = (x[1] - x[0]) * 1000
    half_width = int(round((period_nm / 2) * (1 + margin) / dx_nm))

    for center in groove_centers:
        L = min(half_width, center, len(x) - 1 - center)
        if L > min_half_width:
            x_g, y_g = extract_single_groove(x, y, center, period_nm, margin,
                                             allow_asymmetric=False)
            # Re-create the exactly symmetric axis the original used, so that
            # floating point matches to the last bit.
            x_g = np.linspace(-L * dx_nm, L * dx_nm, len(y_g)) / 1000
            grooves.append((x_g, y_g))

    if len(grooves) == 0:
        raise ValueError("No grooves could be extracted.")

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


def normalize_profile(x_avg, y_avg, period_nm, apply_smoothing=True,
                      smoothing_window=5):
    """
    Turn an averaged groove into a period-normalised, seamlessly tiling profile.

    Four steps, in the order the original performed them:
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
