"""
Core data processing functions for AFM analysis
"""
import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.signal import find_peaks


def load_afm_data(filename, default_scan_size=2.0, settings=None):
    """
    Load AFM data, whatever format it arrived in.

    The format is detected from the content, so the caller need not know it.
    Raw Nanoscope files (.spm, and the extensionless companions Nanoscope also
    writes) are read directly; everything else takes the text path used for
    Gwyddion exports.

    The dispatch happens here instead of at each call site by design: there
    are six callers, and none of them should need to handle the format.

    Parameters:
        settings: optional AnalysisSettings, consulted for which .spm channel and
            scan direction to read. Ignored for text files.

    Returns:
        data: 2D array of height data, in metres
        scan_x_size: scan width in microns
    """
    from .. import io as _io  # noqa: F401 - keeps the package import explicit
    from ..io.ibw import is_igor_wave_file, read_ibw
    from ..io.spm import (DEFAULT_CHANNEL, DEFAULT_DIRECTION, is_nanoscope_file,
                          read_spm)

    if is_nanoscope_file(filename):
        channel = getattr(settings, 'spm_channel', None) or DEFAULT_CHANNEL
        direction = getattr(settings, 'spm_direction', None) or DEFAULT_DIRECTION
        data, scan_x_size = read_spm(filename, channel=channel,
                                     direction=direction)
        print(f"  Nanoscope file: {channel} / {direction}, "
              f"{data.shape[0]}x{data.shape[1]} px, {scan_x_size:g} µm")
        return data, scan_x_size

    if is_igor_wave_file(filename):
        # The scan direction is shared with .spm; the channel is not, since
        # Asylum names its channels differently ("ZSensor", "Height"), so the
        # reader picks the topography channel itself (io/ibw.py).
        direction = getattr(settings, 'spm_direction', None) or DEFAULT_DIRECTION
        data, scan_x_size = read_ibw(filename, direction=direction)
        print(f"  Asylum file: {direction}, "
              f"{data.shape[0]}x{data.shape[1]} px, {scan_x_size:g} µm")
        return data, scan_x_size

    return _load_text(filename, default_scan_size)


def load_scan_y_size(filename, settings=None):
    """
    The slow-axis length of the recorded lines, in microns, or None.

    Separate from `load_afm_data` so that function's two-value return, which
    six callers unpack, stays as it is. Only groove-orientation measurement
    needs this number, because it is the one step that compares a length
    along x with a length along y.

    None means the file does not say, and the caller must not substitute the
    fast-axis width: a scan of 2 µm by 1 µm read as square would report half
    its real groove rotation, with no sign of the assumption.

    Sources:
        .spm: the second ``Scan Size`` value, scaled to the lines recorded
        .ibw: the Y scaling, or the note's ScanSize, SlowRatio and ScanLines
        text: a ``Height:`` line in the header, as Gwyddion writes it
    """
    from ..io import ibw as _ibw
    from ..io import spm as _spm

    if _spm.is_nanoscope_file(filename):
        channel = getattr(settings, 'spm_channel', None) or _spm.DEFAULT_CHANNEL
        direction = getattr(settings, 'spm_direction', None) or _spm.DEFAULT_DIRECTION
        return _spm.scan_y_size_um(filename, channel=channel, direction=direction)
    if _ibw.is_igor_wave_file(filename):
        return _ibw.scan_y_size_um(filename)
    return _text_scan_y_size(filename)


def _text_scan_y_size(filename):
    """The ``# Height: 0.333 µm`` line of a Gwyddion text export, or None."""
    import re

    try:
        with open(filename, encoding="utf-8", errors="replace") as handle:
            header = [handle.readline() for _ in range(10)]
    except OSError:
        return None
    for line in header:
        match = re.match(r"\s*#?\s*Height:\s*(\d+\.?\d*(?:[eE][+-]?\d+)?)\s*"
                         r"(µm|um|micron|nm)", line, re.IGNORECASE)
        if match:
            value = float(match.group(1))
            return value / 1000.0 if match.group(2).lower() == 'nm' else value
    return None


def _load_text(filename, default_scan_size=2.0):
    """
    Load a text export, such as one from Gwyddion.

    The scan size is read from the header comments when they carry it. The
    behaviour is unchanged: this is the original load_afm_data body, moved so
    the dispatcher above stays readable. The stored baselines depend on it.
    """
    import re

    profile = _load_text_profile(filename)
    if profile is not None:
        return profile

    scan_x_size = default_scan_size
    
    # Try to read scan size from header
    try:
        with open(filename, 'r', encoding="utf-8", errors="replace") as f:
            # Read first 10 lines to look for scan size info
            header_lines = [f.readline() for _ in range(10)]
            
            for line in header_lines:
                # Common AFM file formats include scan size info
                # Look for patterns like "Width: 2.000 µm"
                if 'width' in line.lower():
                    match = re.search(r'(\d+\.?\d*)\s*(um|µm|micron)', line, re.IGNORECASE)
                    if match:
                        scan_x_size = float(match.group(1))
                        print(f"  Detected scan width from header: {scan_x_size} µm")
                        break
                    
                    match = re.search(r'(\d+\.?\d*)\s*nm', line, re.IGNORECASE)
                    if match:
                        scan_x_size = float(match.group(1)) / 1000  # Convert nm to µm
                        print(f"  Detected scan width from header: {scan_x_size} µm ({match.group(1)} nm)")
                        break
                        
                elif 'scan' in line.lower() and 'size' in line.lower():
                    # Extract number
                    match = re.search(r'(\d+\.?\d*)\s*(um|µm|micron)', line, re.IGNORECASE)
                    if match:
                        scan_x_size = float(match.group(1))
                        print(f"  Detected scan size from header: {scan_x_size} µm")
                        break
    except Exception:
        pass
    
    # Load the data
    try:
        # Another program's file: only the numbers are parsed, so an
        # undecodable byte in the skipped header must not stop the load.
        with open(filename, encoding="utf-8", errors="replace") as handle:
            data = np.genfromtxt(handle, skip_header=4, skip_footer=0)
    except Exception as e:
        raise Exception(f"Error loading {filename}: {e}") from e
    
    if scan_x_size == default_scan_size:
        print(f"  Using default scan size: {scan_x_size} µm")
    
    return data, scan_x_size


_PROFILE_UNITS = {'m': 1e6, 'mm': 1e3, 'um': 1.0, 'µm': 1.0, 'nm': 1e-3}


def _load_text_profile(filename):
    """
    A two-column ``x z`` line profile, as one 1-row image, or None.

    Gwyddion's profile export writes one line per sample: position along the
    profile and height, both in metres unless the header names a unit in
    brackets (``x [nm]``). Every caller downstream expects a 2-D height image
    and a scan width, so the profile is handed back in that shape: one scan
    line, and the width it actually spans, ``x[-1] - x[0]``. There is no
    header width to trust here; the x column is the measurement.

    The x column must be evenly spaced, because ``raw_data`` rebuilds the axis
    from the width and the sample count. A profile that is not would be
    silently stretched, so it is refused instead.

    Returns None for anything that is not exactly two numbers per data line,
    so the image path below sees every file it saw before.

    Ported from Aditi's version of the groove-averaging script (2026-10-27),
    which handled both formats and took the width from the x column.
    """
    import re

    try:
        with open(filename, encoding="utf-8", errors="replace") as handle:
            lines = handle.readlines()
    except OSError:
        return None     # the image path reports it, as it always has

    header, rows = [], []
    for line in lines:
        fields = line.split()
        if not fields:
            continue
        try:
            values = [float(f) for f in fields]
        except ValueError:
            if rows:
                return None     # text after the data: not this format
            header.append(line)
            continue
        if len(values) != 2:
            return None
        rows.append(values)

    if len(rows) < 3:
        return None

    xz = np.asarray(rows)
    x, z = xz[:, 0], xz[:, 1]

    # Units from a bracketed header label, x's first and z's second; metres
    # when the header says nothing, which is what Gwyddion writes.
    units = re.findall(r'\[\s*(m|mm|um|µm|nm)\s*\]', ''.join(header), re.IGNORECASE)
    x_unit = units[0].lower() if units else 'm'
    z_unit = units[1].lower() if len(units) > 1 else x_unit
    x_um = x * _PROFILE_UNITS[x_unit]
    z_m = z * _PROFILE_UNITS[z_unit] * 1e-6

    step = np.diff(x_um)
    if not (np.all(step > 0) and np.allclose(step, step.mean(), rtol=1e-3)):
        raise ValueError(
            f"{filename}: the x column of this profile is not evenly spaced "
            f"(steps {step.min():.4g} to {step.max():.4g} µm), and the groove "
            "averaging assumes it is. Re-export the profile from the image.")

    scan_x_size = float(x_um[-1] - x_um[0])
    print(f"  Line profile: {len(z)} points over {scan_x_size:g} µm "
          f"(x in {x_unit}, height in {z_unit})")
    return z_m[np.newaxis, :], scan_x_size


def raw_data(data, scan_x_size):
    """
    Extract displacement and height profile from 2D AFM data (averages all rows)

    NaN pixels (masked, or on a dropped scar line) are left out of each
    column's mean; a column NaN on every line stays NaN. See core/region_mask.py.
    """
    from .region_mask import nan_mean_rows
    disp_um = scan_x_size * np.arange(0, len(data[0])) / (len(data[0]) - 1)
    mean = nan_mean_rows(data)
    groov_nm = (10**9) * (mean - np.nanmin(mean))
    return disp_um, groov_nm


def raw_data_multi_group(data, scan_x_size, n_groups=10):
    """
    Extract multiple height profiles from different row groups in 2D AFM data
    
    This function divides the AFM image into n_groups horizontal bands and extracts
    an averaged profile from each band. This allows for many more measurements per
    image and enables assessment of within-image variation.
    
    Parameters:
        data: 2D array of AFM height data (rows × columns)
        scan_x_size: scan width in microns
        n_groups: number of row groups to create (default: 10)
        
    Returns:
        disp_um: 1D array of x-positions in microns (same for all groups)
        profiles_nm: list of n_groups 1D arrays, each containing a height profile in nm
        group_info: dict with metadata about the grouping
    """
    from .region_mask import nan_mean_rows

    # Lines that are NaN from end to end (masked rows, dropped scar lines) are
    # left out before the bands are cut, so they contribute to no group and the
    # groups stay equal in size over what remains. Their indices still refer to
    # the image as loaded.
    all_rows = data.shape[0]
    kept_rows = np.flatnonzero(~np.all(np.isnan(data), axis=1))
    rows_excluded = all_rows - len(kept_rows)
    if rows_excluded:
        data = data[kept_rows]

    n_rows = data.shape[0]
    n_cols = data.shape[1]
    if n_rows == 0:
        raise ValueError("every scan line is masked or dropped; nothing left to analyse")
    
    # Calculate x-displacement (same for all profiles)
    disp_um = scan_x_size * np.arange(0, n_cols) / (n_cols - 1)
    
    # Determine row ranges for each group
    rows_per_group = n_rows // n_groups
    
    if rows_per_group < 3:
        print(f"Warning: Only {rows_per_group} rows per group with n_groups={n_groups}")
        print(f"  Reducing to n_groups={n_rows // 3} to ensure at least 3 rows per group")
        n_groups = max(1, n_rows // 3)
        rows_per_group = n_rows // n_groups
    
    profiles_nm = []
    group_ranges = []
    
    for i in range(n_groups):
        # Calculate row range for this group
        start_row = i * rows_per_group
        
        # For the last group, include any remaining rows
        if i == n_groups - 1:
            end_row = n_rows
        else:
            end_row = (i + 1) * rows_per_group
        
        # Extract and average rows in this group
        group_data = data[start_row:end_row, :]
        profile_avg = nan_mean_rows(group_data)
        
        # Convert to nm and subtract minimum
        profile_nm = (10**9) * (profile_avg - np.nanmin(profile_avg))
        
        profiles_nm.append(profile_nm)
        # In the image's own line numbers, so a masked gap shows as a wider band
        group_ranges.append((int(kept_rows[start_row]), int(kept_rows[end_row - 1]) + 1))
    
    # Package metadata
    group_info = {
        'n_groups': n_groups,
        'n_rows': n_rows,
        'rows_per_group': rows_per_group,
        'group_ranges': group_ranges,
        'rows_excluded': rows_excluded,
    }
    
    print(f"  Extracted {n_groups} row-group profiles")
    print(f"    Rows per group: {rows_per_group}")
    print(f"    Total rows used: {n_groups * rows_per_group} of {n_rows}")
    if rows_excluded:
        print(f"    Lines excluded (masked or dropped): {rows_excluded} of {all_rows}")
    
    return disp_um, profiles_nm, group_info


@dataclass(frozen=True)
class FlatteningReport:
    """
    What `flatten_profile` did, as against what it was asked to do.

    Flattening sets the plane every blaze angle on a scan is quoted against, so
    a row group whose background was fitted by a fallback is not comparable with
    one that got the method asked for. The warnings say so at the time; this is
    the same information in a form the caller can store beside the result.

    ``residual_rms_nm`` is the scatter of the fitted features about the
    background, in nanometres: how far from a smooth low-order surface the lands
    actually lie, and so how much the reference can be trusted. It is ``None``
    for `linear` and `polynomial`, which fit the whole trace rather than
    features -- their residual would be dominated by the groove depth and would
    say nothing about flatness.

    Read it with ``n_features`` beside it. A fit through exactly
    ``poly_order + 1`` points passes through them all and reports a residual of
    zero, which is arithmetic rather than a flat surface; ``dof`` is how many
    points beyond that the fit had, and a residual with ``dof == 0`` means
    nothing.
    """

    method_requested: str
    method_used: str
    poly_order_requested: int
    poly_order_used: int
    n_features: int
    edge_exclusion_skipped: bool
    residual_rms_nm: Optional[float]

    @property
    def dof(self) -> int:
        """Features beyond the minimum the fitted order needed."""
        if not self.n_features:
            return 0
        return max(0, self.n_features - (self.poly_order_used + 1))

    @property
    def fell_back(self) -> bool:
        """True when the background is not the one that was asked for."""
        return (self.method_used != self.method_requested
                or self.poly_order_used != self.poly_order_requested
                or self.edge_exclusion_skipped)

    def summary(self) -> str:
        """One line for a log or a summary file."""
        text = self.method_used
        if self.method_used != self.method_requested:
            text += f" (asked for {self.method_requested})"
        if self.poly_order_used != self.poly_order_requested:
            text += f", degree {self.poly_order_used} not {self.poly_order_requested}"
        if self.edge_exclusion_skipped:
            text += ", edge exclusion skipped"
        if self.residual_rms_nm is not None:
            text += (f", {self.n_features} features, residual "
                     f"{self.residual_rms_nm:.3f} nm")
            if not self.dof:
                text += " (exact fit, meaningless)"
        return text


def _feature_residual_rms(feature_x, feature_y, coeffs):
    """RMS of the fitted features about the background they produced, in nm."""
    residual = np.asarray(feature_y, dtype=float) - np.polyval(coeffs, feature_x)
    return float(np.sqrt(np.mean(residual ** 2)))


def flatten_profile(x, y, method='linear', poly_order=2, exclude_edges=0.0,
                   feature='both', period_nm=None, return_report=False,
                   exclude=None):
    """
    Remove background tilt/curvature from AFM profile
    
    Parameters:
        method: 'linear', 'polynomial', 'groove_peaks', or 'level_grooves'
        poly_order: polynomial order (for 'polynomial' and 'level_grooves' methods)
        exclude_edges: fraction of data to exclude from edges (0.0 to 0.5).
            For 'level_grooves' it is skipped, with a warning, when excluding
            the edges would leave too few features to fit the requested order.
        feature: for 'level_grooves': 'peaks', 'troughs', or 'both'
        period_nm: estimated period for finding grooves (for 'level_grooves')
        return_report: also return a :class:`FlatteningReport` describing what
            was actually done, following the `return_n_edge_rejected` keyword
            on `find_groove_positions`. Off by default so existing two-value
            unpacking keeps working.
        exclude: boolean array, True where ``y`` is not a measurement: the
            straight bridge across a masked gap (core/region_mask.py). Those
            points, and any peak or trough within two pixels of them, take no
            part in any fit; the background is still evaluated everywhere. A
            bridge left in pulled every method's background, by 1.3 degrees of
            blaze angle for `linear` on a scan with 2.5 periods masked.
        
    Returns:
        y_flat: flattened profile
        background: the removed background
        report: a FlatteningReport, only when `return_report` is true
    """
    # What was actually done, for the report. Each branch narrows these.
    method_used = method
    order_used = poly_order
    n_features = 0
    edge_exclusion_skipped = False
    residual_rms = None

    # Points the fits may use. Everything, unless part of the trace is a bridge.
    usable = np.ones(len(y), dtype=bool)
    if exclude is not None and np.any(exclude):
        gap = np.asarray(exclude, dtype=bool)
        near_gap = np.convolve(gap.astype(float), np.ones(5), mode='same') > 0
        usable = ~near_gap

    def _features(indices):
        """Drop features on or beside a bridge."""
        indices = np.asarray(indices, dtype=int)
        return indices[usable[indices]] if len(indices) else indices

    if method == 'linear':
        poly_coeffs = np.polyfit(x[usable], y[usable], 1)
        background = np.polyval(poly_coeffs, x)
        # `linear` is degree 1 by definition, so poly_order was never a request
        # here and reporting a reduction would make every call look degraded.
        order_used = poly_order
        
    elif method == 'polynomial':
        if exclude_edges > 0:
            n_exclude = int(len(x) * exclude_edges)
            keep = usable[n_exclude:-n_exclude]
            x_fit = x[n_exclude:-n_exclude][keep]
            y_fit = y[n_exclude:-n_exclude][keep]
            poly_coeffs = np.polyfit(x_fit, y_fit, poly_order)
        else:
            poly_coeffs = np.polyfit(x[usable], y[usable], poly_order)
        
        background = np.polyval(poly_coeffs, x)
        
    elif method == 'groove_peaks':
        peaks, _ = find_peaks(y, distance=len(y)//20)
        peaks = _features(peaks)
                
        if len(peaks) < 3:
            warnings.warn(
                f"groove_peaks found {len(peaks)} peak(s), too few to fit; "
                f"using a degree-{poly_order} polynomial through the whole profile",
                UserWarning, stacklevel=2)
            poly_coeffs = np.polyfit(x[usable], y[usable], poly_order)
            background = np.polyval(poly_coeffs, x)
            method_used = 'polynomial'
        else:
            poly_coeffs = np.polyfit(x[peaks], y[peaks], poly_order)
            background = np.polyval(poly_coeffs, x)
            n_features = len(peaks)
            residual_rms = _feature_residual_rms(x[peaks], y[peaks], poly_coeffs)
            
    elif method == 'level_grooves':
        if period_nm is None:
            raise ValueError("period_nm must be provided for level_grooves method")
        
        feature_x = []
        feature_y = []
        
        if feature in ['peaks', 'both']:
            dx_nm = (x[1] - x[0]) * 1000
            min_distance = int(0.5 * period_nm / dx_nm)  # Changed from 0.3 to 0.5 for more reliable detection
            # Use a lower prominence to catch more peaks
            height_range = np.max(y) - np.min(y)
            min_prominence = 0.2 * height_range  # Require peaks to be at least 20% of height range
            peaks, _ = find_peaks(y, distance=min_distance, prominence=min_prominence)
            peaks = _features(peaks)
            if len(peaks) > 0:
                feature_x.extend(x[peaks])
                feature_y.extend(y[peaks])
        
        if feature in ['troughs', 'both']:
            y_inv = -y
            dx_nm = (x[1] - x[0]) * 1000
            min_distance = int(0.5 * period_nm / dx_nm)  # Changed from 0.3 to 0.5
            height_range = np.max(y) - np.min(y)
            min_prominence = 0.2 * height_range
            troughs, _ = find_peaks(y_inv, distance=min_distance, prominence=min_prominence)
            troughs = _features(troughs)
            if len(troughs) > 0:
                feature_x.extend(x[troughs])
                feature_y.extend(y[troughs])
        
        if len(feature_x) < 2:
            warnings.warn(
                f"level_grooves found {len(feature_x)} feature(s), too few to fit; "
                "using linear",
                UserWarning, stacklevel=2)
            poly_coeffs = np.polyfit(x[usable], y[usable], 1)
            background = np.polyval(poly_coeffs, x)
            method_used = 'linear'
            order_used = 1
        else:
            feature_x = np.array(feature_x)
            feature_y = np.array(feature_y)
            
            if exclude_edges > 0:
                x_min = x[0] + exclude_edges * (x[-1] - x[0])
                x_max = x[-1] - exclude_edges * (x[-1] - x[0])
                mask = (feature_x >= x_min) & (feature_x <= x_max)
                kept = int(np.count_nonzero(mask))
                # Edge exclusion is meant to improve the fit, so it must not be
                # allowed to starve it. This used to keep the remnant and fit
                # that while printing "using all": with two features left and a
                # quadratic asked for, the background was rank-deficient and
                # still looked like a fit.
                if kept >= poly_order + 1:
                    feature_x = feature_x[mask]
                    feature_y = feature_y[mask]
                else:
                    warnings.warn(
                        f"edge exclusion leaves {kept} feature(s), too few for a "
                        f"degree-{poly_order} background; using all of them",
                        UserWarning, stacklevel=2)
                    edge_exclusion_skipped = True

            # A degree-k polynomial needs k+1 points to be determined. Reaching
            # here guarantees two, so the order is reduced rather than letting
            # polyfit return a rank-deficient answer behind a RankWarning.
            order = min(poly_order, len(feature_x) - 1)
            if order < poly_order:
                warnings.warn(
                    f"{len(feature_x)} feature(s) support only a degree-{order} "
                    f"background, not degree-{poly_order}",
                    UserWarning, stacklevel=2)

            poly_coeffs = np.polyfit(feature_x, feature_y, order)
            background = np.polyval(poly_coeffs, x)
            order_used = order
            n_features = len(feature_x)
            residual_rms = _feature_residual_rms(feature_x, feature_y, poly_coeffs)
    else:
        raise ValueError(f"Unknown flattening method: {method}")
    
    y_flat = y - background

    if return_report:
        report = FlatteningReport(
            method_requested=method,
            method_used=method_used,
            poly_order_requested=poly_order,
            poly_order_used=order_used,
            n_features=n_features,
            edge_exclusion_skipped=edge_exclusion_skipped,
            residual_rms_nm=residual_rms,
        )
        return y_flat, background, report

    return y_flat, background


def find_groove_positions(x, y, period_nm, prominence_factor=0.3, distance_factor=0.7,
                          edge_exclusion=0.0, return_n_edge_rejected=False):
    """
    Find positions of groove centers (minima)

    Parameters:
        edge_exclusion: reject centers lying within this many periods of either
            end of the profile. A groove that close to the edge has no complete
            facet in the data (the scan starts or stops part-way through it),
            so any angle fitted to it is meaningless. 0.0 disables the check.
        return_n_edge_rejected: also return how many centers were rejected

    Returns:
        peaks, or (peaks, n_rejected) if return_n_edge_rejected
    """
    y_inv = -y  # Invert so grooves become peaks
    dx_nm = (x[1] - x[0]) * 1000
    period_idx = int(period_nm / dx_nm)
    height_range = np.max(y) - np.min(y)
    min_prominence = prominence_factor * height_range

    peaks, properties = find_peaks(y_inv,
                                   distance=int(distance_factor * period_idx),
                                   prominence=min_prominence)

    n_rejected = 0
    if edge_exclusion > 0 and len(peaks) > 0:
        margin_idx = edge_exclusion * period_idx
        last = len(x) - 1
        keep = (peaks >= margin_idx) & (peaks <= last - margin_idx)
        n_rejected = int((~keep).sum())
        peaks = peaks[keep]

    if return_n_edge_rejected:
        return peaks, n_rejected
    return peaks


def extract_single_groove(x, y, center_idx, period_nm, margin=0.1, allow_asymmetric=True):
    """
    Extract a single groove around center_idx
    
    Parameters:
        allow_asymmetric: If True, allows asymmetric extraction near edges
    """
    dx_nm = (x[1] - x[0]) * 1000
    half_width = int(round((period_nm / 2) * (1 + margin) / dx_nm))
    
    if allow_asymmetric:
        # Allow asymmetric extraction for edge grooves
        left_extent = min(half_width, center_idx)
        right_extent = min(half_width, len(x) - 1 - center_idx)
        
        start = center_idx - left_extent
        end = center_idx + right_extent + 1
    else:
        # Original symmetric extraction
        L = min(half_width, center_idx, len(x) - 1 - center_idx)
        start = center_idx - L
        end = center_idx + L + 1
    
    x_groove = x[start:end]
    y_groove = y[start:end]
    
    return x_groove, y_groove
