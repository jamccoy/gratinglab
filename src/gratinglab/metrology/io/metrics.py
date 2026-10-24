"""
The slope/curvature summary that accompanies a boundary profile.

This used to sit beside a second copy of the ``.ggp`` writer. That copy has
been removed, and :mod:`gratinglab.io.ggp` is now the only implementation, but
the metrics file has no counterpart there and is metrology's own output, so it
stays here.

The file is a human-readable sidecar and is not an interchange format; nothing
parses it back. The numbers that a solver needs travel through
:meth:`~gratinglab.metrology.boundary.BoundaryProfile.to_problem`.
"""

from ..core.tip import format_angles


def write_profile_metrics(path, metrics):
    """Write the slope/curvature summary that accompanies a boundary profile"""
    with open(path, 'w', encoding="utf-8") as f:
        f.write("=== Groove Profile Analysis (Final Phase-Shifted Data) ===\n")
        f.write(f"Period: {metrics['period_nm']:.2f} nm\n")
        f.write(f"Grooves averaged: {metrics['n_grooves']}\n")
        f.write(f"Groove depth: {metrics['groove_depth']:.4f} (fraction of period)\n")
        f.write(f"Peak-to-valley: {metrics['peak_to_valley']:.4f} (fraction of period)\n")
        f.write(f"RMS slope: {metrics['rms_slope']:.4f} (normalized units)\n")
        f.write(f"Max slope magnitude: {metrics['max_slope']:.4f}\n")
        f.write(f"Max curvature: {metrics['max_curvature']:.6f} (normalized units)\n")
        # Present only when the pipeline ran a tip correction. Stated even so:
        # a corrected depth and an uncorrected one are different measurements,
        # and this sidecar is where a reader learns which this file holds.
        # Present when the caller measured the groove orientation. A period
        # along the scan line and one along the groove normal differ by
        # 1/cos(theta), and this line says which the file holds.
        if 'groove_orientation' in metrics:
            angle = metrics.get('groove_rotation_deg')
            if angle is None:
                f.write(f"Groove rotation: not measured "
                        f"({metrics.get('groove_rotation_note', '')})\n")
            else:
                done = ("corrected" if metrics.get('groove_rotation_corrected')
                        else f"not corrected, {metrics.get('groove_rotation_note', '')}")
                f.write(f"Groove rotation: {angle:+.3f} +/- "
                        f"{metrics['groove_rotation_stderr_deg']:.3f} deg "
                        f"({done})\n")
        if metrics.get('tip_correction', 'none') != 'none':
            f.write(f"Tip correction: {metrics['tip_correction']} "
                    f"(R = {metrics['tip_radius_nm']:g} nm, "
                    f"{format_angles(metrics['tip_half_angle_deg'], metrics.get('tip_half_angle_right_deg'))})\n")
            f.write(f"Tip-certain pixels: "
                    f"{100.0 * metrics['tip_certain_fraction']:.1f}% "
                    f"(the rest are upper bounds on the surface)\n")
        else:
            f.write("Tip correction: none\n")
        # Only when used, so an unrepaired, unmasked export is byte-identical
        # to what this file always held.
        if metrics.get('scar_repair', 'none') != 'none':
            f.write(f"Scar repair: {metrics['scar_repair']} "
                    f"({metrics.get('scar_lines', 0)} line(s))\n")
        if metrics.get('region_mask'):
            f.write(f"Region mask: {metrics['region_mask']} "
                    f"({metrics.get('masked_pixels', 0)} pixels, "
                    f"{metrics.get('mask_rejected_grooves', 0)} groove(s) "
                    f"rejected)\n")
    return path
