"""
The slope/curvature summary that accompanies a boundary profile.

This used to sit beside a second copy of the ``.ggp`` writer. That copy is
gone -- :mod:`gratinglab.io.ggp` is now the only implementation -- but the
metrics file has no counterpart there and is metrology's own output, so it
stays here.

It is a human-readable sidecar, not an interchange format: nothing parses it
back. The numbers that a solver actually needs travel through
:meth:`~gratinglab.metrology.boundary.BoundaryProfile.to_problem`.
"""


def write_profile_metrics(path, metrics):
    """Write the slope/curvature summary that accompanies a boundary profile"""
    with open(path, 'w') as f:
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
        if metrics.get('tip_correction', 'none') != 'none':
            f.write(f"Tip correction: {metrics['tip_correction']} "
                    f"(R = {metrics['tip_radius_nm']:g} nm, "
                    f"half angle = {metrics['tip_half_angle_deg']:g} deg)\n")
            f.write(f"Tip-certain pixels: "
                    f"{100.0 * metrics['tip_certain_fraction']:.1f}% "
                    f"(the rest are upper bounds on the surface)\n")
        else:
            f.write("Tip correction: none\n")
    return path
