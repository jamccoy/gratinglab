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
    return path
