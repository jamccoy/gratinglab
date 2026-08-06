"""
PCGrate boundary profile (.ggp) output

A .ggp file is two header lines followed by one "x y" pair per line, both
normalised: x runs 0 to 1 across exactly one period, y is height as a fraction
of that period.

    3 0 - Polygonal type
    Period: 1 PSC: 1
    0.000000 0.000000
    0.000500 0.000001
    ...

The header must NOT be commented. Writing these files with np.savetxt(header=...)
prepends "# ", which PCGrate does not accept - existing .ggp files in this project
carry three different header variants, and only the uncommented two-line form is
correct. write_ggp always emits that form.
"""
import os

GGP_HEADER = "3 0 - Polygonal type\nPeriod: 1 PSC: 1"


def write_ggp(path, x_normalized, y_normalized, fmt='%f'):
    """
    Write a PCGrate .ggp boundary profile.

    Parameters:
        x_normalized: positions in [0, 1] spanning one period
        y_normalized: heights as a fraction of the period
    """
    if len(x_normalized) != len(y_normalized):
        raise ValueError(f"x and y differ in length: "
                         f"{len(x_normalized)} vs {len(y_normalized)}")

    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)

    with open(path, 'w') as f:
        f.write(GGP_HEADER + "\n")
        for xi, yi in zip(x_normalized, y_normalized):
            f.write(f"{fmt % xi} {fmt % yi}\n")

    return path


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
