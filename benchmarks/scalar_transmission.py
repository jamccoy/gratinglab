"""Two scalar transmission phases against the RCWA, for the #83 design.

Produces the tables in ``docs/theory/scalar-transmission.md`` §5. Run from the
repository root::

    python benchmarks/scalar_transmission.py            # both tables, ~1 min
    python benchmarks/scalar_transmission.py --normal   # normal incidence only

The grating is a relief cut into a lossless n = 1.5 table, vacuum above, and
the light is transmitted into the glass at 500 nm, TE, in plane. The table is
written to a temporary folder named by ``GRATINGLAB_MATERIALS``, so your own
library is not read or changed.

Neither phase is in the package yet. Both are computed here:

``thin``
    order-free, k g(t) (n cos b0 - cos a) with b0 the zero order's Snell angle.
``refr``
    order-dependent, k g(t) (n cos b_m - cos a), the transmissive twin of
    ScalarSolver's reflective k g (cos a + cos b_m).

Both carry the flat-interface TE Fresnel flux factor
4 cos a n cos b_m / (cos a + n cos b_m)^2. The phase enters with a negative
sign, which puts a blazed grating's peak on the same order as the RCWA's; that
was matched by which order peaks, not derived from ``docs/conventions.md``.
"""

from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

import numpy as np

N = 1.5
WL = 500.0
ORDERS = (-2, -1, 0, 1, 2)
SAMPLES = 4096


def scalar_efficiency(problem, alpha_deg: float, m: int, model: str) -> float:
    """Transmitted efficiency of order ``m`` under one of the two phases."""
    t = (np.arange(SAMPLES) + 0.5) / SAMPLES
    g = np.asarray(problem.height_nm(t), float)
    k = 2 * np.pi / WL
    ca, sa = np.cos(np.radians(alpha_deg)), np.sin(np.radians(alpha_deg))
    s = (sa + m * WL / problem.period) / N
    if abs(s) >= 1:
        return 0.0  # evanescent in the glass
    cb = np.sqrt(1 - s * s)
    cb0 = np.sqrt(1 - (sa / N) ** 2)
    phase = k * g * ((N * cb0 - ca) if model == "thin" else (N * cb - ca))
    amplitude = np.mean(np.exp(-1j * phase) * np.exp(-2j * np.pi * m * t))
    return float(abs(amplitude) ** 2 * 4 * ca * N * cb / (ca + N * cb) ** 2)


def klein_cook(problem) -> float:
    """Q = 2 pi lambda h / (n period^2): thin well below 1, thick above."""
    return 2 * np.pi * WL * problem.depth / (N * problem.period**2)


def truncation(period: float) -> int:
    return 40 if period < 5000 else 90


def blazed(period: float):
    """A blaze deep enough that a thin prism sends normal light into one order."""
    from gratinglab import Blazed, Problem

    theta = float(np.degrees(np.arctan(WL / (N - 1) / period)))
    return Problem(period=period, profile=Blazed(blaze_angle=theta), coating="Glass15")


def normal_incidence(rcwa) -> None:
    from gratinglab import Illumination, Lamellar, Problem

    light = Illumination(alpha_deg=0.0, gamma_deg=90.0, polarization="TE")
    print("Normal incidence: transmitted efficiency, rcwa / thin / refr")
    for period in (2500.0, 10000.0):
        cases = [
            (
                f"lamellar h={h:g}",
                Problem(
                    period=period,
                    coating="Glass15",
                    profile=Lamellar(depth_fraction=h / period, duty_cycle=0.5),
                ),
            )
            for h in (100.0, 500.0, 1000.0)
        ]
        cases.append((f"blazed h={WL / (N - 1):g}", blazed(period)))
        for label, problem in cases:
            scan = rcwa.solve(problem, light, [WL], truncation=truncation(period))
            tr = scan.transmission
            print(
                f"\n  {label}, period {period:g} nm, Q = {klein_cook(problem):.3f}"
                f"  (rcwa R {scan.total[0]:.4f}, T {tr.total[0]:.4f})"
            )
            for m in ORDERS:
                print(
                    f"    {m:+d}  {tr.order(m)[0]:.4f}"
                    f"  {scalar_efficiency(problem, 0.0, m, 'thin'):.4f}"
                    f"  {scalar_efficiency(problem, 0.0, m, 'refr'):.4f}"
                )


def oblique(rcwa) -> None:
    from gratinglab import Illumination

    print("\nOblique incidence: the blaze peak, rcwa / thin / refr")
    for period in (2500.0, 10000.0):
        problem = blazed(period)
        for alpha in (-30.0, 30.0):
            light = Illumination(alpha_deg=alpha, gamma_deg=90.0, polarization="TE")
            tr = rcwa.solve(
                problem, light, [WL], truncation=truncation(period)
            ).transmission
            peak = max(range(-3, 4), key=lambda m: tr.order(m)[0])
            print(
                f"  period {period:g} nm, alpha {alpha:+g}, order {peak:+d}:"
                f"  {tr.order(peak)[0]:.4f}"
                f"  {scalar_efficiency(problem, alpha, peak, 'thin'):.4f}"
                f"  {scalar_efficiency(problem, alpha, peak, 'refr'):.4f}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--normal", action="store_true", help="normal incidence only")
    args = parser.parse_args()

    folder = Path(tempfile.mkdtemp(prefix="gratinglab-bench-"))
    (folder / "Glass15.nk").write_text(
        "# constant index, for benchmarks/scalar_transmission.py\n"
        "wavelength (nm), n, k\n300, 1.5, 0.0\n800, 1.5, 0.0\n",
        encoding="utf-8",
    )
    os.environ["GRATINGLAB_MATERIALS"] = str(folder)

    from gratinglab import get_solver

    rcwa = get_solver("rcwa")
    normal_incidence(rcwa)
    if not args.normal:
        oblique(rcwa)


if __name__ == "__main__":
    main()
