"""The scalar solver's reflectivity models against the rigorous solvers.

Produces the tables in ``docs/findings.md``, "The groove-resolved reflectivity
model against the rigorous solvers". Run from the repository root::

    python benchmarks/coated_scalar_reflectivity.py              # all three, ~3 min
    python benchmarks/coated_scalar_reflectivity.py --case critical

Each model is isolated by dividing out the groove's diffraction:

    scalar    eta(coated, model) / eta(uncoated)
    rigorous  eta(coated) / eta(perfect conductor), averaged over TE and TM

so most of Kirchhoff's shape error cancels and what is left is how each
handles the coating across the groove. The coating is ``typical:Au``, the
built-in model, so no table is needed; both sides see the same index.

Cases:

``grazing``
    In plane, alpha = 80 deg, shallow sinusoid, lambda/p <= 0.04. The local
    graze runs from about 1 to 19 deg, well past gold's critical angle. No
    point is shadowed. Rigorous side: the C-method.
``critical``
    In plane, alpha = 87 deg, lambda/p = 0.016 to 0.032. The local graze, 1.2
    to 4.8 deg, straddles gold's critical angle at 1.2 nm (about 3.6 deg).
    Rigorous side: the C-method; truncation 80, 120 and 180 agree to three
    decimals.
``reference``
    ``docs/theory/scalar.md`` section 9's geometry, Blazed 29.5/70.5, gamma =
    1.25 deg, alpha = 19.99 deg. Outside the scalar regime: the reduced ratio
    lambda/(p sin gamma) is 0.22 to 0.65. Rigorous side: the integral solver;
    512 and 1536 boundary points agree to about 1e-4. Both scalar runs use
    visibility="horizon", so shadowing cancels; ``facet`` has no masks and is
    divided by the unmasked bare run. ``average`` is omitted here, because its
    exit side sees no horizon and the ratio would not isolate the reflectivity.
"""

from __future__ import annotations

import argparse

import numpy as np

from gratinglab import Blazed, Illumination, Problem
from gratinglab.compare import sweep
from gratinglab.profiles import Sinusoidal
from gratinglab.solvers import scalar

COATING = "typical:Au"
STRONGEST = 4  # orders reported per wavelength, by perfect-conductor efficiency


def rigorous_ratio(coated, bare, illumination_for, wavelengths, method, coated_opts, bare_opts):
    """Coated / perfect-conductor efficiency, mean over TE and TM."""
    eff = {}
    for pol in ("TE", "TM"):
        (c,) = sweep(coated, illumination_for(pol), wavelengths, [method], {method: coated_opts})
        (b,) = sweep(bare, illumination_for(pol), wavelengths, [method], {method: bare_opts})
        eff[pol] = (c, b)
    orders = eff["TE"][1].orders
    coated_mean = 0.5 * (eff["TE"][0].efficiency + eff["TM"][0].efficiency)
    bare_mean = 0.5 * (eff["TE"][1].efficiency + eff["TM"][1].efficiency)
    return orders, coated_mean, bare_mean


def report(label, wavelengths, orders, coated, bare, columns):
    """Print, per wavelength, the strongest orders and each model's ratio."""
    print(f"\n== {label}")
    names = "  ".join(f"{name:>8s}" for name in columns)
    for i, wl in enumerate(wavelengths):
        total = coated[i].sum() / bare[i].sum()
        print(f"\n  lambda {wl:g} nm, rigorous total coated/PEC {total:.4f}")
        print(f"     m   eta_PEC  rigorous  {names}")
        for j in sorted(np.argsort(bare[i])[::-1][:STRONGEST], key=lambda j: orders[j]):
            m = int(orders[j])
            values = "  ".join(f"{column(i, m):8.4f}" for column in columns.values())
            print(f"   {m:+3d}  {bare[i, j]:7.4f}  {coated[i, j] / bare[i, j]:8.4f}  {values}")


def scalar_columns(coated, bare, wavelengths, illumination, models):
    """{name: f(i, m)} giving each scalar model's coated/uncoated ratio."""
    columns = {}
    for name, (model, visibility, bare_visibility) in models.items():
        kwargs = {"reflectivity_model": model}
        if visibility:
            kwargs["visibility"] = visibility
        s = scalar.solve(coated, illumination, wavelengths, **kwargs)
        b = scalar.solve(bare, illumination, wavelengths,
                         **({"visibility": bare_visibility} if bare_visibility else {}))

        def ratio(i, m, s=s, b=b):
            k = int(np.flatnonzero(s.orders == m)[0])
            return s.efficiency[i, k] / b.efficiency[i, k]

        columns[name] = ratio
    return columns


def in_plane(label, period, depth_fraction, alpha, wavelengths, truncation):
    profile = Sinusoidal(depth_fraction=depth_fraction)
    coated = Problem(period=period, profile=profile, coating=COATING)
    bare = Problem(period=period, profile=profile)
    orders, c, b = rigorous_ratio(
        coated, bare, lambda pol: Illumination.classical(alpha=alpha, polarization=pol),
        wavelengths, "cmethod",
        {"conductivity": "tabulated", "truncation": truncation},
        {"conductivity": "perfect", "truncation": truncation},
    )
    models = {name: (name, None, None) for name in ("local", "average", "facet")}
    columns = scalar_columns(
        coated, bare, wavelengths,
        Illumination.classical(alpha=alpha, polarization="unpolarized"), models,
    )
    report(label, wavelengths, orders, c, b, columns)


def reference(wavelengths, points):
    profile = Blazed(blaze_angle=29.5, antiblaze_angle=70.5)
    coated = Problem(period=315.15, profile=profile, coating=COATING)
    bare = Problem(period=315.15, profile=profile)
    orders, c, b = rigorous_ratio(
        coated, bare,
        lambda pol: Illumination(alpha_deg=19.99, gamma_deg=1.25, polarization=pol),
        wavelengths, "integral",
        {"conductivity": "tabulated", "boundary_points": points},
        {"conductivity": "perfect", "boundary_points": points},
    )
    models = {"local": ("local", "horizon", "horizon"), "facet": ("facet", None, None)}
    columns = scalar_columns(
        coated, bare, wavelengths,
        Illumination(alpha_deg=19.99, gamma_deg=1.25, polarization="unpolarized"), models,
    )
    report("reference: Blazed 29.5/70.5, gamma 1.25 deg (outside the scalar regime)",
           wavelengths, orders, c, b, columns)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--case", choices=["grazing", "critical", "reference", "all"],
                        default="all")
    args = parser.parse_args()
    if args.case in ("grazing", "all"):
        in_plane("grazing: alpha 80 deg, past the critical angle",
                 100.0, 0.02, 80.0, np.array([2.0, 4.0]), 60)
    if args.case in ("critical", "all"):
        in_plane("critical: alpha 87 deg, straddling the critical angle",
                 50.0, 0.01, 87.0, np.array([0.8, 1.2, 1.6]), 120)
    if args.case in ("reference", "all"):
        reference(np.array([1.5, 3.0, 4.5]), 512)


if __name__ == "__main__":
    main()
