"""Wall time of the rigorous solvers, under the BLAS threading they inherit.

Run from the repository root::

    python benchmarks/solver_speed.py              # every case, both envs
    python benchmarks/solver_speed.py --cases rcwa-lamellar-80
    python benchmarks/solver_speed.py --env default --repeat 5

Each case runs in a fresh interpreter, because OpenBLAS reads
``OPENBLAS_NUM_THREADS`` once, when it is loaded. ``default`` inherits the
caller's environment with every BLAS thread variable removed; ``single`` sets
them to 1. Before the solvers capped their own BLAS threads the two columns
differed by up to two orders of magnitude on a loaded machine; afterwards they
should agree, since the solve no longer depends on the environment.

The times are the best of ``--repeat`` solves after one warm-up, so the
figure is the solve's cost and not the first call's imports and caches.
Machine load matters a great deal to the ``default`` column: record it
alongside any number taken from here.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

_THREAD_VARS = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)

#: ``name -> (solver, problem kind, illumination kind, wavelengths, options)``.
CASES: dict[str, tuple] = {
    "integral-512": ("integral", "blazed", "inplane", [500.0], {"boundary_points": 512}),
    "integral-1024": ("integral", "blazed", "inplane", [500.0], {"boundary_points": 1024}),
    "integral-tab-512": (
        "integral", "blazed-metal", "conical", [500.0],
        {"boundary_points": 512, "conductivity": "tabulated"},
    ),
    "rcwa-lamellar-30": ("rcwa", "lamellar", "inplane", [500.0], {"truncation": 30}),
    "rcwa-lamellar-80": ("rcwa", "lamellar", "inplane", [500.0], {"truncation": 80}),
    "rcwa-conical-30": ("rcwa", "lamellar", "conical", [500.0], {"truncation": 30}),
    "rcwa-conical-80": ("rcwa", "lamellar", "conical", [500.0], {"truncation": 80}),
    "rcwa-sine-30": (
        "rcwa", "sine", "inplane", [500.0], {"truncation": 30, "slices": 16},
    ),
    # Scans, for the workers option.
    "integral-scan-8": (
        "integral", "blazed", "inplane", list(range(400, 720, 40)),
        {"boundary_points": 512},
    ),
    "rcwa-scan-8": (
        "rcwa", "sine", "conical", list(range(400, 720, 40)),
        {"truncation": 30, "slices": 16},
    ),
}


def _child(spec: dict) -> dict:
    """Run one case in this interpreter and return its timings."""
    import time

    import numpy as np

    from gratinglab.illumination import Illumination
    from gratinglab.materials.optical import OpticalConstants, write_ari
    from gratinglab.problem import Problem
    from gratinglab.profiles import Blazed, Lamellar, Sinusoidal
    from gratinglab.solvers import get_solver

    library = Path(spec["library"])
    for name, index in (("glass", 1.5 + 0j), ("metal", 0.2 + 3.0j)):
        write_ari(
            OpticalConstants(
                name=name,
                wavelength_nm=(200.0, 1000.0),
                decrement=(1.0 - index.real,) * 2,
                absorption=(index.imag,) * 2,
                source="solver_speed",
            ),
            library / f"{name}.ari",
        )
    os.environ["GRATINGLAB_MATERIALS"] = str(library)

    solver_name, kind, mount, wavelengths, options = CASES[spec["case"]]
    options = {**options, **spec.get("extra", {})}
    if kind == "lamellar":
        problem = Problem(period=600.0, profile=Lamellar(depth_fraction=0.3), coating="glass")
    elif kind == "sine":
        problem = Problem(period=600.0, profile=Sinusoidal(depth_fraction=0.3), coating="glass")
    elif kind == "blazed-metal":
        problem = Problem(
            period=600.0, profile=Blazed(blaze_angle=20.0, antiblaze_angle=70.0),
            coating="metal",
        )
    else:
        problem = Problem(period=600.0, profile=Blazed(blaze_angle=20.0, antiblaze_angle=70.0))
    if mount == "conical":
        illumination = Illumination.offplane(graze=40.0, azimuth=25.0, polarization="unpolarized")
    else:
        illumination = Illumination.classical(alpha=20.0, polarization="unpolarized")

    solver = get_solver(solver_name)
    solver.solve(problem, illumination, wavelengths, **options)  # warm-up
    times = []
    for _ in range(spec["repeat"]):
        started = time.perf_counter()
        scan = solver.solve(problem, illumination, wavelengths, **options)
        times.append(time.perf_counter() - started)
    return {
        "best": min(times),
        "times": times,
        "checksum": float(np.nansum(scan.efficiency)),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cases", nargs="*", default=list(CASES), choices=list(CASES))
    parser.add_argument("--env", nargs="*", default=["default", "single"],
                        choices=["default", "single"])
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--workers", type=int, default=None,
                        help="pass workers=N to the solve (scan cases)")
    parser.add_argument("--src", type=Path, default=None,
                        help="the gratinglab source tree to time (default: this checkout's src/)")
    parser.add_argument("--child", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    if args.child:
        print(json.dumps(_child(json.loads(args.child))))
        return 0

    src = args.src or Path(__file__).resolve().parents[1] / "src"
    print(f"{'case':<20}" + "".join(f"{e:>12}" for e in args.env) + "   checksum")
    with tempfile.TemporaryDirectory() as library:
        for case in args.cases:
            row, checksum = [], None
            for env_name in args.env:
                env = {k: v for k, v in os.environ.items() if k not in _THREAD_VARS}
                if env_name == "single":
                    env.update({k: "1" for k in _THREAD_VARS})
                env["PYTHONPATH"] = str(src)
                spec = {"case": case, "repeat": args.repeat, "library": library}
                if args.workers is not None:
                    spec["extra"] = {"workers": args.workers}
                out = subprocess.run(
                    [sys.executable, __file__, "--child", json.dumps(spec)],
                    env=env, capture_output=True, text=True, check=False,
                )
                if out.returncode != 0:
                    print(out.stderr, file=sys.stderr)
                    row.append("failed")
                    continue
                result = json.loads(out.stdout.strip().splitlines()[-1])
                row.append(f"{result['best']:.3f}s")
                checksum = result["checksum"]
            print(f"{case:<20}" + "".join(f"{r:>12}" for r in row) + f"   {checksum!r}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
