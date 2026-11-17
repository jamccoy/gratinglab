# Examples: using gratinglab from your own code

gratinglab is a Python library, and the windows are built on top of it.
Everything they compute is a call you can make yourself. These scripts are
short, run as they are, and are meant to be copied from.

```bash
pip install gratinglab        # the library alone; no GUI needed
python examples/01_efficiency_scan.py
```

None of them needs an optical-constant table: they use `typical:` materials,
the built-in approximate X-ray model, so they run on a fresh install.

| Script | What it shows |
|---|---|
| [`01_efficiency_scan.py`](01_efficiency_scan.py) | The three objects every calculation starts from (`Problem`, `Illumination`, a solver) and how to read the result. |
| [`02_compare_solvers.py`](02_compare_solvers.py) | The same grating through the scalar solver and the RCWA, and the difference between them. |
| [`03_your_own_material.py`](03_your_own_material.py) | Reading your own `wavelength, n, k` table, keeping a library, and the four ways to name a material. |
| [`04_mirror_and_coating.py`](04_mirror_and_coating.py) | A coating stack as a flat mirror: reflectivity, penetration depth, and when a layer is thick enough to treat as bulk. |
| [`05_parameter_sweep.py`](05_parameter_sweep.py) | Something the windows have no slider for: sweeping the blaze angle and the coating for the best first order. |

## Where things are

What a script usually needs comes from the package itself:

```python
from gratinglab import Problem, Blazed, Illumination, solve, materials, Mirror
```

| You want | Name |
|---|---|
| The grating, and layers on it | `Problem`, `ConformalLayer` |
| Groove shapes | `Blazed`, `Lamellar`, `Sinusoidal`, `Trapezoidal` |
| The light | `Illumination` |
| To solve in one call | `solve(problem, illumination, wavelengths, method="scalar", ...)` |
| A solver object, and their names | `get_solver`, `available_solvers` |
| Several solvers at once | `sweep`, `align` |
| Materials | `materials` |
| A flat mirror | `Mirror` |

Two more are imported from their own modules:

| You want | Import |
|---|---|
| Whether a result has converged | `from gratinglab.convergence import check_convergence` |
| A grating from an AFM scan | `from gratinglab.metrology.boundary import build_boundary_profile` |

The windows can write a script for you: **File, Copy as Python** in the
efficiency window copies the problem on screen and the solve its tab would
run, and in the coatings window it copies the Mirror tab's calculation.
[`docs/python.md`](../docs/python.md) is the longer tour.

Lengths are nanometres and angles are degrees throughout. Conventions for
angles, orders and signs are in [`docs/conventions.md`](../docs/conventions.md),
and each solver's theory page under [`docs/theory/`](../docs/theory/) says what
it assumes.

The interface is not yet declared stable: the project is before 1.0, and a
minor release may change it, with the change stated in the changelog. These
scripts are run by the test suite (`tests/test_examples.py`), so they match the
version they ship with.

`metrology/uncertainty_demo.py` is older and stands apart: it illustrates how
the blaze-angle uncertainty is decomposed, and does not call the package.
