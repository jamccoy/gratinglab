# Using gratinglab from Python

gratinglab is a library, and its windows are built on it. Anything a window
computes is a call you can make from a script, a notebook or a program of your
own. This page is the tour; the runnable versions are in
[`examples/`](../examples/README.md).

```bash
pip install gratinglab
```

That installs the library with numpy, scipy and pydantic, and no GUI toolkit.
Lengths are nanometres and angles are degrees throughout;
[`conventions.md`](conventions.md) defines the signs and the orders.

## The quickest way in: let a window write it

Set a problem up in the efficiency window and choose **File, Copy as Python**
(`Ctrl+Shift+C`). The clipboard then holds a script that rebuilds the problem
on screen and solves it as the showing tab would, with that tab's settings. On
the Compare tab it compares the ticked solvers; on the Setup and Grating
Geometry tabs it stops at the problem. In the coatings window the same menu
item copies the Mirror tab's calculation.

The script is written from what the solver was given, not from the text in
the form, and it is a starting point: paste it, run it, then change it.

## The three objects

Every efficiency calculation is a grating, the light, and a solver.

```python
import numpy as np
from gratinglab import Problem, Blazed, Illumination, solve

problem = Problem(period=160.0, profile=Blazed(blaze_angle=30.0), coating="typical:Au")
light = Illumination.offplane(graze=1.5, azimuth=30.0)
scan = solve(problem, light, np.linspace(1.0, 5.0, 41), method="scalar")
```

**`Problem`** is the grating: its `period`, its groove `profile`, the material
it is cut into (`coating`), an RMS `roughness`, and any `layers` on top
(`ConformalLayer(material, thickness)`, innermost first). It is immutable;
`problem.model_copy(update={"period": 200.0})` makes a changed one, which is
how a sweep is written.

**Profiles** are `Blazed`, `Lamellar`, `Sinusoidal` and `Trapezoidal`. A
measured groove comes from a file, with
`gratinglab.io.boundary.read_profile_file(path)`, or from an AFM scan through
`gratinglab.metrology`.

**`Illumination`** has a constructor per mount: `Illumination.offplane(graze,
azimuth)` for the extreme off-plane mount of soft X-ray work,
`Illumination.classical(alpha)` in plane, and `Illumination.conical(theta,
phi)` in the angles RCWA codes use.

**Solvers** are asked for by name: `"scalar"` is fast and approximate, and
`"rcwa"`, `"integral"` and `"cmethod"` are rigorous. `solve(...)` is one call;
`get_solver(name)` gives the solver object, with `.capabilities` saying what
it can model and `.solve(...)` to call repeatedly. Each solver's own settings
are keyword arguments, such as `quadrature_points` for the scalar solver and
`truncation` for the RCWA; its theory page under [`theory/`](theory/) says what
they mean. A solver that cannot model a problem refuses it with
`UnsupportedConfiguration` rather than solving something nearby.

## The result

`solve` returns an `EfficiencyScan`: a table with one row per wavelength and
one column per order.

```python
scan.order(1)            # order +1 across the scan, an array
scan.at(2.0)[2]          # order +2 at 2.0 nm
scan.at(2.0).total       # everything reflected into orders at 2.0 nm
scan.absorption          # the absorbed fraction, from solvers that compute it
scan.to_records()        # plain dicts, for pandas or csv
scan.provenance          # which solver, which settings, which material table
```

An order that does not propagate at a wavelength holds 0.0 and is flagged in
`scan.propagating`; it is never dropped, so a column always means one order.

## More than one solver

```python
from gratinglab import sweep, align

comparison = align(sweep(problem, light, wavelengths, ["scalar", "rcwa"]))
comparison.summary("scalar", "rcwa")   # the largest and RMS difference, and where
```

`gratinglab.convergence.check_convergence(solver, problem, light, wavelengths)`
raises a solver's accuracy setting until the answer stops moving, and reports
where it did.

## Materials

A material reaches a calculation by name, as `coating`, a layer's `material`,
or a mirror's `substrate`.

| Name | What it is |
|---|---|
| `"Au"` | A table in your library: `~/.gratinglab/materials`, or the folders in `GRATINGLAB_MATERIALS` |
| `"typical:Au"`, `"typical:C@2.0"` | The built-in approximate X-ray model, with an optional density in g/cm³ |
| `"factors:SiO2@2.2"` | Tabulated atomic scattering factors (`Si.nff`, `O.nff` in your library), a formula and a density |
| `"Au@17.5"` | A CXRO table from your library, rescaled to another density |

```python
from gratinglab import materials

gold = materials.lookup("typical:Au")
gold.n(1.0)                              # complex index at 1 nm
materials.from_nk_file("film.csv")       # your own wavelength, n, k columns
materials.available()                    # what is in your library
```

No tables ship with the package. Asking for a wavelength outside a table's
range raises an error; nothing is extrapolated.

## A coating as a mirror

```python
from gratinglab import Mirror, ConformalLayer

mirror = Mirror("typical:Si", layers=(ConformalLayer(material="typical:Au", thickness=20.0),))
mirror.reflectivity(graze_deg=[1.0, 2.0, 3.0], wavelength_nm=1.0)
mirror.penetration_depth(2.0, 1.0)       # 1/e depth in each material, nm
mirror.bulk_thickness(2.0, 1.0)          # gold needed before the silicon stops mattering
```

[`theory/mirror.md`](theory/mirror.md) says what each of these means.

## What may change

The project is before 1.0 and the interface is not declared stable. A minor
release may change it, and says so in the
[changelog](../CHANGELOG.md). The names on this page are the ones intended to
last, and the scripts in `examples/` are run by the test suite, so they match
the version they ship with.
