# gratinglab

An open platform for diffraction grating groove metrology and efficiency analysis —
from an AFM scan of a real grating to its diffraction efficiency, without leaving the
package.

Rigorous grating efficiency analysis today means either paying for PCGrate (integral
method) or assembling a patchwork of open RCWA codes that each carry their own
conventions, normalizations, and blind spots. Nothing open lets you run the *same*
physical problem through several methods and compare honestly — and nothing open
connects a *measured* surface to any of them.

`gratinglab` is built around one idea: **the physical problem and the numerical method
are separate, and the problem spec is serializable.** One `Problem` goes to every solver,
and the disagreement between methods is made visible rather than hidden.

```
AFM scan  ->  BoundaryProfile   (gratinglab.metrology)
                      |          measured period, averaged groove, facet fit
                      v
Problem       (period + profile + materials — no solver fields, ever)
Illumination  (wavelength + direction + polarization)
                      |
                      v
Solver (plugin)  ->  Result (+ provenance)
```

The two halves are one package because the seam between them loses information when
it is a file. A PCGrate `.ggp` boundary carries a *shape* and not a grating: it cannot
hold the period, so the period gets retyped, guessed, or recovered from a sidecar
later. `BoundaryProfile.to_problem()` carries the period the scan measured, and the
fitted facet angle with it.

```python
from gratinglab.metrology.boundary import build_boundary_profile

profile = build_boundary_profile(data, scan_x_size, settings)
problem = profile.to_problem(coating="Au", blaze_angle=29.5)
```

## Status

Pre-alpha. The scalar solver, the comparison harness, the convergence harness and the
materials layer work end to end — efficiencies are absolute when a coating is named. RCWA
is contributed and not yet integrated.

## Scope

Two halves. **Metrology** (`gratinglab.metrology`) reads AFM scans — Nanoscope `.spm`
or Gwyddion text exports — and produces blaze angles per groove with row-group
statistics and an ICC correction, plus one averaged groove normalised to a period.
**Efficiency** is everything else: the problem spec, the solvers, the comparison and
convergence harnesses, and the materials layer.

Both mounts are first-class:

- **In-plane / classical** — UV–Vis–NIR optical engineering.
- **Extreme off-plane (conical), grazing incidence** — soft X-ray spectroscopy, where
  PCGrate is the community standard and general-purpose open RCWA codes are weakest.

Planned method backends, in order:

| Method | Status |
|---|---|
| Scalar (Kirchhoff, X-ray blaze form) | **done**, validated |
| RCWA (in-plane, conical, off-plane; Li factorization, S-matrix) | contributed, not yet integrated |
| C-method (Chandezon) | planned |
| Integral method | planned |

An already-computed efficiency table from an external code is a first-class "method" in
the comparison harness, so scalar-vs-RCWA-vs-integral plots are possible without a
licence for the code that produced the reference data:

```python
reference = read_scan("run.txt", method="integral")
align(sweep(problem, illumination, wavelengths, ["scalar", reference]))
```

The caller declares what physics produced an imported file, because the file itself does
not say. Comparison legends then read `scalar vs integral` rather than naming a product,
while `Provenance.version` still records exactly which program and version generated it.

## Install

```bash
pip install -e .                      # solvers only: numpy, scipy, pydantic
pip install -e ".[metrology]"         # + AFM readers and the boundary pipeline
pip install -e ".[gui]"               # + the Qt windows
```

The core stays deliberately light — a solver-only user never installs a plotting
stack — and each subpackage says which extra it needs rather than failing with a bare
`ModuleNotFoundError`.

### Where the measurement data lives

Scans are **not** in the repository. They are research data, and the same rule already
applies to the PCGrate reference corpus:

| Data | Environment variable | Default |
|---|---|---|
| AFM scans | `GRATINGLAB_AFM_DIR` | `~/Documents/afm_scans` |
| PCGrate reference corpus | `GRATINGLAB_REF_DIR` | `~/Documents/diffraction_efficiency` |

Tests that need either **skip cleanly** when it is absent, so the suite is green on a
fresh clone. Most of the metrology suite does not need it at all: it runs against a
committed synthetic scan of a *known* grating
(`tools/metrology/make_synthetic_scan.py`), which is also the control that lets a test
assert the pipeline recovers the period and blaze angle it was given.

## GUI

```bash
.venv/bin/python -m pip install -e ".[gui,metrology]"
gratinglab-gui                 # geometry, solvers, comparison
gratinglab-metrology-gui       # AFM import, blaze angles, boundary profiles
```

Qt (PySide6), so `pip install "gratinglab[gui]"` is self-sufficient — no system toolkit
to install first. PySide6 is LGPLv3 and an optional, dynamically-linked extra, which
leaves this package BSD-3.

Grating geometry lives in a dock on the left (Ctrl+G closes it and hands the width back),
and the tabs are one per concern:

- **Grating Geometry** — a rotatable 3D view of the diffraction cone, with every wave
  vector drawn as an exact unit vector. Because the orders separate in *azimuth* rather
  than polar angle — 93.8° of fan at 0.0° of polar spread for the reference off-plane
  geometry — the "down the cone axis" preset plus a uniform 27× zoom makes them readable
  at true scale. γ is never exaggerated: uniform zoom is a similarity transform, so every
  angle in the view survives it exactly, and the magnification is stated on the canvas.
  Beside it, the dispersion-plane cross-section, a `sin β` ladder carrying evanescent
  orders too, and the groove profile. The picture follows the form as you type.
- **Scalar (Kirchhoff)** — the solve, efficiency per order with per-order visibility, and
  the **provenance panel**: every validity warning verbatim, the summed-efficiency check,
  and whether convergence was actually demonstrated. Nothing is hidden there; a result
  that violates energy conservation says so. **Check convergence…** runs the sweep and
  reports the coarsest defensible setting, with a progress bar and a Cancel that really
  stops it — solvers report per wavelength, so cancellation is not a euphemism for
  looking away.
- **Setup** — the material library: which optical constants are installed, what range
  each covers, where they came from, and how to add your own. The coating *control*
  lives in the geometry dock, because it is a `Problem` field like the period and one
  `Problem` feeds every tab.

Help carries the theory pages, rendered with typeset math and real tables.

A macOS `.app` can be built with `python tools/make_app.py` (macOS 13+, which is the
PySide6 wheel's own floor); the icon is regenerated by `python tools/make_icon.py`.

## Documentation

| Document | What it holds |
|---|---|
| [`docs/roadmap.md`](docs/roadmap.md) | **Where the project is, what is next, and the open questions.** Start here |
| [`docs/conventions.md`](docs/conventions.md) | **Normative.** Signs, units, normalizations, and the reference errata |
| [`docs/theory/scalar.md`](docs/theory/scalar.md) | Full derivation for the scalar solver, its closed forms and its limits |
| [`docs/theory/metrology.md`](docs/theory/metrology.md) | What the AFM pipeline assumes, what it measures, and what it does not |
| [`docs/findings.md`](docs/findings.md) | Empirical results, each with the evidence that established it |
| [`docs/mutation-testing.md`](docs/mutation-testing.md) | What is mutated, what is not, and how to read a survivor |
| [`benchmarks/corpus.toml`](benchmarks/corpus.toml) | Geometry for the reference corpus, most of it recovered from the data |

`conventions.md` is worth singling out: read it before writing a solver or filing a bug
about a sign. The short version is time convention `exp(-iωt)`, lossy *n* = n′ + ik,
grating equation `sin α + sin β_m = mλ/(p sin γ)`, efficiencies **absolute**, lengths in
nm, angles in degrees at the API boundary.

A theory page accompanies each solver as it lands, stating its assumptions, its
derivation in the code's own notation, and the conditions under which it stops being
trustworthy.

## Development

```bash
.venv/bin/python -m pip install -e ".[dev,gui,metrology]"
.venv/bin/python -m pytest
.venv/bin/ruff check .
```

Three practices this project relies on:

- **Physics self-checks over reference data.** `checks.check_reciprocity` and
  `checks.check_energy_balance` constrain the *model* rather than comparing against a
  formula derived the same way the solver computes it. They need no reference data and
  apply to every backend.
- **Convergence demonstrated, not assumed.** `convergence.check_convergence` sweeps
  whatever a solver declares as its accuracy knob and reports the coarsest setting it
  can defend — or reports that it found none, which is a result rather than a missing
  one. It requires *three* consecutive values to agree, because quadrature error is
  measurably non-monotone for a blazed profile and one agreement can be an accident
  the next refinement contradicts
  ([findings](docs/findings.md#quadrature-error-is-not-monotone-in-the-number-of-points)).

  ```python
  report = check_convergence(scalar, problem, illumination, wavelengths)
  report.converged_at           # 4096 -- what a production run should use
  report.scan.provenance.converged   # True, at last
  ```
- **Mutation testing.** Deliberately corrupting the physics and confirming the suite
  notices — coverage says a line ran, this asks whether anything would have noticed if
  it were wrong. It has found real gaps both times it was run; see
  [findings](docs/findings.md#mutation-testing-found-a-real-gap) and
  [docs/mutation-testing.md](docs/mutation-testing.md) for how to run it.

  ```bash
  .venv/bin/python -m pip install -e ".[mut]"
  .venv/bin/mutmut run
  ```

## License

BSD-3-Clause.
