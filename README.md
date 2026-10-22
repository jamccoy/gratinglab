# gratinglab

[![PyPI](https://img.shields.io/pypi/v/gratinglab?cacheSeconds=3600)](https://pypi.org/project/gratinglab/)
[![Python](https://img.shields.io/pypi/pyversions/gratinglab?cacheSeconds=3600)](https://pypi.org/project/gratinglab/)
[![tests](https://github.com/jamccoy/gratinglab/actions/workflows/test.yml/badge.svg)](https://github.com/jamccoy/gratinglab/actions/workflows/test.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue)](https://github.com/jamccoy/gratinglab/blob/master/LICENSE)

**[PyPI](https://pypi.org/project/gratinglab/)** · **[Documentation](https://github.com/jamccoy/gratinglab/wiki)** · **[Releases](https://github.com/jamccoy/gratinglab/releases)** · **[Issues](https://github.com/jamccoy/gratinglab/issues)**

An open platform for diffraction grating groove metrology and efficiency analysis,
from an AFM scan of a real grating to its diffraction efficiency, in one package.

Rigorous grating efficiency analysis today means either paying for PCGrate (integral
method) or assembling a patchwork of open RCWA codes that each carry their own
conventions, normalizations, and blind spots. No open code runs the *same* physical
problem through several methods for a fair comparison, and none connects a *measured*
surface to any of them.

`gratinglab` is built around one idea: **the physical problem and the numerical method
are separate, and the problem spec is serializable.** One `Problem` goes to every solver,
so disagreement between methods shows up in the comparison.

```
AFM scan  ->  BoundaryProfile   (gratinglab.metrology)
                      |          measured period, averaged groove, facet fit
                      v
Problem       (period + profile + materials; no solver fields)
Illumination  (wavelength + direction + polarization)
                      |
                      v
Solver (plugin)  ->  Result (+ provenance)
```

The two halves are one package because a file between them loses information. A
PCGrate `.ggp` boundary carries a *shape* and not a grating: it cannot hold the period,
so the period gets retyped, guessed, or recovered from a sidecar later.
`BoundaryProfile.to_problem()` carries the period the scan measured, and the fitted
facet angle with it. To save a profile without losing that, write it as gratinglab's
own `.json` (`gratinglab.io.boundary.write_boundary_json`): reading it back gives the
same problem, and the efficiency window's "Load boundary profile" fills in the
period from it. The metrology window's Boundary tab exports `.ggp`, `.json` or
`.csv`, and its "Send to efficiency window" button opens the groove there directly.

```python
from gratinglab.metrology.boundary import build_boundary_profile

profile = build_boundary_profile(data, scan_x_size, settings)
problem = profile.to_problem(coating="Au", blaze_angle=29.5)
```

## Status

Pre-alpha. The scalar and integral solvers, the comparison harness, the convergence
harness and the materials layer work end to end. Efficiencies are absolute when a
coating is named, using an optical constant table of your own (see Install). RCWA is
contributed and not yet integrated.

## Scope

Two halves. **Metrology** (`gratinglab.metrology`) reads AFM scans (Nanoscope `.spm`
or Gwyddion text exports) and produces blaze angles per groove with row-group
statistics and an ICC correction, plus one averaged groove normalised to a period.
Tip convolution can be undone where the tip reached the surface, and the region it
never reached is marked rather than smoothed over. **Efficiency** is everything else:
the problem spec, the solvers, the comparison and convergence harnesses, the materials
layer, and the resolving power the measured groove count implies.

Both mounts are supported:

- **In-plane / classical:** UV–Vis–NIR optical engineering.
- **Extreme off-plane (conical), grazing incidence:** soft X-ray spectroscopy, where
  PCGrate is the community standard and general-purpose open RCWA codes are weakest.

Planned method backends, in order:

| Method | Status |
|---|---|
| Scalar (Kirchhoff, X-ray blaze form) | **done**, validated |
| Integral method (boundary integral, conical, undercut-capable) | **done**: perfect conductivity validated orderwise against PCGrate, and finite conductivity (`conductivity="tabulated"`, Goray & Schmidt 2010) with absolute efficiencies and recorded absorption |
| RCWA (in-plane, conical, off-plane; Li factorization, S-matrix) | contributed, not yet integrated |
| C-method (Chandezon) | planned |

The comparison harness accepts an efficiency table already computed by an external code
as one of its methods, so scalar-vs-RCWA-vs-integral plots need no licence for the code
that produced the reference data:

```python
reference = read_scan("run.txt", method="integral")
align(sweep(problem, illumination, wavelengths, ["scalar", reference]))
```

The caller declares what physics produced an imported file, because the file itself does
not say. Comparison legends then read `scalar vs integral` rather than naming a product,
while `Provenance.version` still records which program and version generated it.

### From a measured surface to a spectrograph number

The scan carries its own groove count, so resolving power needs nothing retyped:

```python
problem = profile.to_problem(coating="Au")     # n_grooves = the grooves averaged
resolving_power(problem, illumination, wavelength=2.0, order=3).r
```

`R = |m|·N` by the Rayleigh criterion, with `line_profile` giving the finite-*N* line
shape it is derived from. A problem whose `n_grooves` is unset is **refused** rather
than treated as infinite: "R is undefined" and "R is very large" are different answers.
Resolving power degraded by real groove-placement error is planned next, and the ideal
formula is not used to approximate it.

### The tip is not the surface

An AFM image is the surface dilated by the tip. `tip_correction='erosion'` inverts that
where the inversion exists (a cone-plus-spherical-apex tip, Villarrubia 1997) and
reports a **certainty map** for where it does not. No algorithm can recover a trough
the tip never entered, and a correction that filled it in silently would be worse than
no correction.

It has already been useful. On the one real scan whose depth and facet angle disagree
by 7.6°, erosion with the nominal probe recovers *nothing*, and reproducing the deficit
forward takes an apex worn to about 80 nm, which rules out the nominal tip as the
explanation ([findings](https://github.com/jamccoy/gratinglab/blob/master/docs/findings.md)).

## Install

```bash
pip install gratinglab                       # solvers only: numpy, scipy, pydantic
pip install "gratinglab[metrology]"          # + AFM readers and the boundary pipeline
pip install "gratinglab[gui]"                # + the Qt windows
pip install "gratinglab[refractiveindex]"    # + refractiveindex.info tables
```

With [uv](https://docs.astral.sh/uv/), the same extras work:

```bash
uv pip install "gratinglab[refractiveindex]"    # into the active environment
uv add gratinglab                               # into a uv project
uv tool install "gratinglab[gui]"               # gratinglab-gui as a standalone command
```

From a checkout, `pip install -e ".[dev,gui,metrology,refractiveindex]"` instead.

No optical constant tables ship with the package, so absolute efficiencies need
your own. A library is a directory of them, one material per file, named for
the material:

| File | From | Range |
|---|---|---|
| `Pt.txt` | a [CXRO](https://henke.lbl.gov/optical_constants/) export | X-ray |
| `Pt.yml` | a [refractiveindex.info](https://refractiveindex.info) file, renamed from its author's name | ultraviolet to infrared |
| `Pt.ari` | a PCGrate table | whatever it covers |

Save one as `~/.gratinglab/materials/Pt.yml` and `coating="Pt"` works.
`GRATINGLAB_MATERIALS` points elsewhere, and `materials.from_cxro_file(path)`
or `from_refractiveindex_file(path)` reads a file without a library. The YAML
files need the `refractiveindex` extra. See
[`SOURCES.md`](https://github.com/jamccoy/gratinglab/blob/master/src/gratinglab/materials/data/SOURCES.md)
for why none ship.

The core install has no plotting stack, so a user who only runs solvers never installs
one. Each subpackage that needs an extra names it when it is missing, rather than
failing with a bare `ModuleNotFoundError`.

Besides numpy, scipy and pydantic, the one dependency is
[endstation](https://github.com/jamccoy/endstation), shared scaffolding for desktop
science tools by the same author. gratinglab uses its optional-extra guards (the
messages above), form validation, Qt layout sizing, test helpers, and the tools that
build the icon and the macOS app and publish the wiki. endstation's own core has no
dependencies, so it adds nothing to a solver-only install; its Qt parts arrive only
with the `gui` extra.

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
(`tools/metrology/make_synthetic_scan.py`), which also lets a test check that the
pipeline recovers the period and blaze angle it was given.

## GUI

```bash
.venv/bin/python -m pip install -e ".[gui,metrology]"
gratinglab-gui                 # geometry, solvers, comparison
gratinglab-metrology-gui       # AFM import, blaze angles, boundary profiles
```

Qt (PySide6), so `pip install "gratinglab[gui]"` is self-sufficient, with no system
toolkit to install first. PySide6 is LGPLv3 and an optional, dynamically-linked extra,
which leaves this package BSD-3.

Grating geometry lives in a dock on the left (Ctrl+G closes it and hands the width back),
and there is one tab per concern:

- **Grating Geometry:** a rotatable 3D view of the diffraction cone, with every wave
  vector drawn as an exact unit vector. The orders separate in *azimuth* rather than
  polar angle (93.8° of fan at 0.0° of polar spread for the reference off-plane
  geometry), so the "down the cone axis" preset plus a uniform 27× zoom makes them
  readable at true scale. γ is never exaggerated: uniform zoom is a similarity
  transform, so every angle in the view is preserved exactly, and the magnification is
  stated on the canvas. Beside it are the dispersion-plane cross-section, a `sin β`
  ladder that includes evanescent orders, and the groove profile. The picture updates
  as you type.
- **Scalar (Kirchhoff):** the solve, efficiency per order with per-order visibility, and
  the **provenance panel**: every validity warning verbatim, the summed-efficiency check,
  and whether convergence was demonstrated. A result that violates energy conservation
  says so. **Check convergence…** runs the sweep and reports the coarsest defensible
  setting, with a progress bar and a Cancel button. Solvers report per wavelength, so
  Cancel stops the computation, not only the display.
- **Integral method:** the rigorous backend, with the same panel and the same
  provenance, plus its own boundary condition control (perfect conductivity or
  tabulated optical constants). It takes minutes where scalar takes milliseconds, so the
  window opens on the fast solver and this tab solves when you ask it to.
- **Setup:** the material library: which optical constants are installed, what range
  each covers, where they came from, and how to add your own. The coating *control*
  lives in the geometry dock, because it is a `Problem` field like the period and one
  `Problem` feeds every tab.

Each registered solver gets a tab generated from the registry, so a contributed backend
gets its panel by registering rather than by editing the window.

Help carries the theory pages, rendered with typeset math and real tables.

A macOS `.app` can be built with `python tools/make_app.py` (macOS 13+, which is the
PySide6 wheel's own floor); the icon is regenerated by `python tools/make_icon.py`.

## Documentation

| Document | What it holds |
|---|---|
| [`docs/roadmap.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/roadmap.md) | **Where the project is, what is next, and the open questions.** Start here |
| [`docs/conventions.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/conventions.md) | **Normative.** Signs, units, normalizations, and the reference errata |
| [`docs/theory/scalar.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/theory/scalar.md) | Full derivation for the scalar solver, its closed forms and its limits |
| [`docs/theory/integral.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/theory/integral.md) | The boundary-integral solver: conical formulation, perfect and finite conductivity |
| [`docs/theory/metrology.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/theory/metrology.md) | What the AFM pipeline assumes, what it measures, and what it does not, including what tip erosion can and cannot give back |
| [`docs/findings.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/findings.md) | Empirical results, each with the evidence that established it |
| [`docs/references.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/references.md) | The literature each part is built from, mapped to the part it serves |
| [`docs/mutation-testing.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/mutation-testing.md) | What is mutated, what is not, and how to read a survivor |
| [`benchmarks/corpus.toml`](https://github.com/jamccoy/gratinglab/blob/master/benchmarks/corpus.toml) | Geometry for the reference corpus, most of it recovered from the data |

Read `conventions.md` before writing a solver or filing a bug about a sign. The short
version is time convention `exp(-iωt)`, lossy *n* = n′ + ik, grating equation
`sin α + sin β_m = mλ/(p sin γ)`, efficiencies **absolute**, lengths in nm, angles in
degrees at the API boundary.

Each solver has a theory page stating its assumptions, its derivation in the code's own
notation, and the conditions under which it stops being trustworthy.

The [wiki](https://github.com/jamccoy/gratinglab/wiki) is a published copy of `docs/`,
not a second source. Edit `docs/` and republish, or the copy people read becomes the
one nobody edits:

```bash
git clone https://github.com/jamccoy/gratinglab.wiki.git ../gratinglab.wiki
.venv/bin/python tools/publish_wiki.py --wiki ../gratinglab.wiki   # then review and push
```

## Development

```bash
.venv/bin/python -m pip install -e ".[dev,gui,metrology,refractiveindex]"
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
  can defend, or reports that it found none, which is itself a result. It requires
  *three* consecutive values to agree, because quadrature error is measurably
  non-monotone for a blazed profile and one agreement can be an accident the next
  refinement contradicts
  ([findings](https://github.com/jamccoy/gratinglab/blob/master/docs/findings.md#quadrature-error-is-not-monotone-in-the-number-of-points)).

  ```python
  report = check_convergence(scalar, problem, illumination, wavelengths)
  report.converged_at           # 4096 -- what a production run should use
  report.scan.provenance.converged   # True
  ```
- **Mutation testing.** The suite is run against deliberately broken copies of the
  physics to check that it notices. Coverage says a line ran; this asks whether any test
  would fail if the line were wrong. It has found real gaps both times it was run; see
  [findings](https://github.com/jamccoy/gratinglab/blob/master/docs/findings.md#mutation-testing-found-a-real-gap)
  and [docs/mutation-testing.md](https://github.com/jamccoy/gratinglab/blob/master/docs/mutation-testing.md)
  for how to run it.

  ```bash
  .venv/bin/python -m pip install -e ".[mut]"
  .venv/bin/mutmut run
  ```

## License

BSD-3-Clause.
