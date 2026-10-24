# gratinglab

[![PyPI](https://img.shields.io/pypi/v/gratinglab?cacheSeconds=3600)](https://pypi.org/project/gratinglab/)
[![Python](https://img.shields.io/pypi/pyversions/gratinglab?cacheSeconds=3600)](https://pypi.org/project/gratinglab/)
[![tests](https://github.com/jamccoy/gratinglab/actions/workflows/test.yml/badge.svg)](https://github.com/jamccoy/gratinglab/actions/workflows/test.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue)](https://github.com/jamccoy/gratinglab/blob/master/LICENSE)

**[PyPI](https://pypi.org/project/gratinglab/)** · **[Documentation](https://github.com/jamccoy/gratinglab/wiki)** · **[Changelog](https://github.com/jamccoy/gratinglab/blob/master/CHANGELOG.md)** · **[Releases](https://github.com/jamccoy/gratinglab/releases)** · **[Issues](https://github.com/jamccoy/gratinglab/issues)**

gratinglab is a Python package for measuring the groove profiles of diffraction
gratings from atomic force microscope (AFM) scans and for computing the diffraction
efficiency of a grating from such a profile.

The AFM metrology is the part of the package in regular use. It measures blaze angles
from a scan, averages the grooves into a single profile, reports how much that profile
varies down the scan, corrects for the shape of the AFM tip where the data allow it,
and exports the averaged groove as a boundary profile. The efficiency solvers are
**experimental**: they pass their test suite but have not yet been used for real
analysis, and their results should be treated as provisional.

Rigorous methods for computing grating efficiency are well established, and open codes
exist for several of them, although each carries its own conventions and
normalizations. gratinglab is intended to fill two gaps that remain: running the same
physical problem through several methods so that their results can be compared
directly, and starting that problem from a measured surface rather than an idealised
one. To this end, the physical problem and the numerical method are kept separate. A
`Problem` describes the grating and can be saved to a file, and the same `Problem` is
passed to every solver, so any difference between their results is attributable to
the methods themselves.

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

The metrology and the solvers are distributed as one package so that no information is
lost between them. A plain boundary file records the groove shape but not the grating
period, which must then be entered again by hand. `BoundaryProfile.to_problem()`
carries the period measured from the scan, together with the fitted facet angle. A
profile can be saved with that information intact in gratinglab's own `.json` format
(`gratinglab.io.boundary.write_boundary_json`); reading it back gives the same
problem, and "Load boundary profile" in the efficiency window takes the period from
it. The Boundary tab of the metrology window exports `.json`, `.csv`, or a `.ggp`
boundary file for PCGrate, and its "Send to efficiency window" button opens the groove
in the efficiency window directly.

```python
from gratinglab.metrology.boundary import build_boundary_profile

profile = build_boundary_profile(data, scan_x_size, settings)
problem = profile.to_problem(coating="Au", blaze_angle=29.5)
```

## Status

The project is at a pre-alpha stage. The AFM metrology is the more mature of its two
parts. On the efficiency side, the scalar and integral solvers, the comparison and
convergence harnesses and the materials layer work end to end. They are checked
against closed-form results, published benchmark tables, energy balance and
reciprocity, but remain experimental until they have been used on real problems.
Efficiencies are absolute when a coating is named, using an optical constant table
supplied by the user (see Install). A rigorous coupled-wave analysis (RCWA) backend
is planned, written in-house.

## Scope

The metrology (`gratinglab.metrology`) reads AFM scans, either Nanoscope `.spm` files
or Gwyddion text exports. It produces blaze angles for each groove with statistics
over groups of scan rows and a correction for the correlation between those groups,
an averaged groove normalised to one period, and statistics describing how that
groove varies between horizontal swaths of the scan. It can also undo the broadening
introduced by the AFM tip where the tip reached the surface, and it marks the regions
the tip could not reach instead of estimating them.

The efficiency side comprises the problem specification, the solvers, the comparison
and convergence harnesses, the materials layer, and a calculation of the resolving
power implied by the number of grooves measured.

The aim is gratings in general: reflection and transmission, metallic and
dielectric, bare and coated, from X-rays to the infrared. Every mount is
represented, from in-plane (classical) through conical to extreme off-plane at
grazing incidence, the soft X-ray spectroscopy regime the project began with, where
the integral method is the established approach and general-purpose RCWA codes are
least reliable. Which method serves which regime, and what is not yet covered
(coated and multilayer gratings, for one), is mapped in
[`docs/roadmap.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/roadmap.md#regimes-and-the-methods-that-serve-them).

The method backends, in order of development:

| Method | Status |
|---|---|
| Scalar (Kirchhoff, x-ray blaze form) | experimental; checked against closed-form results |
| Integral method (boundary integral, conical, undercut-capable) | experimental; perfect and finite conductivity (`conductivity="tabulated"`, Goray & Schmidt 2010), absolute efficiencies with the absorbed fraction recorded; a corner-graded mesh for profiles with corners (`corner_grading=4`), under either boundary condition; checked against published benchmark tables (Maystre & Popov; Goray & Schmidt), energy balance and reciprocity |
| RCWA (in-plane, conical, off-plane; Li factorization, S-matrix, layer stacks) | planned, written in-house |
| C-method (Chandezon) | planned |

The comparison harness also accepts an efficiency table already computed by an
external program as one of its methods, so that comparison plots can be made without
a licence for the program that produced the reference data:

```python
reference = read_scan("run.txt", method="integral")
align(sweep(problem, illumination, wavelengths, ["scalar", reference]))
```

Because the file does not record the physics behind it, the caller states which method
produced an imported table. Comparison legends then name the method, for example
`scalar vs integral`, while `Provenance.version` records the program and version that
generated the file.

### Resolving power from a measured surface

A scan records the number of grooves it covers, so resolving power can be computed
without entering anything by hand:

```python
problem = profile.to_problem(coating="Au")     # n_grooves = the grooves averaged
resolving_power(problem, illumination, wavelength=2.0, order=3).r
```

The resolving power follows from the Rayleigh criterion as `R = |m|·N`, and
`line_profile` gives the line shape for a finite number of grooves from which that
result is derived. A problem whose `n_grooves` is not set is refused rather than
treated as having infinitely many grooves. Resolving power degraded by errors in groove
placement is planned; until it is implemented, the ideal formula is not used as an
approximation for it.

### Correcting for the AFM tip

An AFM image is the true surface dilated by the shape of the tip.
`tip_correction='erosion'` reverses this dilation where it can be reversed (Villarrubia
1997) and reports a certainty map showing where it cannot. The tip is modelled as a
spherical apex on two flanks whose angles may differ, so a pyramidal probe such as
Bruker's ScanAsyst-Air can be represented. No method can recover the bottom of a
trough the tip did not enter, and the correction leaves such regions as upper bounds
instead of filling them in.

On one real scan, the groove depth and the fitted facet angle disagree by 7.6°.
Erosion with the probe at its nominal radius recovers none of that difference, and
reproducing it by simulation requires an apex worn to about 80 nm, which excludes the
nominal tip as the cause
([findings](https://github.com/jamccoy/gratinglab/blob/master/docs/findings.md)).

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

From a checkout, use `pip install -e ".[dev,gui,metrology,refractiveindex]"` instead.

No optical constant tables are distributed with the package, so absolute efficiencies
require tables supplied by the user. A library is a directory of such tables, one
material per file, each named after its material:

| File | From | Range |
|---|---|---|
| `Pt.txt` | a [CXRO](https://henke.lbl.gov/optical_constants/) export | x-ray |
| `Pt.yml` | a [refractiveindex.info](https://refractiveindex.info) file, renamed from its author's name | ultraviolet to infrared |
| `Pt.ari` | a PCGrate table | whatever it covers |

Saving one as `~/.gratinglab/materials/Pt.yml` makes `coating="Pt"` available.
`GRATINGLAB_MATERIALS` points to a different directory, and
`materials.from_cxro_file(path)` or `from_refractiveindex_file(path)` reads a single
file without a library. The YAML files require the `refractiveindex` extra.
[`SOURCES.md`](https://github.com/jamccoy/gratinglab/blob/master/src/gratinglab/materials/data/SOURCES.md)
explains why no tables are distributed.

The core install does not include a plotting stack, so a user who only runs the solvers
does not install one. A subpackage that needs an optional extra names it when it is
missing, instead of failing with a bare `ModuleNotFoundError`.

Apart from numpy, scipy and pydantic, the only dependency is
**[endstation](https://pypi.org/project/endstation/)**
([source](https://github.com/jamccoy/endstation)), a library of shared components for
desktop scientific applications by the same author, which gratinglab is built on.
gratinglab uses its guards for optional extras (the messages described above), form
validation, Qt layout sizing, test helpers, and the tools that build the icon and the
macOS application and publish the wiki. The core of endstation has no dependencies of
its own, so it adds nothing to a solver-only install; its Qt components are installed
only with the `gui` extra.

### Where the measurement data lives

AFM scans are not included in the repository, because they are research data. The
same applies to the private reference corpus of efficiency tables:

| Data | Environment variable | Default |
|---|---|---|
| AFM scans | `GRATINGLAB_AFM_DIR` | `~/Documents/afm_scans` |
| Reference efficiency tables | `GRATINGLAB_REF_DIR` | `~/Documents/diffraction_efficiency` |

Tests that need either are skipped when it is absent, so the suite passes on a fresh
clone. Most of the metrology tests do not need it: they run against a committed
synthetic scan of a grating with known geometry
(`tools/metrology/make_synthetic_scan.py`), which also allows a test to check that the
pipeline recovers the period and blaze angle used to generate it.

## GUI

```bash
.venv/bin/python -m pip install -e ".[gui,metrology]"
gratinglab-app                 # the chooser: one place to start, both windows
gratinglab-gui                 # straight to geometry, solvers, comparison
gratinglab-metrology-gui       # straight to AFM import, blaze angles, boundary profiles
```

There are two windows, and `gratinglab-app` is the front door that asks which one
you want. Both halves are this one package — the metrology measures a groove and the
solvers compute an efficiency from the result — but the desktop layer used to ship
them as two commands and two application bundles, the AFM one still named after the
repository it was merged from, which made the half in regular use the half that was
hardest to find. The two direct commands above are unchanged, for a shortcut or a
habit.

The windows are built with Qt (PySide6), so `pip install "gratinglab[gui]"` needs no
system toolkit to be installed first. PySide6 is licensed under the LGPLv3 and is an
optional, dynamically linked extra, so this package remains under the BSD-3 licence.

The grating geometry is set in a dock on the left of the efficiency window (Ctrl+G
hides it), and the window has one tab for each task:

- **Grating Geometry:** a rotatable 3D view of the diffraction cone, in which every wave
  vector is drawn as an exact unit vector. In the reference off-plane geometry the
  orders separate in azimuth rather than in polar angle (93.8° of azimuthal spread at
  0.0° of polar spread), so the view looking down the cone axis, combined with a uniform
  27× zoom, makes them distinguishable at true scale. The grazing angle γ is not
  exaggerated: a uniform zoom preserves every angle in the view, and the magnification
  is stated on the canvas. Alongside it are the cross-section in the dispersion plane, a
  `sin β` diagram that includes the evanescent orders, and the groove profile. The view
  updates as values are typed.
- **Scalar (Kirchhoff):** the scalar solution, the efficiency of each order with
  per-order visibility, and a provenance panel that lists every validity warning, the
  summed-efficiency check, and whether convergence was demonstrated. A result that
  violates energy conservation is flagged. "Check convergence…" runs a convergence
  sweep and reports the coarsest setting that meets the criterion, with a progress bar
  and a Cancel button. Because solvers report progress for each wavelength, Cancel
  stops the computation itself.
- **Integral method:** the rigorous backend, with the same panel and provenance, its
  own boundary-condition control (perfect conductivity or tabulated optical
  constants) and a **Corner mesh** control: off, or a grading order that clusters
  the boundary nodes at the profile's corners, which a cornered profile needs to
  converge quickly (4 is recommended). It takes minutes where the scalar solver takes milliseconds, so the
  window opens on the scalar solver and this tab computes only when asked to.
- **Setup:** the material library, showing which optical constants are installed, the
  range each covers, where each came from, and how to add others. The coating itself is
  chosen in the geometry dock, because it is a `Problem` field, like the period, and the
  same `Problem` is used by every tab.

Each registered solver receives a tab generated from the solver registry, so a
new backend appears in the window once it is registered.

The Help menu contains the theory pages, rendered with typeset equations and tables.
**View → Theme** switches between System (following the operating system), Light,
Dark, Solarized Light, Solarized Dark and High contrast; the choice is remembered and
applies to the metrology window too. Diffraction orders are coloured by sign, cool for
negative and warm for positive, with each order's colour shown beside it in the order
list.

A macOS `.app` is built with `python tools/apps.py bundle --app gratinglab`, on macOS
13 or later, the minimum supported by the PySide6 wheels. That is one `GratingLab.app`
opening the chooser, which is what you want installed; `--app metrology` and `--app
efficiency` build bundles that skip the chooser and open one window directly, for
anyone who would rather have that in the Dock. `--dest ~/Applications` puts it where Finder will list it. The bundle
is a stub that execs the console script in this checkout's virtualenv rather than a
frozen application, so it is not portable to another machine; `python tools/apps.py
check --app gratinglab` reports when that virtualenv has moved, which is the one
failure a bundle cannot show from the outside. `python tools/apps.py icon --app
gratinglab` regenerates the icon from the artwork in `gratinglab.gui.icon`.

## Documentation

| Document | Contents |
|---|---|
| [`docs/roadmap.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/roadmap.md) | The current state of the project, planned work and open questions; the place to start |
| [`docs/conventions.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/conventions.md) | The normative definitions: signs, units, normalizations, and errata in the references |
| [`docs/theory/scalar.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/theory/scalar.md) | The derivation of the scalar solver, its closed-form results and its limits |
| [`docs/theory/integral.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/theory/integral.md) | The boundary-integral solver: conical formulation, perfect and finite conductivity |
| [`docs/theory/metrology.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/theory/metrology.md) | The assumptions of the AFM pipeline, what it measures and what it does not, including what tip erosion can recover |
| [`docs/findings.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/findings.md) | Empirical results, each with the evidence for it |
| [`docs/references.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/references.md) | The literature each part of the package is based on |
| [`docs/mutation-testing.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/mutation-testing.md) | What is mutated, what is not, and how to interpret a surviving mutant |
| [`benchmarks/corpus.toml`](https://github.com/jamccoy/gratinglab/blob/master/benchmarks/corpus.toml) | Geometry for the reference corpus, most of it recovered from the data |

`conventions.md` should be read before writing a solver or reporting a bug involving a
sign. In brief: the time convention is `exp(-iωt)`, the complex refractive index is
*n* = n′ + ik, the grating equation is `sin α + sin β_m = mλ/(p sin γ)`, efficiencies
are absolute, lengths are in nm, and angles are in degrees at the API boundary.

Each solver has a theory page stating its assumptions, its derivation in the notation
of the code, and the conditions under which its results cease to be reliable.

The [wiki](https://github.com/jamccoy/gratinglab/wiki) is a published copy of `docs/`
and of the metrology window's own help pages, which appear there with an "AFM" prefix.
Changes are made to the sources in the repository and then republished:

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

The test suite relies on three practices in particular:

- **Physical self-checks.** `checks.check_reciprocity` and `checks.check_energy_balance`
  test properties that any correct model must satisfy, independently of how the
  solver computes its result. They require no reference data and apply to every
  backend.
- **Demonstrated convergence.** `convergence.check_convergence` varies whatever a
  solver declares as its accuracy parameter and reports the coarsest setting that
  meets the criterion, or reports that none did. It requires three consecutive values
  to agree, because the quadrature error of a blazed profile does not decrease
  monotonically, so a single agreement can be coincidental
  ([findings](https://github.com/jamccoy/gratinglab/blob/master/docs/findings.md#quadrature-error-is-not-monotone-in-the-number-of-points)).

  ```python
  report = check_convergence(scalar, problem, illumination, wavelengths)
  report.converged_at           # 4096 -- what a production run should use
  report.scan.provenance.converged   # True
  ```
- **Mutation testing.** The suite is run against deliberately altered copies of the
  physics code to check that the changes are detected. Coverage shows that a line was
  executed; mutation testing shows whether any test would fail if that line were wrong.
  Both runs so far have found real gaps in the tests; see
  [findings](https://github.com/jamccoy/gratinglab/blob/master/docs/findings.md#mutation-testing-found-a-real-gap)
  and [docs/mutation-testing.md](https://github.com/jamccoy/gratinglab/blob/master/docs/mutation-testing.md)
  for how to run it.

  ```bash
  .venv/bin/python -m pip install -e ".[mut]"
  .venv/bin/mutmut run
  ```

## License, copyright and reuse

gratinglab is released under the
[BSD 3-Clause License](https://github.com/jamccoy/gratinglab/blob/master/LICENSE),
copyright (c) 2026 Jake McCoy. What follows is a plain summary for convenience; the
`LICENSE` file is what governs, and nothing here is legal advice.

**What you may do.** Use, modify and redistribute the code, in source or binary form,
for any purpose, including commercial use and inside closed-source work. There is no
copyleft obligation: software that imports gratinglab does not become BSD-licensed,
and you are not required to publish your modifications.

**What you must do.** Retain the copyright notice, the list of conditions and the
disclaimer when you redistribute — in the source for a source distribution, and in the
documentation or other accompanying materials for a binary one. The third clause
separately forbids using the author's or contributors' names to endorse or promote a
derived product without prior written permission.

**What you do not get.** The software is provided as is, without warranty of any kind,
and the author is not liable for any use of it. That disclaimer carries more weight
here than it does for most packages, because the efficiency solvers are experimental
(see Status): a number this software produces is evidence to be checked against
another method, not a result to be relied on because a computer printed it.

Contributions are accepted under the same licence. A pull request is taken as agreement
that the contribution may be distributed under BSD-3-Clause.

### Citing it

[`CITATION.cff`](https://github.com/jamccoy/gratinglab/blob/master/CITATION.cff) is
the citation of record, and GitHub's "Cite this repository" renders it as BibTeX or
APA. It names the version, which is checked against `gratinglab.__version__` on every
release, so a citation generated from it points at the code that produced the number
rather than at whatever is current. If a result depends on a solver backend or on an
optical constant table, cite those too: the method's own reference is in
[`docs/references.md`](https://github.com/jamccoy/gratinglab/blob/master/docs/references.md),
and `Provenance` records what a run actually used.

### The licences you inherit by installing it

| Installed with | Package | Licence |
|---|---|---|
| core | numpy, scipy | BSD-3-Clause |
| core | pydantic | MIT |
| core | [endstation](https://pypi.org/project/endstation/) | BSD-3-Clause |
| `[gui]`, `[metrology]` | matplotlib | PSF-based, BSD-compatible |
| `[gui]` | PySide6-Essentials | LGPL-3.0 |
| `[refractiveindex]` | PyYAML | MIT |

A solver-only install is permissive throughout. The `gui` extra is the one worth
knowing about: Qt for Python is LGPL-3.0, which is the reason this project depends on
PySide6 and not on PyQt6, whose licence is GPL and would reach your own code.
Importing an LGPL library dynamically does not change the licence of what you write,
but redistributing a bundled application that ships Qt inside it brings LGPL
obligations for that component — chiefly your users' right to substitute their own
build of it. Running the window yourself, or installing either extra from PyPI, raises
no such question.

### What this licence does not cover

**Optical constant tables.** None are distributed with the package, and a table you
add is governed by the terms of whoever published it, not by this licence. CXRO and
refractiveindex.info each state their own, and a table exported from another program
carries whatever that program's terms attach to it.
[`SOURCES.md`](https://github.com/jamccoy/gratinglab/blob/master/src/gratinglab/materials/data/SOURCES.md)
records where each kind comes from and why none ships here.

**Measurement data.** AFM scans and the reference efficiency corpus are the research
group's data and are deliberately absent from the repository, as described above under
Install. The code is open; the measurements it was validated against are not ours to
publish, and the documentation says which results are therefore not reproducible from
a bare clone.
