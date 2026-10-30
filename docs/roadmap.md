# Roadmap

This page describes the current state of the project, the planned work and the
reasons for its ordering. It is the recommended starting point.

---

## The goal

The project aims to be an open platform for grating groove metrology and
efficiency analysis, for gratings in general: reflection and transmission,
metallic and dielectric, bare and coated, from X-rays to the infrared, in any
mount. It began with one application, soft X-ray spectroscopy in the extreme
off-plane mount, and that application shaped the first solvers. It is no
longer the scope (decided 2026-10-24); the regimes it now covers, and which
method serves each, are mapped below.

The product is not any one solver. The following do not exist openly:

1. A neutral cross-method harness with one problem spec and rigorous conventions.
2. An open C-method (Chandezon) solver.
3. An open integral-method solver.
4. A path from a measured surface to a rigorous efficiency, open end to end.

The fourth arrived with M17 and has no open equivalent at all. It was also
useful immediately, because it falsified an earlier assumption; see the TASTE
entry under Open questions.

The organising principle is that the physical problem and the numerical method
are separate, and that the problem specification is serializable. One `Problem`
reaches every solver, so disagreement between methods becomes visible instead of
being hidden behind incompatible conventions.

RCWA is a backend here and not the product. Open RCWA is already a crowded
field (S4, grcwa, nannos, torcwa, inkstone, RETICOLO), and nothing in this
project competes with those as a standalone code. An RCWA backend is
nevertheless planned, written in-house (§1): it is the natural method for the
dielectric, lamellar, coated and transmission cases a general platform has to
cover, and it is the independent rigorous cross-check on the integral solver
wherever both apply.

### Regimes and the methods that serve them

| Regime | Typical use | Methods | State |
|---|---|---|---|
| Soft X-ray and EUV, grazing incidence, off-plane | X-ray spectroscopy, the founding application | integral (rigorous), scalar (fast) | integral checked order by order against an independent code; corner-graded mesh |
| Soft X-ray and EUV, grazing incidence, in-plane | grazing monochromators, beamline optics | integral, scalar | served by the same solvers, less exercised |
| UV to near-IR, in-plane or conical, metallic | spectrometers, monochromators, pulse compression | integral (finite conductivity), C-method, RCWA | integral reproduces the published benchmark tables; the three rigorous solvers agree with one another to 2e-4 RMS on the TASTE grooves |
| Dielectric, lamellar and transmission gratings | beam splitters, telecom, transmission spectroscopy | RCWA; integral for one interface into a lossless substrate | RCWA agrees with the integral solver on lamellar gratings; transmitted orders are recorded in `EfficiencyScan.transmission` |
| Coated and multilayer gratings | EUV multilayer gratings, protected and overcoated metals | RCWA | `Problem.layers` (conformal, rounded corners) solved by the RCWA; checked against the integral solver where it can say anything, slow in metallic TM |
| Crossed (2D), volume, radial and VLS gratings | photonic structures, imaging spectrographs | later | not started |

The validation logic of the first milestones generalises rather than changes.
Nearly all published benchmark data (Moharam & Gaylord, Li, Popov/Nevière) is
in-plane in the UV to near-IR, so a method is validated there first and then
trusted in regimes the literature does not reach, such as an off-plane mount at
2 nm. Measured efficiencies (the als_eff_meas split under "Later") are what
take a regime from validated against the literature to validated against the
bench.

---

## Current state

| Component | Status |
|---|---|
| Conventions ([`conventions.md`](conventions.md)) | normative, settled |
| Geometry: grating equation, orders, blaze relations | done |
| Illumination: direction cosines, three mount constructors | done |
| Profiles: 4 representations for 4 methods, undercut support | done. `Trapezoidal` added for the shape an etched grating actually has: a sidewall angle and a corner radius, reducing to `Lamellar` exactly when vertical and sharp. `Lamellar` is untouched, being the idealisation the benchmark tables are computed for |
| `Problem`, `Solver` protocol, registry | done |
| **Scalar solver** ([theory](theory/scalar.md)) | done, validated. M19 added the reduced-wavelength validity guard, cast-shadow visibility (`visibility="horizon"`, opt-in, sub-cell boundaries) and the cost of the s/p degeneracy |
| Efficiency-table importer (external codes) | done |
| `.ggp` boundary profiles | done, now as a single implementation |
| **Groove metrology** ([theory](theory/metrology.md)): AFM scans, blaze angles, boundary profiles | done, absorbed at M17 |
| Measured groove → solver, in-process | done, `BoundaryProfile.to_problem` |
| Comparison harness | done. `tests/test_cross_method.py` commits the scalar-vs-integral agreement and divergence ladders and needs no corpus |
| Physics self-checks: reciprocity, energy balance | done |
| Convergence harness ([`convergence.py`](../src/gratinglab/convergence.py)) | done |
| Progress + cancellation on the `Solver` protocol | done |
| Mutation testing ([`mutation-testing.md`](mutation-testing.md)) | done, run on demand; `checks.py` and `geometry.py` at 100% |
| Materials: readers, user libraries, Fresnel, roughness, absolute efficiency | done. Three formats are read: CXRO/Henke exports (X-ray), refractiveindex.info YAML (ultraviolet to infrared, all nine dispersion formulas, read over the whole 3,720-file database) and PCGrate `.ari`. A coating comes from the user's own library, `~/.gratinglab/materials` or `GRATINGLAB_MATERIALS`. No tables ship with the package until CXRO's redistribution terms are settled ([`SOURCES.md`](../src/gratinglab/materials/data/SOURCES.md)) |
| Resolving power ([`resolution.py`](../src/gratinglab/resolution.py)) | done: `R = |m|·N` together with the finite-N line profile from ISSI eq. (8), each of which tests the other; measured `N` arrives from metrology via `to_problem`. Groove-error-degraded R is the named successor |
| **Tip deconvolution** ([`core/tip.py`](../src/gratinglab/metrology/core/tip.py), [theory §4b](theory/metrology.md)) | done at M22: Villarrubia erosion with the certainty map, opt-in via `tip_correction='erosion'`, validated pointwise on a committed dilated fixture. Applied to the real TASTE scan, it narrowed the rounded-groove open question by excluding a nominal tip ([`findings.md`](findings.md)) |
| **SEM micrographs** ([`metrology/io/sem.py`](../src/gratinglab/metrology/io/sem.py)) | done: Zeiss SmartSEM TIFFs read in-house, with the instrument's own pixel size, stage tilt and foreshortening. A second metrology input beside the AFM, and the only one that can see an undercut or measure the width of a crest. Validated on 379 of the group's micrographs and by a synthetic fixture the suite builds for itself |
| **Cross-section to profile** ([`boundary/cross_section.py`](../src/gratinglab/metrology/boundary/cross_section.py), [theory §6](theory/metrology.md)) | first slice done: `trace_cross_section` traces the groove off a cleave edge (first crossing of a half level, sub-pixel down the columns and across the walls), fits the ridge lattice, averages the grooves into a `FromProfileData` and fits a `Trapezoidal` by orthogonal distance with groove-to-groove uncertainties. Foreshortening is an explicit input. Validated on rendered cross-sections of a known groove: depth to 0.1 nm, walls to 0.02 nm, angle to 0.05°, errors matching the scatter between images. Not corrected: the beam's rounding of the radius and a bright rim's outward shift of each wall, both measured. Not yet: undercut walls, a capped cross-section, and a run on the real OGRE micrographs, which are not on MC1 or MAC1 |
| GUI (Qt/PySide6): tabs, geometry dock, 3D conical view | done |
| CI (Linux py3.11/3.12, macOS and Windows py3.12) | green. CI runs with `PYTHONWARNDEFAULTENCODING=1`, so a text `open()` that leaves its encoding to the platform fails on every platform, Windows included |
| **RCWA** ([theory](theory/rcwa.md)) | **first slice done** (0.2.0, 2026-10-27): Li's rules, S-matrix, vectorial and conical, lamellar exact and other profiles staircased, transmitted orders, conformal layers (`Problem.layers`). Validated on lamellar gratings against Fresnel, Airy, reciprocity, the integral solver (to 8e-4 in metallic TM) and three published values (Granet's metallic TM, Lalanne & Hugonin's dielectric, Gralak's silicon rods); two published TE/s values are not reproduced and are recorded. A GUI tab (0.2.1): truncation, slices, and the form refusing what the solver would; conformal layers on the geometry panel. Not yet: staircase accuracy in metallic TM, layers of a second material against an independent code. Not experimental since 0.3 |
| **C-method** ([theory](theory/cmethod.md)) | first slice done (0.2.0): perfect conductor in any mount and one interface into a material, absorbing or lossless, conical and coupled. Validated against Goray & Schmidt Table 3 (Li's C-method), Fresnel, Maréchal–Stroke and the integral solver. Refuses vertical facets and undercuts; slow on sawtooth corners and unstable at high truncation on profiles deeper than about P/2, both reported in provenance. Refuses walls too steep for the truncation since 0.3 (#7). Not experimental since 0.3. No coatings or layer stacks. GUI tab since 0.2.1 |
| **Integral method** ([theory](theory/integral.md)) | perfect-conductivity and finite-conductivity milestones done. Conical from the start, undercut-capable, with a GUI tab; validated against M&P Table 4.1 and, order by order, against an independent integral-method code on the private corpus (a comparison that also identified the TASTE profile). Not experimental since 0.3, after the TASTE measured comparison. M20 added `conductivity="tabulated"`, the coupled conical system of Goray & Schmidt (2010): absolute efficiencies with absorption recorded, a flat limit equal to Fresnel to 6e-5, their Table 3 reproduced per order, and both boundary conditions selectable in the GUI. A graded corner mesh (`corner_grading`) is done for both boundary conditions |

The suite holds roughly 2,500 tests (2,539 passed and 52 skipped on
2026-10-30); the skipped tests are those that need the
private reference corpus or the AFM scans. The figure is deliberately rounded:
this line once claimed 328 and was wrong by 400, so `pytest -q` is the authority
and the figure here gives only the order of magnitude.

---

## What is next, in order

**Order from here, as of 2026-10-30.** The numbered sections below keep their
numbers, because other documents cite them by number; this list is the order of
work. The first three items of the 2026-10-24 order are done:

1. ~~Measured efficiency in the comparison harness~~: the gratinglab half of
   the als_eff_meas split shipped in 0.1.12 (2026-10-24); see "Later". The
   als_eff_meas half, its `to_measured()` export, remains.
2. ~~RCWA, written in-house (§1)~~: first slice in 0.2.0 (2026-10-27), with
   layer stacks from the start; §1 says what remains.
3. ~~The C-method~~: first slice in 0.2.0 (2026-10-27), refusing walls too
   steep for its truncation since 0.3.0.

What is left, as it becomes useful: groove-error-degraded resolving power
(§2), the analytic profiles in the native format (§5), the PyXFocus formulas
(§3) and the JOSS submission (§6). The solvers' experimental status, which the
first three served, was lifted in 0.3.0 on the TASTE campaign; see
"Versioning".

**Run (2026-10-28): the comparison with the ported ALS campaign.** This was
deferred on 2026-10-24 for want of the beamline tree, which has since been
found on a third machine. The comparison runs, and it does not clear the gate.
The rigorous solver reaches an RMS of 0.0706 in absolute efficiency over the
251 measured points, and the discrepancy survives every recorded uncertainty,
both profile families available and the removal of the absolute level. The
obstacle is that the campaign does not record the grating's groove, so the only
comparison available is a two-parameter fit, which cannot separate a profile
error from a solver error. What would close it is metrology on this grating,
not more solving. The solvers therefore stayed experimental on this campaign;
see [`findings.md`](findings.md), "The measured comparison ran". The TASTE
campaign, which does record its groove, cleared the gate the next day, and
0.3.0 shipped on it; see the versioning note under 0.3.

### 0. Integral method, moved ahead of RCWA

An earlier revision of this file said "do not attempt the integral method before
the first release" and placed it as phase 5. That ordering was reversed for three
reasons that the earlier revision could not have known:

1. **The corpus can only validate an integral solver.** The one rigorous
   reference dataset on disk is perfect-conductivity integral-method output from
   a commercial code. The scalar solver could test no more than the diffraction
   pattern against it, and a like-for-like comparison required a solver of the
   same kind.
2. **The off-plane X-ray application was then the purpose of the project, and
   the integral method is the method of record there.** The conical decoupling
   theorem makes the perfectly conducting milestone cheaper at grazing
   incidence, because the mesh has to resolve only the reduced wavelength
   λ/sin γ.
3. **RCWA was then expected as contributed work, and was blocked on the
   contributor.** That expectation was dropped on 2026-10-24: RCWA will be
   written in-house (§1). The second-implementation stress test of the solver
   protocol (Capabilities, progress, tabs, convergence knob) was carried out with
   the integral backend instead, and the protocol held; see the current-state
   table.

The decision was justified immediately. The first rigorous-vs-rigorous comparison
matched the reference order by order, within the reference's own energy-balance
error, and it also identified the previously unconfirmed groove profile of the
TASTE run ([`findings.md`](findings.md)).

The remaining integral milestones, in order, are:

- ~~**Finite conductivity**~~: done at M20. It implements the coupled conical
  system of Goray & Schmidt (2010) as `conductivity="tabulated"`, giving absolute
  efficiencies with the absorbed fraction recorded and `R + A = 1` as a checkable
  theorem ([`theory/integral.md`](theory/integral.md) §8). The Leontovich rung
  was deliberately skipped: it requires |n| ≫ 1 and is therefore invalid in the
  soft-X-ray regime the milestone exists for, and it builds no operators that the
  full system reuses. The mode was validated without a finite-conductivity
  reference run (see the narrowed open question below): a flat interface
  reproduces Fresnel exactly, the energy theorems hold two-sided, the
  perfect-conductor limit agrees with the integral-method agreement already
  established, and Goray & Schmidt's own Table 3 is reproduced within 6e-5 per
  order. The last of these also showed their Table 4 to be a publisher's
  duplication ([`findings.md`](findings.md)).
- ~~**Graded corner mesh**~~ (M&P §4.6.6): done for perfect conductivity, as
  `corner_grading=p` (opt-in, 4 recommended). Cornered TM went from first order
  to an energy-balance deviation of ~1e-5 at P = 256, and reciprocity improved
  by the same factor; see [`theory/integral.md`](theory/integral.md) §5. **With
  finite conductivity it is done too** (2026-10-24): `R + A` reaches 2e-6 to
  3e-5 at P = 512 on a visible-light metal and on conical X-ray gratings. It took
  a corner singularity subtraction of `K⁺` and the tangential derivative taken on
  the assembled trace, not Kress's near-corner modification, and the negative
  permittivity turned out not to matter ([`findings.md`](findings.md)).
- ~~A solver-contract conformance test~~: done; see the note under RCWA below.

### 1. RCWA: an in-house backend

**First slice done, 2026-10-27**: `solvers/rcwa/`, with the design note and
its validation in [`theory/rcwa.md`](theory/rcwa.md). Every requirement below
is met in it. What remains is listed in that page's §9. The three decisions
it raised for the package were settled the same day: a coated grating is
`Problem.layers` of conformal layers with rounded convex corners, a staircase
is allowed by default with its TM warning, and transmitted orders are recorded
in `EfficiencyScan.transmission`.

RCWA will be written within this project. Until 2026-10-24 it was expected as
contributed code (`fgrise` has an implementation), and this section was a list
of requirements to state to a contributor; it is now the plan of work. It comes
after the measured-efficiency validation (see the order above), so that the
integral solver, the method of record for the founding application, is trusted
against the bench before a second rigorous method is built beside it.

Generalising the project is what makes this backend necessary rather than a
cross-check: dielectric, lamellar, transmission and coated gratings are its
natural territory, and in some of them the integral solver as it stands cannot
help at all.

The physics requirements are unchanged and non-negotiable:

- Li's Fourier factorization rules must be used from the start; otherwise TM
  metallic cases converge badly and the comparison plots will misrepresent the
  method.
- S-matrix / enhanced transmittance-matrix propagation is required, never plain
  T-matrix.
- The formulation must be fully vectorial and conical, since off-plane incidence
  is a primary use case from the outset.
- Layer stacks above and below the grating region from the start: a coated or
  multilayer grating is a row of the regime table that nothing else here can
  serve, and a layer stack is a large retrofit into an RCWA code written
  without one.

The framework asks for the following:

| Requirement | Why, and what reads it |
|---|---|
| A `Capabilities` declaration | `conical`, `polarizations`, `handles_undercut`. `Capabilities.check` refuses an out-of-scope case before solving instead of approximating it |
| `accuracy_knob` naming the truncation parameter | `convergence.check_convergence` sweeps it. Without one, a result can never report `converged=True` |
| `accuracy_ceiling`, if the knob has a largest accepted value | `convergence.default_ladder` stops there, so a sweep that has not settled reports that rather than failing on a rung the solver refuses. The C-method declares 256 |
| `reports_progress`, if `solve` takes a `progress` keyword | Opt-in, so that a backend predating the hook is never handed a keyword its signature cannot take. It also makes the solver cancellable, since that callback is the only check point a NumPy-bound solve has |
| `UnsupportedConfiguration` rather than a quiet approximation | This is the rule on which the whole comparison rests: a method that silently smooths a vertical facet returns a plausible but wrong number |
| One entry in `_TAB_FACTORIES` | Tabs are generated from `available_solvers()`; the GUI needs nothing else |

The interfaces have now met their second implementation, the integral backend.
Most of what the earlier warning here predicted did not happen: the
`Capabilities` declaration, the progress and cancellation contract, the
convergence knob and the per-solver tab pattern all accommodated a slow,
rigorous, polarization-resolving backend without changing shape. Two things did
change. The window's "solve on open" behaviour had to prefer the fast solver once
a solver taking minutes existed, and the tab layer now consists of two
near-identical files, which is the evidence the base-class extraction was
waiting for.

A conformance test that every registered solver must pass has now been
[written](../tests/test_solver_conformance.py) and is parametrised over
`available_solvers()`. A new backend, the in-house RCWA included, is therefore
held to the declaration, refusal, scan and progress contracts from the moment it
registers, before any test has been written for it specifically.

### 2. Resolving power: the ideal grating done, the measured grooves to follow

`gratinglab.resolution` now answers the ideal-grating question. It provides
`R = |m|·N` by the Rayleigh criterion and the finite-N line profile from which
that result is derived (`interference_factor`, moved from the scalar solver to
`geometry`, where non-solver consumers can reach it). The measured groove count
arrives through `BoundaryProfile.to_problem`, so an AFM scan leads to a
spectrograph resolving power in-process. A problem without `n_grooves` is
refused and is not treated as infinite.

The successor is groove-error-degraded resolving power: the line spread of a real
groove ensemble arising from period jitter and placement errors (quantities that
metrology already partly holds, `y_std_nm` and `groove_centers`), with the
Rayleigh criterion applied to the degraded profile. This would connect the
metrology of groove errors to spectrograph performance, and it is deliberately
not approximated by the ideal formula now in the API.

### 3. ~~PyXFocus analytics: porting the formulas~~: done, unreleased (2026-10-30)

`geometry.littrow_alpha` and `Illumination.littrow` are in, and
`tests/test_littrow.py` holds them to the transcribed PyXFocus functions on six
mounts from 1.5° off-plane to in-plane. Two things the plan below did not know:

- **The frame mapping it proposed was wrong.** PyXFocus's `inc` is the graze
  onto the grating surface and its yaw is measured in the surface from the
  grooves, so `cos inc` is not `sin γ`; at γ = 1.5° and α = 30° they differ by
  0.97. The only quantity the frames share is the dispersion cosine,
  `cos(inc)·sin ψ = sin γ·sin α`. Conventions §10, row 8.
- **There was nothing separate to port for the blaze angle.** In our frame the
  facet angle that blazes at Littrow is the Littrow azimuth itself, δ = α,
  because `blaze_direction` gives 2δ − α. `blazeAngle2` equals it to 1e-12 and
  `blazeAngle` is the same angle with the opposite sign. `blazeYaw` is singular
  in-plane, where its yaw is exactly 90°.

The plan as it stood:

PyXFocus (the sequential X-ray raytracer that this project's author forked)
stays where it is. It is unpackaged, pinned to numpy<2 and built with Fortran, so
absorbing it would add a liability without adding a capability. The part worth
porting is small and analytic, and comes from its `grating.py`:

- `blazeYaw`, the Littrow condition, as `geometry.littrow_alpha` together with an
  `Illumination.littrow(order, wavelength, period, graze)` constructor.
- `blazeAngle`/`blazeAngle2`, the facet angle that Littrow requires at (λ, m).
  The first step is to verify numerically that the two agree and to pin the frame
  mapping (its `cos inc` ↔ our `sin γ`) against `Illumination.offplane`,
  following the discipline of §10 of the conventions: transcribe, test, then port
  one form.
- `litBetaAlpha`, `eta`, `yaw`: already representable with the `Illumination`
  machinery, so these are covered with tests instead of being ported.
- `radialApproxEffect`: not ported. It implements radial/variable-line-space
  physics that `Problem` cannot represent (and the source contains a live
  `pdb.set_trace()`). Radial/VLS gratings are a "Later" item.

Its raytraced resolving-power estimator (order-centroid spacing over spot width)
is the instrument-level cross-check of the grating-level figure from
`resolution.py`. That comparison belongs with any future instrument layer and
falls outside this item.

### 4. ~~AFM tip deconvolution~~: done at M22

Villarrubia (1997) morphological erosion with the certainty map was implemented
as planned. It uses a parametric cone-plus-spherical-apex tip and runs as a stage
in `metrology/core/tip.py` on the 2-D array after image flattening, in both
branches. It is opt-in through `AnalysisSettings` (`tip_correction`,
`tip_radius_nm`, `tip_half_angle_deg`), and the certain fraction is recorded in
the metrics, the sidecar and the log. It was validated on a committed dilated
fixture, namely the sharp synthetic scan imaged through a known worn tip with the
identical noise field, so that recovery can be asserted pointwise. Recovery
reaches machine precision on clean data and the noise floor (0.18 nm mean) on the
fixture, against a 9 nm mean error without correction.
[`theory/metrology.md`](theory/metrology.md) §4b carries the three load-bearing
statements: erosion gives a bound and cannot resurrect the surface; the
certainty map is the deliverable; and a facet steeper than the flank is lost
entirely.

The run on the real TASTE scan produced the milestone's finding
([`findings.md`](findings.md), "A nominal tip does not explain the rounded
groove"). Erosion with the nominal probe (R ≈ 1–2 nm, 18°) recovers no depth,
because the measured surface is already reachable by the tip, and reproducing the
depth deficit in the forward direction requires an apex worn to ~80 nm. The open
question therefore narrows from "grating or tip?" to "grating, or a severely worn
tip?". That question can be answered outside the software, by the probe's
condition at scan time, a re-scan with a fresh tip, or a cross-section.

Extended 2026-10-22: the tip model now covers asymmetric pyramid tips. Each
flank has its own half angle (`tip_half_angle_right_deg`, `None` for the old
symmetric cone), with Bruker ScanAsyst-Air presets across and along the
cantilever in `core.tip.TIP_PRESETS` and a Tip shape picker in the metrology
window. The erosion now uses the reflected tip, which a symmetric tip never
required. Run forward again on the ideal TASTE groove, the model narrows the gap
without closing it: with the 65° back face facing the anti-blaze facet, it
implies 26.4° at the nominal 2 nm and 25.1° at the worst case the specification
allows, against the measured 20.3° ([`findings.md`](findings.md)).

Two things deliberately did not land. Measured tip images are not accepted as
input (blind tip estimation is the natural successor if a worn tip becomes the
suspect), and the certainty map is not coupled into `to_problem`: a profile with
an uncertain trough still converts, and the defensibility warning remains the
depth-vs-facet-angle diagnostic named in the finding.

### 5. Native boundary format

The format is complete for measured profiles. `gratinglab.io.boundary` writes and
reads a versioned `.json` that carries the period in nm, the groove count, the
metrics and the provenance, and its `to_problem` equals
`BoundaryProfile.to_problem` field for field
(`tests/metrology/test_boundary_pipeline.py`). Its `t` may double back, so an
undercut groove is representable, although nothing measured produces one yet. The
metrology window exports the format, and the efficiency window loads it together
with the period. Two things remain: the analytic profiles (`Blazed`, `Lamellar`,
`Sinusoidal`) are not written to it, and `benchmarks/corpus.toml` still holds
periods for `.ggp` files written before the format existed.

The original motivation for the format still applies, and is as follows.

`.ggp` cannot carry the period in nm, provenance, undercut boundaries or a format
version. The missing period is exactly why
[`benchmarks/corpus.toml`](../benchmarks/corpus.toml) has to exist as a sidecar.
JSON was chosen because the profile classes are already pydantic models with
tested round-tripping.

M17 reduced the urgency of the format without removing the need for it. With
metrology in the package, a measured groove reaches a solver through
`BoundaryProfile.to_problem` without touching disk, carrying its measured period
and optionally a fitted facet angle. That closes the gap completely for profiles
that never leave the process. Until the format existed, everything that had to be
saved, sent, archived or diffed passed through `.ggp`, which drops the period on
export, cannot record where the profile came from, and cannot represent an undercut
boundary in either direction.

### 6. First release + JOSS

0.1.0 was the first release and 0.3.2 (2026-10-30) the latest. Every release
is published to PyPI by Trusted Publishing from
`.github/workflows/release.yml` on a `v*` tag. That workflow refuses a tag that
does not match `pyproject.toml`, checks that the wheel and sdist carry no optical
constant table, and installs the wheel into a clean environment to check that a
missing coating is refused with a helpful message and that a user library
resolves.

Two items remain: a Zenodo DOI (Zenodo archives only releases made after it is
switched on, so the archive starts from the next one) and the JOSS submission.
The release workflow's header says Zenodo listens for the GitHub Release, but a
search of Zenodo's public records on 2026-10-30 found no gratinglab entry, even
after 0.3.2. So either the integration is not switched on for this repository,
or its records are not public yet; the Zenodo GitHub settings page is where to
check.
Shipping something citable and correct early is what recruits collaborators, and
the reordering above strengthens the case: an open integral-method solver checked
order by order against an established code gives a reviewer something to trust,
and the scalar/integral comparison plots need no PCGrate licence to reproduce.

Bundled optical constants are an open question and are not counted as a gap.
refractiveindex.info is CC0 and could ship; CXRO's terms are unconfirmed, so for
now nothing ships and every coating is the user's own table.

### Later

Measured-profile fitting, crossed/2D gratings and radial/VLS gratings are
planned for later. (The C-method was on this list and shipped in 0.2.0.) See [`references.md`](references.md) for the literature
mapped to each.

#### Measured efficiency and the als_eff_meas split

**The gratinglab half is done (0.1.12, 2026-10-24).** `MeasuredEfficiency`
lives in `measured.py`, and `compare.align`/`records` accept it beside
`EfficiencyScan`, with `Comparison.measured_order` for error-bar plots and
`Comparison.residuals` for the model at each measured point. The model is swept
on the measured wavelengths and never interpolated onto them. The PCGrate XML
dialect is read by `read_xml_table`. The XML reader keeps the original's
positional layout and takes wavelengths from the caller, because no real export
was in reach to confirm tag names or a scan attribute; those are the two things
to check when one is.

Still to do on the als_eff_meas side: its `to_measured()` export, whose
conversion is written out in the worked-comparison fixture in
`tests/test_measured.py`, and the deletion of its duplicated optics. The
conversion is now confirmed against real data: the beamline tree was found on a
third machine (not MC1 or MAC1) and the fixture passes, 251 points across orders
−2 to 7, so `to_measured()` has a specification that has been run rather than
only written.

Two of the three things that comparison was to settle are settled. The
recovered `alpha` is already in `Illumination.alpha_deg`'s convention, with no
flip, pinned by
`test_als_eff_meas_recovers_alpha_in_this_projects_sign_convention`; and the
unrecorded polarization turns out not to matter at this graze, TE and TM
differing by at most 9.1e-4 over all 251 points. The placeholder groove depth
remains, and is now the thing blocking the 0.3 gate rather than a detail of the
test: the campaign does not record the groove, and no rectangle or trapezoid
reproduces the measurement. See [`findings.md`](findings.md), "The measured
comparison ran".

The design it was built to, kept for the reasoning:

The sibling `als_eff_meas` package reduces ALS beamline scans to absolute
measured efficiency per order, carrying uncertainties end to end. It is the
measured-data half that this platform still lacks. It is also specific to
X-ray/EUV testing at one facility, so the integration is planned as a split,
in which each package keeps the parts suited to it:

- **Moves here, when taken up:** a measured-efficiency data model
  (`MeasuredEfficiency`, in a new `measured.py`) that `compare.align`/`records`
  accept beside `EfficiencyScan`, so that measured, scalar and integral results
  appear on one plot with error bars. It will not be a second `EfficiencyTable`
  (that name is taken by the modelled-table importer in
  `io/efficiency_table.py`), nor uncertainty fields grafted onto
  `EfficiencyScan`, whose dense-grid invariants are shaped by models; measured
  data is sparse per order, with a σ per point. The model uses parallel
  `value`/`sigma` float arrays, so the `uncertainties` package never becomes a
  dependency here. The PCGrate XML dialect from its `models.py`, which
  `io/efficiency_table.py` (text dialects) lacks, also moves here.
- **Stays there:** everything specific to the beamline, namely file pairing, dark
  subtraction, three-point flux sums, the ODR arc fit and cone-angle recovery,
  and TOML campaigns. `als_eff_meas` gains a gratinglab dependency, deletes its
  duplicated Fresnel/Névot–Croce/CXRO code in favour of `gratinglab.materials`,
  and exports a `to_measured()`.
- **The scipy pin stays out:** `als_eff_meas` pins `scipy<1.19` because
  `scipy.odr` backs its arc fit and is scheduled for removal. The arc fit does
  not move, so the pin does not propagate. That package must eventually migrate
  to `scipy.optimize.least_squares` in any case, and this project should not
  inherit the debt in the meantime.

The other siblings are settled: `afm_blaze_meas` is already absorbed (M17),
`endstation` remains the shared scaffolding dependency, `ebl_fracture` remains a
separate fabrication-layout tool (it contains no optics, and its own roadmap
excludes blazed and VLS geometry), and PyXFocus is covered by §3 above, which
ports its formulas and leaves the raytracer where it is. The subsection below
proposes the one interface `ebl_fracture` would need, and records why the
separation should stay a separation.

#### A design description, and where the fabrication boundary falls

**Proposed, not started (2026-10-24).** Nothing is written and nothing is
committed to; this records the shape and the reasoning so that the question is
settled deliberately rather than by whoever next feels the gap.

The package describes a grating that *exists*. `BoundaryProfile` carries a
measured groove, and `Problem` carries the period, the profile and the
materials of a grating in hand. There is no object for the grating one intends
to *make*: period, duty cycle or blaze angle, patterned area, groove
orientation, and the substrate it goes on. That is the gap a design section
would fill, and the point worth noticing is that it is *physics and not
fabrication* — vendor-neutral, with no fields, no subfields, no dose and no
file format in it.

What it buys is the loop. Design → predicted efficiency from the solvers →
fabricate → measure with the metrology → measured efficiency → compare against
the prediction. Every stage of that already exists in this package or its
siblings, including the last one now that `MeasuredEfficiency` and
`Comparison.residuals` are in. What is missing is a single shared description
of the intended grating, so that the prediction and the thing actually written
to the resist are known to be the same grating rather than two sets of numbers
typed twice.

**The interface, if it is built.** `ebl_fracture`'s `txl2gpf` already accepts a
grating as loose arguments (`--grating 1.0 --duty 0.5 --area 1000x1000`). It
would accept the design object instead, with `gratinglab` an *optional*
dependency of `ebl_fracture` and never the reverse. This package would learn
nothing about GPF, fields, dose or any instrument.

**Why not absorb the fracture tool**, recorded so the question is not reopened
by accident:

- **No shared physics.** Nothing in fracture needs the grating equation;
  nothing in efficiency needs a field grid or a dose ladder. The overlap is the
  word "grating".
- **Different validation regimes.** `ebl_fracture` is validated by byte-exact
  round-trips of a vendor binary and by area conservation over ~135,000
  randomised shapes. This package is validated against closed-form results,
  published benchmark tables, energy balance and reciprocity. Both are sound
  and neither transfers.
- **It is pre-instrument.** By its own README nothing it has written has been
  near the tool. Absorbing unexposed, instrument-specific code would make this
  package's releases wait on beam time.
- **The blazed case is its deferral, not ours to take on.** Its roadmap names
  multi-dose and grayscale — that is, blazed profiles — as the *first* widening
  of scope, with the reference material already identified, and gates it on a
  plain uniform grating coming off the machine correctly. That sequencing is
  right for a reverse-engineered writer and would be inherited wholesale by a
  merge. A `txl`-to-`gpf` flow for *blazed* gratings is therefore a later step
  in that project, not an early one here.

**What this is not.** General pattern design beyond a grating — arrays,
alignment marks, dicing lanes, whole-wafer layout — belongs with the layout
tool and not here. The line to hold is that this package describes gratings as
optics; the layout tool describes patterns as exposures.

**The four questions, answered (2026-10-24).** Three of them the existing code
already decides; the fourth is the one that makes the section worth building.

**A new type, with a `to_problem()`.** Not a constructor on `Problem`.
`Problem` is frozen with `extra="forbid"` and is the contract every solver
reads; its fields beyond period and profile — `roughness`, `n_grooves` —
describe a grating in hand under illumination, not an intention. A patterned
area, a groove angle and a tolerance mean nothing to a solver, and putting them
there would oblige every backend to ignore them. The precedent settles the
duplication objection: `BoundaryProfile.to_problem()` is already how a
*measurement* becomes a problem, and it is not accused of duplicating
`Problem`. A design enters the same way, so there are two producers and one
contract rather than two overlapping descriptions.

**The profile union already settles the duty-cycle question.** `Blazed`
(blaze and antiblaze angles), `Lamellar` (`depth_fraction`, `duty_cycle`) and
`Sinusoidal` (`depth_fraction`) exist as a discriminated union, so no struct
ever carries both a duty cycle and a blaze angle and the ignored-field problem
is structurally impossible. Reusing that union, rather than inventing a parallel
vocabulary of intents, is also what makes the comparison below like-for-like:
the design and the measurement describe the shape in the same terms. One
narrowing — `FromProfileData` is excluded, because a measured groove is not an
intention, so a design's profile is a strictly narrower union than a problem's.

**Geometry, plus the two material identifiers `Problem` already carries, and
nothing else.** `coating` and `substrate` are in `Problem` because absolute
efficiency needs them, and a design that cannot name them cannot predict the
efficiency it is being designed for, which is most of the point. Resist, dose,
proximity correction, beam step and field size stop at the boundary: they are
instrument-specific, the fracture tool's settings file already carries them, and
nothing here reads them. Depth intent is not part of this question at all — it
lives in the profile. The borderline case is resist thickness, which bounds the
achievable depth: that is a constraint a *process* places on a design, not a
property of the design, so it belongs to whatever checks one against the other.

**A tolerance, and it is the reason to build this.** Without it a design is a
tidier way to type four numbers; with it, the metrology can answer whether the
grating that came back is the grating that was asked for — and the AFM pipeline
already measures every input that question needs: the period from groove
spacing, the blaze angle from the facet fits, the depth from the averaged
groove, and now the groove angle from the orientation fit. Two conditions on
`check`, both learned the hard way today:

- It must report the measurement's own uncertainty beside the deviation. Per-scan
  SEM on the blaze angle is a median 0.77°, and the 0.05° tolerance pinned into
  `test_both_routes_give_the_same_blaze_angle` is what a threshold set an order
  of magnitude below the uncertainty looks like: it fails on 62 of 69 archive
  pairs and means nothing.
- Where a quantity has no uncertainty attached, it should refuse rather than
  return a verdict, per "Refuse rather than approximate" below. A conformance
  claim without an error bar is the most authoritative-looking wrong answer this
  package could produce.

The per-side facet angle that the avg lineage dropped at rev3 (see the corpus
survey) would give `check` a second blaze estimate independent of groove
detection, trimming and the row-group partition, which is what a conformance
claim actually wants.

**The type, concretely:**

```python
class Tolerance(BaseModel):            # frozen, extra="forbid"
    period_nm: float | None = None     # half-widths, not ranges
    blaze_angle_deg: float | None = None
    depth_nm: float | None = None
    duty_cycle: float | None = None
    groove_angle_deg: float | None = None

class Design(BaseModel):               # frozen, extra="forbid"
    period: float                      # nm, as everywhere (conventions.md §1)
    profile: Blazed | Lamellar | Sinusoidal   # the intent, in Problem's vocabulary
    area_nm: tuple[float, float] | None = None  # patterned extent, fast then slow
    groove_angle: float = 0.0          # deg from the substrate's reference edge
    coating: str | None = None
    substrate: str | None = None
    tolerance: Tolerance | None = None

    def to_problem(self, *, roughness=0.0, n_grooves=None) -> Problem: ...
    def check(self, measured: BoundaryProfile) -> Conformance: ...
```

Every field has a consumer here, which is the test of whether it belongs: the
profile and period feed `to_problem`, `coating`/`substrate` feed absolute
efficiency, `groove_angle` is checkable against the orientation fit, and
`area_nm` divided by the period is the illuminated groove count `Problem`
already takes as `n_grooves`. It is also what `txl2gpf` needs, which is the
point of the boundary above, but it is not *only* what the fracture tool needs.

What is still undecided, and would be settled by writing it:

- The shape of `Conformance`. Per-quantity deviation, uncertainty and verdict
  is the obvious content; whether it is a type of its own or a `Provenance`-like
  record that a summary file can carry is not obvious.
- What `check` takes, which is harder than it looks. The uncertainties it needs
  do not travel together today. `BoundaryProfile` carries `y_std_nm`, the
  per-point scatter of the averaged groove, but the blaze angle's SEM and its
  ICC correction come from the row-group statistics on the other branch of the
  pipeline, and `docs/theory/metrology.md` §5 records that the two branches
  deliberately do not exchange results. A conformance check wants the shape
  from one and the error bars from the other, so either it takes both, or
  something has to carry them together — which is a decision about the
  metrology's own interfaces and not about this section.

#### The groove residual's two follow-ups

The boundary pipeline measures the groove residual,
`BoundaryProfile.groove_residual_rms_nm`
([theory §4](theory/metrology.md)), by registering each groove to the average
by a sub-pixel shift and a height offset before taking the residual. Two
things it turned up are not done:

- **The exported average is blurred by the same misregistration.** It is the
  plain mean of grooves centred on whole pixels, so the set is misaligned by up
  to half a pixel, which rounds the corners of the averaged groove. Averaging
  the registered grooves instead would sharpen them. It changes every exported
  profile and the byte-level `rev3` equivalence pin, so it wants a setting that
  defaults to today's behaviour, and a check on real scans of how far the depth
  and the facet angles move.
- **A per-scan-line residual.** The residual is taken on the row-averaged
  profile, so roughness that varies along the groove averages away over the
  scan lines and the number is a lower bound. Taking it per line before the
  average would keep that roughness, at the cost of carrying the AFM's full
  noise, which a single scan cannot separate from it. Worth doing only with a
  way to estimate that noise, such as trace and retrace or a repeated scan.

---

## Versioning

The project follows the usual convention for software before 1.0. A patch
release (0.1.x) carries fixes and small additions; a minor release (0.x) marks a
milestone and is also where an interface may change, provided the change is
stated in [`CHANGELOG.md`](../CHANGELOG.md); 1.0 is a commitment that the public
interface is stable until 2.0. Releases 0.1.1 to 0.1.13 were all patch
releases, although two of them (swath statistics and asymmetric tips) added
features that would justify a minor version under a stricter reading. 0.2.0
and 0.3.0 are the two minor releases so far.

The intended milestones are as follows.

- **0.2**: the AFM metrology is complete for the group's own analysis. Its
  components are in place: groove averaging that excludes partial grooves, tip
  correction with real probe shapes, swath statistics, and the boundary JSON
  format that carries the period. The release is due once the graduate students
  using it confirm that it covers their analysis.

  **0.2.0 was released on 2026-10-27 before that sign-off**, because the RCWA
  and C-method backends and a change to what `absorption` means each call for
  a minor release under this policy. The metrology milestone is still open,
  and the students' confirmation is still what closes it; it will be recorded
  here and in the CHANGELOG of whichever release follows it.
- **0.3 and later minor versions**: the efficiency solvers lose their
  experimental status. The criterion is the one stated in the README: use on a
  real problem, with the integral solver compared against measured
  efficiencies, such as the ALS data reduced by `als_eff_meas`, and not only
  against published benchmarks and its own physics checks. **Done in 0.3
  (2026-10-29)**: the TASTE comparison met the criterion, and the two solver
  defects it found, integral absorption with tabulated metals
  ([#6](https://github.com/jamccoy/gratinglab/issues/6)) and the C-method on
  near-vertical walls ([#7](https://github.com/jamccoy/gratinglab/issues/7)),
  are dealt with: #7 is refused, and #6 is fixed by the opt-in
  `profile_smoothing="auto"`. The scalar solver is an approximation and was never on this
  path. Two things stay open and are not solver defects: the 4.6% in
  `period · sin γ` on the TASTE campaign, and the integral solver's absorption
  on a tightly rounded valley, recorded as a limitation in
  `docs/theory/integral.md`.

  Attempted 2026-10-28 against the ALS October 2021 OGRE campaign and not met:
  that campaign records no groove, so the comparison was a two-parameter fit
  leaving an RMS of 0.0706 that no recorded uncertainty explained, and a fit
  cannot separate a profile error from a solver error. The gate therefore
  needed a measured groove for the grating that was measured.

  The December 2018 TASTE campaign has one. 296 absolute per-order efficiencies
  from 80 to 800 eV over orders −4 to 8, against the grating's own AFM scan
  through `metrology.boundary`, with every input measured and none fitted:
  RMS 0.0383, or 0.0363 with the roughness factor the paper applies. The
  separation the gate was for is now measured rather than argued — the integral
  method, the C-method and RCWA agree with one another to 2e-4 RMS while all
  three sit 0.038 from the grating, so the residual is not the solver. The
  absolute level is accounted for exactly by gold's reflectivity and the AFM's
  own 1.5 nm roughness. What remains is one scalar, `period · sin γ`, which the
  measurement wants 4.6% smaller than the published numbers allow; that is a
  question for the mount, and it is an open question below. See
  [`findings.md`](findings.md), "The verdict on 0.3, from the TASTE campaign".

  Integrating the in-house RCWA backend was to be a further minor release; it
  arrived first, in 0.2.0.
- **1.0**: the public interface is declared stable. That covers `Problem`,
  `BoundaryProfile`, `AnalysisSettings`, the boundary JSON format and the solver
  plugin protocol. The preconditions are that the solvers have been validated
  against measurements, that these interfaces have stopped changing between
  releases, and that the Zenodo DOI and the JOSS submission (§6) are complete.
  JOSS does not require 1.0, but a reviewed paper is a natural point at which to
  fix the interface it describes.

## Open questions

Answers to the following questions would unblock real work.

### The panter1 grating period

The following values were recovered from the reference data with no assumptions:

- `period × sin γ = 9.8450 nm` (±0.005)
- `α = +32.92°` (±0.1)

These two cannot be separated without one further fact. The candidates are:

| γ | period |
|---|---|
| 1.5000° | 376.1 nm |
| 1.7902° | 315.15 nm (same grating as TASTE) |
| 2.0000° | 282.1 nm |

The title string in the original plotting script reads "η = 1.5° & φ = 0.97°".
The 1.5° is consistent with γ if the period is ~376 nm, but 0.97° is not the
azimuth: the data gives α ≈ 33°, which is also far more plausible given the
~39°/27° facet angles of `panter1.ggp` (near-Littrow). The symbol "φ" therefore
denotes something else in that script, and knowing what it denotes would help in
reading the other legacy scripts.

### Which `.ggp` goes with which run

This question is about the *reference* TASTE grating, the 315 nm one the PCGrate
corpus computes. It is not the 400 nm TASTE grating of the ALS December 2018
campaign, whose groove and whose measured efficiency are both now in hand; see
`findings.md`, "The TASTE campaign". Two gratings from the same programme share
the name, and only one of them is still unidentified.

`AFM_TASTE_test.ggp` and `TASTE_a205_TiPt.ggp` are both candidates for the TASTE
efficiency run, and neither is present in the corpus directory, which holds
`AFM_test.ggp`, `AFM_real_echelle.ggp` and `PCGrateProjects/panter1.ggp`. The
filename in `corpus.toml` therefore does not resolve, which also needs to be
settled. Until it is, the comparison uses an idealised `Blazed(29.5°, 70.5°)`,
which is probably why the scalar total sits near 0.55 instead of closer to unity.

M16 remeasured this: the total is 0.5431 for the idealised sawtooth and 0.6138
for `AFM_real_echelle.ggp`, against the reference's unit total. The flux
obliquity of M16-C moved it in the right direction, but only slightly, because
this mount is close enough to Littrow that the factor is ~0.99 on the orders
carrying the power. The residual is due to profile mismatch and not to
normalisation; see `findings.md`, "The corpus can test the diffraction but not
the reflectivity".

M17 tested that explanation, and the test failed. The metrology merge brought
into the repository an AFM scan of what is plausibly the same grating
(`data/TASTE_ALS_A205_Ti_Pt_flatten.txt`), so the idealised sawtooth could be
replaced by a measured groove. The total moved in the wrong direction, from
0.5420 idealised to 0.3362 measured, against the reference's unit total. The
groove's facet fit (27.91° ± 2.13°) and its depth (which implies 20.33°) disagree
by 7.6°: the groove is rounded and not faceted, and the solver responds to the
depth. A flat land would explain the same deficit and was excluded (it would need
25–30% of the period, whereas 2.9% is within 2° of flat), as was smearing from the
groove averaging (0.4%). The full evidence is in `findings.md`, "The measured
groove is rounded, not faceted".

The open question is therefore no longer only "which `.ggp`". Whether the
rounding belongs to the grating or comes from the tip is undetermined, and in
either case an absolute efficiency from a raw AFM boundary is not yet defensible.
Nothing in the pipeline warns about it. The arbiter has now run (§4 above, M22):
erosion with the nominal tip explains none of the 27.91°/20.33° gap, and
reproducing the gap forward needs an ~80 nm apex. The question is narrowed but
not closed, to "grating, or a severely worn tip?", and what settles the remainder
lies outside the software: the probe's condition at scan time, a re-scan with a
fresh tip, or a cross-section. See `findings.md`, "A nominal tip does not explain
the rounded groove".

The software is now ready for the third of those: `metrology/io/sem.py` reads a
cross-section with its scale attached. The measurement is not. The collection's
cross-sections are of laminar gratings and date from 2021, the TASTE scans from
2022, and the two have never been taken of one sample — so the experiment that
would close this is a cleave and an afternoon on the microscope, not more code.
What the cross-sections there do establish is recorded under "What a
cross-section shows that a scan cannot".

### A finite-conductivity reference run

This question was narrowed at M20 but not closed. The corpus still has no usable
finite-conductivity table (`panter1_finite` stays `usable = false`), but the
rigorous solver no longer depends on one: `conductivity="tabulated"` is validated
on the Fresnel flat limit, the two-sided energy theorems, the perfect-conductor
limit and Goray & Schmidt's published Table 3. A fresh finite-conductivity
reference run would still provide a like-for-like corpus item. The recommended
recipe is TE, Au, the TASTE geometry and the identified `AFM_test.ggp` groove, on
a wavelength grid that avoids ±0.05 nm around every order-passing-off point (the
documented failure locus of `panter1_finite`), with the same `.ari` optical
constants exported by `write_ari` so that both codes use identical n(λ). The run
should be gated with `check_energy_balance` before it is admitted to
`corpus.toml`. It would also finally arbitrate the scalar solver's
groove-resolved reflectivity model, which still rests on internal consistency
plus, more recently, cross-checks against the rigorous absolute mode.

### The integral solver's absorption at large (1 − |n|)

`conductivity="tabulated"` returns an absorbed fraction that is too large when
the coating's complex index departs from unity, and the error reaches 8% of the
incident power. Measured at the TASTE geometry (Au, 400 nm period, conical,
γ = 1.73°), against the C-method and RCWA, which both satisfy `R + A = 1` to
1e-5 there:

| eV | \|n\| | A, integral | A, C-method | R + A, integral |
|---:|---:|---:|---:|---:|
| 800 | 0.9971 | 0.34891 | 0.34883 | 1.00007 |
| 200 | 0.9888 | 0.22232 | 0.22012 | 1.00220 |
| 120 | 0.9445 | 0.06810 | 0.05627 | 1.01182 |
| 80 | 0.8880 | 0.16107 | 0.07875 | 1.08266 |

The reflected orders are not affected: they match the C-method to 2e-5 RMS over
all 296 points of that campaign. The error also *grows* with mesh refinement —
3.9% at 512 boundary points, 8.3% at 1024, 8.1% at 2048 — so it is in the
absorption integral's formulation and not in its discretisation, and it is what
stops that configuration reaching a convergence plateau at 1e-5. Dropping the
80 eV point from the probe moves the last two differences from 2.7e-4, 3.4e-4
to 1.3e-4, 1.9e-5.

What would settle it: a derivation check of the absorbed-power integral against
the two-sided theorem in the limit of a flat boundary, where the answer is
Fresnel and `1 − R` is exact, swept over |n| rather than over mesh. The flat
limit is already a test in the suite; it is run at |n| near unity.
`findings.md`, "The TASTE campaign", has the full table.

### The C-method answers confidently on a wall it cannot represent

The C-method solves a boundary through a coordinate transformation, which needs
a single-valued surface of bounded slope. Given the paper's TASTE trapezoid,
whose anti-blaze side is 88.8° from the surface, it returns RMS between 1.6 and
3.3 against the measurement with individual efficiencies up to 31, where the
integral method returns 0.033 to 0.048 on the same grooves.

It does not refuse, and the reason is that its guard is an equality test.
`cmethod/_series.py` raises when a profile has two heights at one position
(`steps == 0.0`), and refuses `Lamellar` and `Blazed(antiblaze_angle=90)` by
type; a wall at 88.8° has a finite step and passes. So the README's "refuses
vertical facets" is true as written and does not cover the case that matters.

**Answered in 0.3.0 ([#7](https://github.com/jamccoy/gratinglab/issues/7),
closed 2026-10-30).** The guard is no longer an equality test. `wall_floor` in
`solvers/cmethod/__init__.py` reads the truncation a profile needs off its own
Fourier series before solving: 0.75 × max(B, 3s), where B is the number of
harmonics holding all but 1e-3 of the power and s the steepest slope. A floor
above 256 raises `UnsupportedConfiguration`, a truncation below the floor
raises `ValueError` naming it, and between the floor and twice it a provenance
warning says the energy balance can understate the error. The ceiling was
measured rather than chosen, against the integral solver on 24 configurations
with walls from 45° to 89.5°; the commit that added it, `3089e95`, carries the
numbers. The paragraphs above are the record of the question as it was asked.

### The 4.6% in `period · sin γ` for the TASTE campaign

Every feature of the TASTE residual is reproduced by shrinking
`period · sin γ` 4.6% below 400 nm × sin 1.73°: the RMS falls from 0.0383 to
0.014, and the modelled blaze order agrees with the measured one at all 37
energies instead of 31. No candidate survives its own measurement — the period
would have to be 382 nm against an AFM reading 396–412 and a dispersion
re-derivation of 401 ± 8, γ would have to be 1.65° against 1.73 ± 0.03, and the
beamline's energy scale cannot be tested because the scan files' fourth column
is a constant 500.3 rather than the energy readback its name suggests.

What would settle it: a re-derivation of the mount angles from the arc fit
without the small-angle assumptions, or a dispersion measurement on this grating
with a grating of known period beside it. Note the shape of this question is the
same as "The panter1 grating period" above — `period × sin γ` recovered without
being separable — and a fact that separates them for one grating may separate
them for both.

### Smaller ones

- Is the vendor console solver runnable (Windows VM)? The ability to generate
  fresh reference cases on demand, instead of only reusing exported tables,
  would justify an afternoon of effort.
- Is there any synchrotron or Panter measured efficiency beyond the current
  corpus that could seed the benchmark suite? This is now answered twice over.
  `als_eff_meas` holds three fully ported ALS campaigns with absolute per-order
  efficiencies and uncertainties -- ogre laminar (Oct 2021) and the two halves
  of the TASTE run (Dec 2018) -- with ~13 legacy per-experiment scripts left in
  its archive. Both campaigns are reached from the comparison harness in
  `tests/test_measured.py`. OGRE could not validate a solver for want of the
  groove; TASTE has one, and `findings.md`, "The TASTE campaign", records what
  that bought.

  The lesson from porting the second campaign is worth carrying to the third:
  the surviving driver script is not the record of the published run. For TASTE
  it disagrees with the saved output in three independent ways -- a different
  number of excluded points, a radius correction that cannot complete the
  reduction, and a throw that is neither the published nor the assumed one --
  and in each case the saved output is the record.

### The groove of OGRE laminar grating I

This was the 0.3 blocker until the TASTE campaign cleared that gate, and it
is still worth settling: a second measured grating would test the solvers on a
laminar groove, which TASTE's is not. It is metrology rather than
software. The
October 2021 ALS campaign gives 251 absolute per-order efficiencies with
uncertainties and a period cross-checked two ways, but no groove depth and no
duty cycle; the original analysis modelled no groove at all, plotting the
measurement against bare Au reflectivity. Fitting the two numbers leaves an RMS
of 0.0706, with even orders over-modelled and odd orders under-modelled, and
sidewall slope trades the high orders against the low ones without improving
either overall ([`findings.md`](findings.md), "The measured comparison ran").

What settles it: an AFM scan of that grating through `FromProfileData`, with
the erosion tip correction the metrology side already has, which turns the fit
into a prediction and makes the residual attributable. The grating is a PSU
part, so the scan is in reach in a way that the TASTE probe's condition at scan
time is not.

---

## Practices worth keeping

- **Refuse rather than approximate.** A solver asked for something outside its
  declared capabilities raises an exception. Silent approximation is what makes
  cross-method plots misleading.
- **Report, never rescale.** Where an approximate theory violates a conservation
  law, the violation is information. Normalising it away destroys the validity
  map.
- **Provenance is mandatory.** A result whose convergence was never demonstrated
  reports itself as not defensible, and `gratinglab.convergence` is the means by
  which it becomes defensible. `converged=False` is itself a result and does not
  indicate a missing one.
- **A control that cannot do what it says must not be offered.** Cancel now
  stops the solver; previously it could only stop the waiting, and the panel
  stated as much without implying more. A backend that declares no
  `reports_progress` cannot be cancelled, and that fact is shown to the user.
- **One agreement is not a plateau.** Refinement error is not always monotone.
  The harness requires three consecutive knob values to agree, because the
  blazed case measurably violates the tempting rule that a single agreement
  suffices
  ([findings](findings.md#quadrature-error-is-not-monotone-in-the-number-of-points)).
- **Commit before mutating.** Mutation scripts revert with `git checkout --`,
  which silently fails on untracked files.
