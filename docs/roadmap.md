# Roadmap

This page describes the current state of the project, the planned work and the
reasons for its ordering. It is the recommended starting point.

---

## The goal

The project aims to be an open platform for grating groove metrology and
efficiency analysis. It is not intended as another RCWA implementation, since
open RCWA is already a crowded field (S4, grcwa, nannos, torcwa, inkstone,
RETICOLO). The following do not exist openly:

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

As a deliberate consequence, RCWA is treated as a *reference backend* and is not
the product.

### Two regimes, two roles

- **Soft X-ray, extreme off-plane (γ ≈ 1.5°): the application.** This is the
  regime in which the commercial integral-method codes dominate and
  general-purpose open RCWA is weakest.
- **UV–Vis–NIR, in-plane: the validation regime.** Nearly all published
  benchmark data (Moharam & Gaylord, Li, Popov/Nevière) lies in this regime. An
  off-plane X-ray solver is difficult to validate against the literature
  directly, but the same solver can be validated in the visible and then applied
  at 2 nm.

---

## Current state

| Component | Status |
|---|---|
| Conventions ([`conventions.md`](conventions.md)) | normative, settled |
| Geometry: grating equation, orders, blaze relations | done |
| Illumination: direction cosines, three mount constructors | done |
| Profiles: 4 representations for 4 methods, undercut support | done |
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
| GUI (Qt/PySide6): tabs, geometry dock, 3D conical view | done |
| CI (Linux py3.11/3.12, macOS and Windows py3.12) | green. CI runs with `PYTHONWARNDEFAULTENCODING=1`, so a text `open()` that leaves its encoding to the platform fails on every platform, Windows included |
| RCWA | **contributed, not yet integrated** |
| C-method | not started |
| **Integral method** ([theory](theory/integral.md)) | perfect-conductivity and finite-conductivity milestones done. Conical from the start, undercut-capable, with a GUI tab; validated against M&P Table 4.1 and, order by order, against an independent integral-method code on the private corpus (a comparison that also identified the TASTE profile). **Experimental**: not yet used for real analysis. M20 added `conductivity="tabulated"`, the coupled conical system of Goray & Schmidt (2010): absolute efficiencies with absorption recorded, a flat limit equal to Fresnel to 6e-5, their Table 3 reproduced per order, and both boundary conditions selectable in the GUI. A graded corner mesh is the named successor |

The suite holds roughly 1450 tests; the skipped tests are those that need the
private reference corpus or the AFM scans. The figure is deliberately rounded:
this line once claimed 328 and was wrong by 400, so `pytest -q` is the authority
and the figure here gives only the order of magnitude.

---

## What is next, in order

### 0. Integral method, moved ahead of RCWA

An earlier revision of this file said "do not attempt the integral method before
the first release" and placed it as phase 5. That ordering was reversed for three
reasons that the earlier revision could not have known:

1. **The corpus can only validate an integral solver.** The one rigorous
   reference dataset on disk is perfect-conductivity integral-method output from
   a commercial code. The scalar solver could test no more than the diffraction
   pattern against it, and a like-for-like comparison required a solver of the
   same kind.
2. **The off-plane X-ray application is the purpose of the project, and the
   integral method is the method of record there.** The conical decoupling
   theorem makes the perfectly conducting milestone cheaper at grazing
   incidence, because the mesh has to resolve only the reduced wavelength
   λ/sin γ.
3. **RCWA is contributed work and is blocked on a contributor.** The
   second-implementation stress test of the solver protocol (Capabilities,
   progress, tabs, convergence knob) was carried out with this backend instead,
   and the protocol held; see the current-state table.

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
- **Graded corner mesh** (M&P §4.6.6): TM on cornered profiles currently
  converges at first order, and with finite conductivity both polarizations do;
  see [`theory/integral.md`](theory/integral.md) §5.
- ~~A solver-contract conformance test~~: done; see the note under RCWA below.

### 1. RCWA: integrating a contributed backend

The RCWA backend is not being written within this project; `fgrise` has RCWA
code. The physics list below is therefore a set of requirements to be stated to a
contributor and is not a plan of work. The more useful half of this section is
the set of conditions that any backend has to satisfy in order to register at
all.

The physics requirements are unchanged and non-negotiable:

- Li's Fourier factorization rules must be used from the start; otherwise TM
  metallic cases converge badly and the comparison plots will misrepresent the
  method.
- S-matrix / enhanced transmittance-matrix propagation is required, never plain
  T-matrix.
- The formulation must be fully vectorial and conical, since off-plane incidence
  is a primary use case from the outset.

The framework asks for the following:

| Requirement | Why, and what reads it |
|---|---|
| A `Capabilities` declaration | `conical`, `polarizations`, `handles_undercut`. `Capabilities.check` refuses an out-of-scope case before solving instead of approximating it |
| `accuracy_knob` naming the truncation parameter | `convergence.check_convergence` sweeps it. Without one, a result can never report `converged=True` |
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
`available_solvers()`. A contributed backend is therefore held to the
declaration, refusal, scan and progress contracts from the moment it registers,
before any test has been written for it specifically.

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

### 3. PyXFocus analytics: porting the formulas

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

0.1.0 is the first release. It is published to PyPI by Trusted Publishing from
`.github/workflows/release.yml` on a `v*` tag. That workflow refuses a tag that
does not match `pyproject.toml`, checks that the wheel and sdist carry no optical
constant table, and installs the wheel into a clean environment to check that a
missing coating is refused with a helpful message and that a user library
resolves.

Two items remain: a Zenodo DOI (Zenodo archives only releases made after it is
switched on, so the archive starts from the next one) and the JOSS submission.
Shipping something citable and correct early is what recruits collaborators, and
the reordering above strengthens the case: an open integral-method solver checked
order by order against an established code gives a reviewer something to trust,
and the scalar/integral comparison plots need no PCGrate licence to reproduce.

Bundled optical constants are an open question and are not counted as a gap.
refractiveindex.info is CC0 and could ship; CXRO's terms are unconfirmed, so for
now nothing ships and every coating is the user's own table.

### Later

C-method, measured-profile fitting, crossed/2D gratings and radial/VLS gratings
are planned for later. See [`references.md`](references.md) for the literature
mapped to each.

#### Measured efficiency and the als_eff_meas split

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
ports its formulas and leaves the raytracer where it is.

---

## Versioning

The project follows the usual convention for software before 1.0. A patch
release (0.1.x) carries fixes and small additions; a minor release (0.x) marks a
milestone and is also where an interface may change, provided the change is
stated in [`CHANGELOG.md`](../CHANGELOG.md); 1.0 is a commitment that the public
interface is stable until 2.0. The releases so far, 0.1.1 to 0.1.7, have all
been patch releases, although two of them (swath statistics and asymmetric tips)
added features that would justify a minor version under a stricter reading.

The intended milestones are as follows.

- **0.2**: the AFM metrology is complete for the group's own analysis. Its
  components are in place: groove averaging that excludes partial grooves, tip
  correction with real probe shapes, swath statistics, and the boundary JSON
  format that carries the period. The release is due once the graduate students
  using it confirm that it covers their analysis.
- **0.3 and later minor versions**: the efficiency solvers lose their
  experimental status. The criterion is the one stated in the README: use on a
  real problem, with the integral solver compared against measured
  efficiencies, such as the ALS data reduced by `als_eff_meas`, and not only
  against published benchmarks and its own physics checks. Integrating the
  contributed RCWA backend is a further minor release.
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

### Smaller ones

- Is the vendor console solver runnable (Windows VM)? The ability to generate
  fresh reference cases on demand, instead of only reusing exported tables,
  would justify an afternoon of effort.
- Is there any synchrotron or Panter measured efficiency beyond the current
  corpus that could seed the benchmark suite? This is partially answered:
  `als_eff_meas` holds one fully ported ALS campaign (ogre laminar, Oct 2021)
  with absolute per-order efficiencies and uncertainties, plus ~15 legacy
  per-experiment scripts in its archive that could be ported in the same way.
  Reaching that data from the comparison harness is the measured-efficiency split
  under "Later".

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
