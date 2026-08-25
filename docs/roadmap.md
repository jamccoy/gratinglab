# Roadmap

Where this project is, where it is going, and why. Start here.

---

## The goal

An open platform for grating groove metrology and efficiency analysis. Not
"another RCWA" — open RCWA is a crowded field (S4, grcwa, nannos, torcwa,
inkstone, RETICOLO). What does *not* exist openly is:

1. A neutral cross-method harness with one problem spec and rigorous conventions.
2. An open C-method (Chandezon) solver.
3. An open integral-method solver.
4. **A path from a measured surface to a rigorous efficiency, open end to end.**

The fourth arrived with M17, and is the one with no open equivalent at all. It
also immediately earned its keep by falsifying something — see the TASTE entry
under Open questions.

The organising idea: **the physical problem and the numerical method are
separate, and the problem spec is serializable.** One `Problem` reaches every
solver, so disagreement between methods becomes visible instead of hidden behind
incompatible conventions.

Deliberate consequence: RCWA is a *reference backend*, not the product.

### Two regimes, two roles

- **Soft X-ray, extreme off-plane (γ ≈ 1.5°) — the application.** Where the
  commercial integral-method codes dominate and general-purpose open RCWA is
  weakest.
- **UV–Vis–NIR, in-plane — the validation regime.** Nearly all published
  benchmark data (Moharam & Gaylord, Li, Popov/Nevière) lives here. You cannot
  easily validate an off-plane X-ray solver against literature; you *can*
  validate the solver in the visible and then apply it at 2 nm.

---

## Current state

| Component | Status |
|---|---|
| Conventions ([`conventions.md`](conventions.md)) | normative, settled |
| Geometry — grating equation, orders, blaze relations | done |
| Illumination — direction cosines, three mount constructors | done |
| Profiles — 4 representations for 4 methods, undercut support | done |
| `Problem`, `Solver` protocol, registry | done |
| **Scalar solver** ([theory](theory/scalar.md)) | done, validated; M19 added the reduced-wavelength validity guard, cast-shadow visibility (`visibility="horizon"`, opt-in, sub-cell boundaries), and the s/p degeneracy price tag |
| Efficiency-table importer (external codes) | done |
| `.ggp` boundary profiles | done, and now a single implementation |
| **Groove metrology** ([theory](theory/metrology.md)) — AFM scans, blaze angles, boundary profiles | done, absorbed M17 |
| Measured groove → solver, in-process | done, `BoundaryProfile.to_problem` |
| Comparison harness | done; `tests/test_cross_method.py` commits the scalar-vs-integral agreement and divergence ladders, no corpus needed |
| Physics self-checks — reciprocity, energy balance | done |
| Convergence harness ([`convergence.py`](../src/gratinglab/convergence.py)) | done |
| Progress + cancellation on the `Solver` protocol | done |
| Mutation testing ([`mutation-testing.md`](mutation-testing.md)) | done, on demand; `checks.py` and `geometry.py` at 100% |
| Materials — CXRO constants, Fresnel, roughness, absolute efficiency | done |
| Resolving power ([`resolution.py`](../src/gratinglab/resolution.py)) | done: `R = |m|·N` plus the finite-N line profile from ISSI eq. (8), each testing the other; measured `N` flows in from metrology via `to_problem`. Groove-error-degraded R is the named successor |
| **Tip deconvolution** ([`core/tip.py`](../src/gratinglab/metrology/core/tip.py), [theory §4b](theory/metrology.md)) | done at M22: Villarrubia erosion with the certainty map, opt-in via `tip_correction='erosion'`, validated pointwise on a committed dilated fixture. Ran on the real TASTE scan and narrowed the rounded-groove open question: a nominal tip is excluded ([`findings.md`](findings.md)) |
| GUI (Qt/PySide6) — tabs, geometry dock, 3D conical view | done |
| CI (Linux + macOS, py3.11/3.12) | green |
| RCWA | **contributed, not yet integrated** |
| C-method | not started |
| **Integral method** ([theory](theory/integral.md)) | **perfect-conductivity and finite-conductivity milestones done**: conical from day one, undercut-capable, GUI tab, validated against M&P Table 4.1 and orderwise against PCGrate to 8e-4 (which also identified the TASTE profile). M20 added `conductivity="tabulated"` — the coupled conical system of Goray & Schmidt (2010), absolute efficiencies with absorption recorded, flat limit = Fresnel to 6e-5, their Table 3 reproduced per order, and both boundary conditions selectable in the GUI. Graded corner mesh is the named successor |

Roughly 1300 tests, 11 skipped (the skips need the private reference corpus).
Rounded on purpose: this line previously claimed 328 and was wrong by 400, so
`pytest -q` is the authority and this is only the order of magnitude. The jump
from ~945 is the metrology suite arriving with its package, not new coverage.

---

## What is next, in order

### 0. Integral method — pulled ahead of RCWA, deliberately

An earlier revision of this file said *"do not attempt the integral method
before the first release"* and slotted it as phase 5. That ordering was
reversed, for three reasons the earlier revision could not have known:

1. **The corpus can only validate an integral solver.** The one rigorous
   reference dataset on disk is perfect-conductivity *integral-method* output
   (PCGrate). Scalar could never test more than the diffraction pattern
   against it; a like-for-like comparison needed a like method.
2. **The off-plane X-ray application is the project's point, and the integral
   method is the method of record there.** The conical decoupling theorem
   makes the perfectly conducting milestone *cheaper* at grazing incidence,
   not harder — the reduced wavelength λ/sin γ is what the mesh resolves.
3. **RCWA is contributed work, blocked on a contributor.** The
   second-implementation stress test of the solver protocol (Capabilities,
   progress, tabs, convergence knob) happened with this backend instead — and
   the protocol held; see the current-state table.

The gamble paid off immediately: the first rigorous-vs-rigorous comparison
reproduced PCGrate orderwise to 8e-4 *and* identified the TASTE run's
previously unconfirmed groove profile ([`findings.md`](findings.md)).

Remaining integral milestones, in order:

- ~~**Finite conductivity**~~ — done at M20: the coupled conical system of
  Goray & Schmidt (2010), `conductivity="tabulated"`, absolute efficiencies
  with the absorbed fraction recorded and `R + A = 1` a checkable theorem
  ([`theory/integral.md`](theory/integral.md) §8). The Leontovich rung was
  deliberately skipped — it requires |n| ≫ 1 and is invalid in the
  soft-X-ray regime the milestone exists for, and it builds no operators the
  full system reuses. Validated without a PCGrate reference run (see the
  narrowed open question below): flat interface = Fresnel exactly, the
  energy theorems two-sided, the perfect-conductor limit against the
  already-banked 8e-4 PCGrate agreement, and Goray & Schmidt's own Table 3
  reproduced within 6e-5 per order — which also exposed their Table 4 as a
  publisher's duplication ([`findings.md`](findings.md)).
- **Graded corner mesh** (M&P §4.6.6) — TM on cornered profiles currently
  converges first-order, and with finite conductivity both polarizations
  do; see [`theory/integral.md`](theory/integral.md) §5.
- ~~A solver-contract conformance test~~ — done, see the note under RCWA below.

### 1. RCWA — integrating a contributed backend

**Not ours to write.** `fgrise` has RCWA code. The physics list below is
therefore a set of requirements to *state to a contributor* rather than a plan
of work, and the more useful half of this section is what any backend has to
satisfy to register at all.

Physics, unchanged and non-negotiable:

- Li's Fourier factorization rules from day one, or TM metallic cases converge
  badly and the comparison plots will misrepresent the method.
- S-matrix / enhanced transmittance-matrix propagation, never plain T-matrix.
- Full vectorial conical formulation — off-plane is not an afterthought.

What the framework asks:

| Requirement | Why, and what reads it |
|---|---|
| A `Capabilities` declaration | `conical`, `polarizations`, `handles_undercut`. `Capabilities.check` refuses an out-of-scope case *before* solving rather than approximating it |
| `accuracy_knob` naming the truncation parameter | `convergence.check_convergence` sweeps it. Without one, a result can never report `converged=True` |
| `reports_progress`, if `solve` takes a `progress` keyword | Opt-in, so a backend predating the hook is never handed a keyword its signature cannot take. It also makes the solver **cancellable** — that callback is the only check point a NumPy-bound solve has |
| `UnsupportedConfiguration` rather than a quiet approximation | The rule the whole comparison rests on: a method silently smoothing a vertical facet returns a plausible, wrong number |
| One entry in `_TAB_FACTORIES` | Tabs are generated from `available_solvers()`; the GUI needs nothing else |

**The seams have now met their second implementation** — the integral
backend. What the earlier warning here predicted mostly did not happen: the
`Capabilities` declaration, the progress/cancellation contract, the
convergence knob, and the per-solver tab pattern all took a slow, rigorous,
polarization-resolving backend without changing shape. Two things did give:
the window's "solve on open" needed to prefer the fast solver once a
minutes-long one existed, and the tab layer is now two near-identical files
— the evidence the base-class extraction was waiting for.

A conformance test every registered solver must pass is now
[written](../tests/test_solver_conformance.py) and parametrised over
`available_solvers()`, so a contributed backend is held to the declaration,
refusal, scan and progress contracts the moment it registers — before anyone
writes a line of test for it specifically.

### 2. Resolving power — the ideal half landed, the measured half named

`gratinglab.resolution` now answers the ideal-grating question: `R = |m|·N`
by the Rayleigh criterion, and the finite-N line profile it is derived from
(`interference_factor`, moved from the scalar solver to `geometry` where
non-solver consumers can reach it). The measured groove count arrives through
`BoundaryProfile.to_problem`, so an AFM scan flows to a spectrograph number
in-process. A problem without `n_grooves` is refused, not treated as infinite.

The successor is **groove-error-degraded resolving power**: the line-spread of
a real groove ensemble — period jitter and placement errors, the quantities
metrology already half-holds (`y_std_nm`, `groove_centers`) — with the
Rayleigh criterion applied to the degraded profile. That is the link between
groove *error* metrology and spectrograph performance, and it is deliberately
not approximated by the ideal formula now in the API.

### 3. PyXFocus analytics — port the formulas, not the raytracer

PyXFocus (the sequential X-ray raytracer this project's author forked) stays
where it is: unpackaged, numpy<2-pinned, Fortran-built — absorbing it would
buy a liability, not a capability. What is worth porting is small and
analytic, from its `grating.py`:

- `blazeYaw` — the Littrow condition, as `geometry.littrow_alpha` plus an
  `Illumination.littrow(order, wavelength, period, graze)` constructor.
- `blazeAngle`/`blazeAngle2` — the facet angle Littrow requires at (λ, m).
  First *verify numerically that the two agree* and pin the frame mapping
  (its `cos inc` ↔ our `sin γ`) against `Illumination.offplane`, the
  §10-of-conventions discipline: transcribe, test, then port one form.
- `litBetaAlpha`, `eta`, `yaw` — already representable with `Illumination`
  machinery; cover with tests rather than porting.
- `radialApproxEffect` — not ported: radial/variable-line-space physics that
  `Problem` cannot represent (and a live `pdb.set_trace()` in the source).
  Radial/VLS gratings are a "Later" item.

Its raytraced resolving-power estimator (order-centroid spacing over spot
width) is the instrument-level cross-check of `resolution.py`'s grating-level
number; that comparison belongs with any future instrument layer, not here.

### 4. ~~AFM tip deconvolution~~ — done at M22, and the arbiter ruled

Villarrubia (1997) morphological erosion **with the certainty map**, exactly
as planned: a parametric cone-plus-spherical-apex tip, a stage in
`metrology/core/tip.py` on the 2-D array after image flattening in both
branches, opt-in through `AnalysisSettings` (`tip_correction`,
`tip_radius_nm`, `tip_half_angle_deg`), the certain fraction recorded in the
metrics, the sidecar and the log. Validated on a committed dilated fixture —
the sharp synthetic scan imaged through a known worn tip with the *identical*
noise field, so recovery is asserted pointwise: machine precision on clean
data, the noise floor (0.18 nm mean) on the fixture, against a 9 nm mean
uncorrected error. [`theory/metrology.md`](theory/metrology.md) §4b carries
the three load-bearing statements (a bound, not a resurrection; the certainty
map is the deliverable; a facet steeper than the flank is gone outright).

The run on the real TASTE scan produced the milestone's finding
([`findings.md`](findings.md), "A nominal tip does not explain the rounded
groove"): erosion with the nominal probe (R ≈ 1–2 nm, 18°) recovers **no**
depth — the measured surface is already tip-reachable — and reproducing the
depth deficit forward requires an apex worn to ~80 nm. The open question
narrows from "grating or tip?" to "grating, or a severely worn tip?", which
is answerable outside the software: the probe's condition at scan time, a
re-scan with a fresh tip, or a cross-section.

What deliberately did not land: measured tip images as input (blind tip
estimation is the natural successor if a worn tip becomes the suspect), and
any coupling from the certainty map into `to_problem` — an uncertain-trough
profile still converts, and the defensibility warning remains the
depth-vs-facet-angle diagnostic named in the finding.

### 5. Native boundary format

`.ggp` cannot carry the period in nm, provenance, undercut boundaries, or a
format version. The missing period is exactly why
[`benchmarks/corpus.toml`](../benchmarks/corpus.toml) has to exist as a sidecar.
JSON, since the profile classes are already pydantic models with tested
round-tripping.

**M17 changed the urgency, not the need.** With metrology in the package, a
measured groove reaches a solver through `BoundaryProfile.to_problem` without
touching disk, carrying its measured period and optionally a fitted facet
angle. That is the whole gap closed — *for profiles that never leave the
process*. Everything that has to be saved, sent, archived or diffed still goes
through a format that drops the period on the way out and cannot say where the
profile came from. Undercut boundaries remain unrepresentable in either
direction.

### 6. First release + JOSS

Shipping something citable and correct early is what recruits collaborators.
The release story is stronger now, not weaker, for the reordering above: a
rigorous open solver validated orderwise against the commercial standard is a
headline, and the scalar/integral comparison plots need no PCGrate licence.

### Later

C-method, measured-profile fitting, crossed/2D gratings, radial/VLS gratings.
See [`references.md`](references.md) for the literature mapped to each.

#### Measured efficiency — the als_eff_meas split

The sibling `als_eff_meas` package reduces ALS beamline scans to absolute
measured efficiency per order, uncertainties carried end-to-end. It is the
measured-data half this platform still lacks — and it is also genuinely
specific to X-ray/EUV testing at one facility, so the integration is a
*split*, not an absorption:

- **Moves here, when taken up:** a measured-efficiency data model
  (`MeasuredEfficiency`, a new `measured.py`) that `compare.align`/`records`
  accept beside `EfficiencyScan` — measured vs scalar vs integral on one
  plot, with error bars. *Not* a second `EfficiencyTable` (that name is taken
  by the modelled-table importer in `io/efficiency_table.py`), and *not*
  uncertainty fields grafted onto `EfficiencyScan`, whose dense-grid
  invariants are model-shaped; measured data is sparse per order with
  per-point σ. Parallel `value`/`sigma` float arrays, so the `uncertainties`
  package never becomes a dependency here. Also the PCGrate **XML** dialect
  from its `models.py`, which `io/efficiency_table.py` (text dialects) lacks.
- **Stays there:** everything beamline-shaped — file pairing, dark
  subtraction, three-point flux sums, the ODR arc fit and cone-angle
  recovery, TOML campaigns. `als_eff_meas` grows a gratinglab dependency,
  deletes its duplicated Fresnel/Névot–Croce/CXRO code in favour of
  `gratinglab.materials`, and exports a `to_measured()`.
- **The scipy pin stays out:** `als_eff_meas` pins `scipy<1.19` because
  `scipy.odr` backs its arc fit and is scheduled for removal. The arc fit
  does not move, so the pin does not propagate — and that package must
  eventually migrate to `scipy.optimize.least_squares` regardless; this
  project should not inherit the debt in the meantime.

The other siblings are settled: `afm_blaze_meas` is already absorbed (M17),
`endstation` stays the shared scaffolding dependency, `ebl_fracture` stays a
separate fabrication-layout tool (no optics; its own roadmap excludes blazed
and VLS geometry), and PyXFocus is §3 above — formulas, not the raytracer.

---

## Open questions

Answers to these unblock real work.

### The panter1 grating period

Recovered from the reference data with no assumptions:

- `period × sin γ = 9.8450 nm` (±0.005)
- `α = +32.92°` (±0.1)

These two cannot be separated without one more fact. Candidates:

| γ | period |
|---|---|
| 1.5000° | 376.1 nm |
| 1.7902° | 315.15 nm (same grating as TASTE) |
| 2.0000° | 282.1 nm |

The title string in the original plotting script reads *"η = 1.5° & φ = 0.97°"*.
The 1.5° is consistent with γ if the period is ~376 nm, but **0.97° is not the
azimuth** — the data says α ≈ 33°, which is also far more plausible given the
~39°/27° facet angles of `panter1.ggp` (near-Littrow). So "φ" there means
something else, and knowing what would help read the other legacy scripts.

### Which `.ggp` goes with which run

`AFM_TASTE_test.ggp` and `TASTE_a205_TiPt.ggp` are both candidates for the TASTE
efficiency run — and **neither is present in the corpus directory**, which holds
`AFM_test.ggp`, `AFM_real_echelle.ggp` and `PCGrateProjects/panter1.ggp`. The
filename in `corpus.toml` therefore does not resolve, which is itself a thing to
settle. Until it is, the comparison uses an idealised `Blazed(29.5°, 70.5°)`,
which is probably why scalar's total sits near 0.55 rather than closer to unity.

M16 remeasured this: the total is 0.5431 for the idealised sawtooth and 0.6138
for `AFM_real_echelle.ggp`, against a PCGrate total of 1.0005. The flux
obliquity of M16-C moved it in the right direction but only slightly, because
this mount is near enough to Littrow that the factor is ~0.99 on the orders
carrying the power. The residual is profile mismatch, not normalisation — see
`findings.md`, "The corpus can test the diffraction but not the reflectivity".

**M17 tested that and it failed.** The metrology merge brought an AFM scan of
what is plausibly the same grating into the repo
(`data/TASTE_ALS_A205_Ti_Pt_flatten.txt`), so the idealised sawtooth could be
swapped for a measured groove. The total moved the *wrong* way: 0.5420 idealised
to 0.3362 measured, against PCGrate's 1.0005. The groove's facet fit (27.91° ±
2.13°) and its depth (implying 20.33°) disagree by 3.6σ: it is rounded rather
than faceted, and the solver reads the depth. A flat land would explain the same
deficit and was excluded (it would need 25–30% of the period; 2.9% is within 2°
of flat), as was smearing from the groove averaging (0.4%). Full evidence in
`findings.md`, "The measured groove is rounded, not faceted".

So the open question is no longer only "which `.ggp`". Whether the rounding is
the grating or the tip is undetermined, and either way an absolute efficiency
from a raw AFM boundary is not defensible yet. Nothing in the pipeline warns
about it. The arbiter has now run (§4 above, M22): erosion with the nominal
tip explains **none** of the 27.91°/20.33° gap, and reproducing it forward
needs an ~80 nm apex. The question is narrowed, not closed — "grating, or a
severely worn tip?" — and what settles the remainder is outside the software:
the probe's condition at scan time, a re-scan with a fresh tip, or a
cross-section. See `findings.md`, "A nominal tip does not explain the
rounded groove".

### A finite-conductivity reference run

Narrowed at M20, not closed. The corpus still has no usable
finite-conductivity table (`panter1_finite` stays `usable = false`), but the
rigorous solver no longer waits on one: `conductivity="tabulated"` is
validated on the Fresnel flat limit, the two-sided energy theorems, the
perfect-conductor limit, and Goray & Schmidt's published Table 3. What a
fresh PCGrate run would still buy is a *like-for-like corpus item* — the
recommended recipe: TE, Au, the TASTE geometry and the identified
`AFM_test.ggp` groove, a wavelength grid avoiding ±0.05 nm around every
order-passing-off point (the documented failure locus of `panter1_finite`),
and the same `.ari` optical constants exported by `write_ari` so both codes
eat identical n(λ). Gate it with `check_energy_balance` before admitting it
to `corpus.toml`. It would also finally arbitrate the scalar solver's
groove-resolved reflectivity model, which still rests on internal
consistency plus (now) cross-checks against the rigorous absolute mode.

### Smaller ones

- Is the vendor console solver runnable (Windows VM)? Generating fresh reference
  cases on demand, rather than only reusing exported tables, is worth an
  afternoon.
- Any synchrotron or Panter measured efficiency beyond the current corpus that
  could seed the benchmark suite? Partially answered: `als_eff_meas` holds one
  fully ported ALS campaign (ogre laminar, Oct 2021) with absolute
  per-order efficiencies and uncertainties, plus ~15 legacy per-experiment
  scripts in its archive that could be ported the same way. Reaching it from
  the comparison harness is the measured-efficiency split under "Later".

---

## Practices worth keeping

- **Refuse rather than approximate.** A solver asked for something outside its
  declared capabilities raises. Silently approximating is what makes
  cross-method plots lie.
- **Report, never rescale.** Where an approximate theory violates a conservation
  law, that is information. Normalising it away destroys the validity map.
- **Provenance is mandatory.** A result whose convergence was never demonstrated
  reports itself as not defensible — and `gratinglab.convergence` is how it
  stops being that. `converged=False` is a result, not a missing one.
- **A control that cannot do what it says must not be offered.** Cancel stops
  the solver now; before it could only stop the waiting, and the panel said so
  rather than implying more. A backend declaring no `reports_progress` is
  simply not cancellable, and that is a fact to surface, not to paper over.
- **One agreement is not a plateau.** Refinement error is not always monotone;
  the harness requires three consecutive knob values to agree, because the
  blazed case measurably violates the tempting one-agreement rule
  ([findings](findings.md#quadrature-error-is-not-monotone-in-the-number-of-points)).
- **Commit before mutating.** Mutation scripts revert with `git checkout --`,
  which silently fails on untracked files.
