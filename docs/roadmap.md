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
| GUI (Qt/PySide6) — tabs, geometry dock, 3D conical view | done |
| CI (Linux + macOS, py3.11/3.12) | green |
| RCWA | **contributed, not yet integrated** |
| C-method | not started |
| **Integral method** ([theory](theory/integral.md)) | **perfect-conductivity milestone done**: conical from day one, undercut-capable, GUI tab, validated against M&P Table 4.1 and orderwise against PCGrate to 8e-4 (which also identified the TASTE profile). Finite conductivity is the named successor |

Roughly 1150 tests, 11 skipped (the skips need the private reference corpus).
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

- **Finite conductivity** — the coupled conical system of Goray & Schmidt
  (2010); optionally the impedance (Leontovich) boundary condition as an
  intermediate rung, which is also what finally tests the
  "ΣE ≈ R_F(ζ) is emergent" reading of `conventions.md` §10 item 5. This is
  what makes soft-X-ray efficiencies *absolute*. Blocked in part on the
  "finite-conductivity reference run" open question below.
- **Graded corner mesh** (M&P §4.6.6) — TM on cornered profiles currently
  converges first-order; see [`theory/integral.md`](theory/integral.md) §5.
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

### 2. Native boundary format

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

### 3. First release + JOSS

Shipping something citable and correct early is what recruits collaborators.
The release story is stronger now, not weaker, for the reordering above: a
rigorous open solver validated orderwise against the commercial standard is a
headline, and the scalar/integral comparison plots need no PCGrate licence.

### Later

Integral finite conductivity (see §0), C-method, measured-profile fitting,
crossed/2D gratings. See [`references.md`](references.md) for the literature
mapped to each.

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
about it.

### A finite-conductivity reference run

The corpus has no usable one. Every perfect-conductivity table has R ≡ 1 by
construction, so nothing in it can validate the reflectivity model that M16-D
made the default, and `panter1_finite` is already marked `usable = false`.
Either a fresh finite-conductivity PCGrate run or the RCWA backend would close
the gap; until one of them lands, the groove-resolved reflectivity rests on
internal consistency alone.

### Smaller ones

- Is the vendor console solver runnable (Windows VM)? Generating fresh reference
  cases on demand, rather than only reusing exported tables, is worth an
  afternoon.
- Any synchrotron or Panter measured efficiency beyond the current corpus that
  could seed the benchmark suite?

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
