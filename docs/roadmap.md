# Roadmap

Where this project is, where it is going, and why. Start here.

---

## The goal

An open comparison platform for grating efficiency analysis. Not "another RCWA" —
open RCWA is a crowded field (S4, grcwa, nannos, torcwa, inkstone, RETICOLO).
What does *not* exist openly is:

1. A neutral cross-method harness with one problem spec and rigorous conventions.
2. An open C-method (Chandezon) solver.
3. An open integral-method solver.

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
| **Scalar solver** ([theory](theory/scalar.md)) | done, validated |
| Efficiency-table importer (external codes) | done |
| `.ggp` boundary profiles | done |
| Comparison harness | done |
| Physics self-checks — reciprocity, energy balance | done |
| Convergence harness ([`convergence.py`](../src/gratinglab/convergence.py)) | done |
| Progress + cancellation on the `Solver` protocol | done |
| Mutation testing ([`mutation-testing.md`](mutation-testing.md)) | done, on demand; `checks.py` and `geometry.py` at 100% |
| Materials — CXRO constants, Fresnel, roughness, absolute efficiency | done |
| GUI (Qt/PySide6) — tabs, geometry dock, 3D conical view | done |
| CI (Linux + macOS, py3.11/3.12) | green |
| RCWA | **contributed, not yet integrated** |
| C-method | not started |
| Integral method | not started |

Roughly 945 tests, 11 skipped (the skips need the private reference corpus).
Rounded on purpose: this line previously claimed 328 and was wrong by 400, so
`pytest -q` is the authority and this is only the order of magnitude.

---

## What is next, in order

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

**Worth saying plainly, because the first contributor is the one who finds
out:** none of those seams has met a second implementation. Each was designed
against exactly one solver, and that solver is fast, closed-form and
single-polarization — about as unlike RCWA as a backend can be. Expect the
protocol to be wrong somewhere, and read the first integration as evidence
about the protocol, not only about the contribution.

A conformance test every registered solver must pass is the natural follow-up.
It belongs *with* that work rather than before it: there is nothing to run it
against until a second solver exists.

### 2. Native boundary format

`.ggp` cannot carry the period in nm, provenance, undercut boundaries, or a
format version. The missing period is exactly why
[`benchmarks/corpus.toml`](../benchmarks/corpus.toml) has to exist as a sidecar.
JSON, since the profile classes are already pydantic models with tested
round-tripping.

### 3. First release + JOSS

Shipping something citable and correct early is what recruits the collaborators
who make the integral method tractable. **Do not attempt the integral method
before the first release.**

### Later

C-method (phase 3), measured-profile fitting (phase 4), integral method
(phase 5), crossed/2D gratings (phase 6). See
[`references.md`](references.md) for the literature mapped to each.

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
efficiency run. Until this is settled the comparison uses an idealised
`Blazed(29.5°, 70.5°)`, which is probably why scalar's total sits near 0.55
rather than closer to unity.

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
