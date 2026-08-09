# Findings

Empirical results established while building this, each with the evidence that
established it. Recorded so they are not rediscovered, and so a claim can be
checked rather than taken on trust.

Findings that constrain how code must behave are also encoded as tests; those
are cross-referenced.

---

## Quadrature error is not monotone in the number of points

**Refining a blazed-grating quadrature can make the answer worse, and a
convergence rule that trusts a single agreement will certify a tolerance the
next refinement violates.**

Scalar solver, reference geometry (period 315.15 nm, δ = 29.5°, anti 70.5°,
α = 25°, γ = 1.5°), λ ∈ [1, 5] nm. Largest change in efficiency per doubling of
`quadrature_points`:

| n | change vs n/2 | ratio to previous change |
|---:|---|---:|
| 256 | 5.613e-04 | 22.8× smaller |
| 512 | 2.265e-04 | 2.5× smaller |
| 1024 | 9.645e-06 | 23.5× smaller |
| **2048** | **2.426e-05** | **2.5× LARGER** |
| 4096 | 4.835e-06 | 5.0× smaller |
| 8192 | 5.964e-07 | 8.1× smaller |
| 16384 | 1.567e-07 | 3.8× smaller |
| 32768 | 3.527e-08 | 4.4× smaller |

Past ~8192 the ratio settles near 4, i.e. second-order convergence — the
expected rate for a rule integrating a function with a slope discontinuity.
Below that it swings between 0.4 and 23.

**The cause is geometric, not numerical noise.** The blazed profile has a kink
at `t = 1/(1 + tan δ / tan δ') = 0.8331`. That is not a dyadic rational, so
doubling `n` changes how near the closest node falls to the kink, and the local
error at the kink dominates. A sinusoid, being analytic, shows none of this: it
is at machine precision (~1.7e-16) from `n = 64` upward and never moves again.
A lamellar profile — a jump, not a kink — converges cleanly at exactly 4× per
doubling with no excursions.

**Consequence, encoded in code.** `gratinglab.convergence` requires a *plateau*:
`DEFAULT_PLATEAU = 2` consecutive differences below tolerance, so three knob
values in a row must agree. At a tolerance of 1e-5 the naive one-agreement rule
stops at n = 512 on the strength of that 9.6e-6, and the next refinement then
moves the answer by 2.4e-5. One extra solve rejects it.

This is why the harness reports the whole ladder rather than a boolean, and why
`converged_at` is the coarsest value of the plateau rather than the first value
that happened to agree.

Pinned by `tests/test_convergence.py::TestTheMeasurementBehindTheRule`, which
fails if refinement ever becomes monotone here — at which point this entry and
the plateau argument would both need revisiting — and by
`TestThePlateauRule::test_and_plateau_1_would_have_stopped_at_the_dip`, which
shows the naive rule failing on a scripted sequence with the same shape.

---

## Corpus geometry recovered from the data alone

**The exported efficiency tables record efficiencies but not the geometry that
produced them, and the vendor project files are binary.** `strings` on them
yields only a build stamp and a path to an optical-constants file.

The geometry is nonetheless recoverable, exactly. Which orders propagate is
fixed by the grating equation, and the constraint **factors**:

$$m_{\max} = \left\lfloor A/\lambda \right\rfloor, \quad
  m_{\min} = -\left\lfloor B/\lambda \right\rfloor$$

$$A = (1 + \sin\alpha)\,p\sin\gamma \qquad B = (1 - \sin\alpha)\,p\sin\gamma$$

$A$ and $B$ decouple, so each is a simple interval intersection over the scan.
Across a 500-point wavelength scan this pins both to about ±0.01 nm.

**Validation.** For the TASTE run the period is independently known —
315.15 nm, hard-coded in the group's own plotting script. Feeding that in gives

$$\gamma = 1.2500°, \qquad \alpha = +19.99°$$

Both clean round numbers, which a mis-derived method would not produce.

A test asserts that our own grating equation, driven by the recovered geometry,
regenerates the exact propagating-order set of the reference file at every
sampled wavelength. It does.

Results live in [`benchmarks/corpus.toml`](../benchmarks/corpus.toml). The
panter1 period remains unresolved — see [roadmap](roadmap.md#the-panter1-grating-period).

> **Sign trap.** The first pass assigned $A$ and $B$ the other way round and got
> α = −19.99°. It looked plausible. What caught it was running the scalar solver
> at both signs: mean summed efficiency was 0.23 at −19.99° versus 0.55 at
> +19.99°. A sign error in a recovered quantity does not announce itself.

---

## A reference file violates energy conservation

`panter1_fix_efftable_finite.txt` — a finite-conductivity run — has **summed
efficiency reaching 3.6**. Verified against the raw file, so it is a defect in
the run, not in the parser. The companion perfect-conductivity run of the same
problem conserves energy to four decimal places, so the setup is sound; it is
the finite-conductivity solve that destabilises.

79 of 561 points (14 %) violate, in five clusters:

| λ range (nm) | points | orders | peak Σℰ |
|---|---|---|---|
| 0.60–0.89 | 30 | 23–33 | 2.94 |
| 0.97–1.03 | 7 | 19–20 | 1.23 |
| 1.11–1.12 | 2 | 18 | 3.60 |
| 1.44–1.49 | 6 | 14 | 2.03 |
| 4.16–4.49 | 34 | 5 | 3.56 |

**The cause is not propagating-order count.** The last cluster carries only five
orders. What predicts the violations is **order passing-off (Rayleigh
anomalies)**: violating points sit a median 0.020 nm from a wavelength where the
propagating-order count changes, against 0.125 nm for clean points, and that
last cluster terminates exactly at the 5→4 passing-off wavelength. The gold
optical constants are smooth across all five bands, so no material edge is
involved.

**Do not use this file as reference data at those wavelengths.** Encoded as
`test_known_defect_finite_conductivity_run_is_unphysical`; if the run is ever
regenerated correctly that test fails and should be deleted.

---

## Scalar theory cannot be reciprocal *and* conserve energy

No scalar formulation satisfies all three of the properties below. This is
structural, not incidental: reciprocity requires the phase to be symmetric under
$\alpha \leftrightarrow \beta_m$, and Parseval requires a single fixed function
that cannot contain $\beta_m$.

| Property | order-dependent $\varphi = kg\sin\gamma(\cos\alpha + \cos\beta_m)$ | $\beta$-free $\varphi = 2kg\sin\gamma\cos\alpha$ |
|---|---|---|
| Energy, $\Sigma \leq 1$ | ✗ up to 1.71 | ✓ 1.0000 exactly (Parseval) |
| Reciprocity $E_m(\alpha) = E_m(\beta_m)$ | ✓ $1\times10^{-17}$ | ✗ 0.44 violation |
| Blaze direction | ✓ exact at all $\alpha$ | ✓ at Littrow only; 7.6° off at 10°, evanescent past ~15° |

Physical optics is reciprocal but not energy-conserving; the pure transmittance
picture is the reverse. **Reciprocity is what the solver keeps.** The
energy-conserving variant is not offered as an option, because silently trading a
well-understood energy defect for a reciprocity violation is a footgun.

Referencing the phase to the blaze direction $\beta_b = 2\delta - \alpha$ looks
like a third way out and is not one: $\beta_b$ is an *output* of the model — it
is where the $\operatorname{sinc}^2$ envelope peaks — so feeding it back in is
circular, and the "zero error at every $\alpha$" it produces is a tautology.

### The implementation is not at fault

In the shallow-groove limit the sum is exactly unity, degrading smoothly with
depth:

| depth / period | $\Sigma_m E_m$ |
|---|---|
| 0.0001 | 1.00000 |
| 0.10 | 0.90472 |

The deviation scales with **phase excursion across the groove** — depth relative
to wavelength, hence working order — not with $\lambda/p$ or the number of
propagating orders. That is why $\Sigma$ looked stubbornly constant (~0.9065)
across wildly different geometries at fixed profile.

### Appendix D's two extra factors are both superseded

The thesis carries $\mathcal{E}_n \propto [\cos\beta_n/\cos\alpha]\|G_n\|^2$ and
then renormalises so $\sum_n \mathcal{E}_n = 1$. Both are **removed**, not made
optional. For $N \to \infty$ grooves the efficiency is the norm-squared Fourier
coefficient and nothing else; ISSI §2.1's groove function
$f(x) = \exp\{ik[n(k)-1]y(x/p)\}$ carries neither factor. The renormalisation is
the more damaging of the two — forcing $\Sigma$ to unity would erase exactly the
validity signal described above.

Full derivation: [theory/scalar.md §5](theory/scalar.md).
The deviation is reported as a provenance warning and never rescaled.

---

## Mutation testing found a real gap

Deliberately corrupting the physics and checking whether the suite noticed.
First sweep: seven mutations, six caught, **one survivor**.

### The obliquity factor was unverified

Inverting $\cos\beta_m/\cos\alpha \to \cos\alpha/\cos\beta_m$ passed the entire
suite. The only test asserted the result *changed* when the option was enabled,
never that it changed **correctly**.

That is exactly the thesis-versus-ISSI discrepancy recorded in
`conventions.md` §5 — a place where being wrong was plausible rather than
hypothetical. It was pinned to the exact ratio, with a second test proving the
two ratios are distinguishable for the geometry used, so the assertion could not
go vacuous.

> **Since superseded.** The obliquity option was removed entirely once the
> factor was established to be an error rather than a convention — see *Scalar
> theory cannot be reciprocal and conserve energy* above. The lesson stands and
> is the reason that section could be written with confidence: a test that only
> asserts *something changed* verifies nothing.

### In-plane tests are blind to γ errors

Breaking $\sin\gamma \to \sin^2\gamma$ left **34 of 37** scalar tests passing.
In-plane cases have $\sin\gamma = 1$, which makes any power of it identical, so
they cannot detect an error in the γ handling at all. Only the off-plane cases
failed — and off-plane is the primary application.

Coverage at the time was 9 off-plane tests against 13 in-plane. Every
closed-form check now runs in **both** mounts.

After both fixes, all 12 mutations are caught.

### Then a real tool ran, and found two more

The hand-rolled sweep was seven mutations chosen by someone who knew where to
look. `mutmut` generated **19 for `geometry.beta` alone**, and four survived.
Two were real:

- **`np.pi - b` → `np.pi + b` for transmitted orders.** `beta(..., transmitted=True)`
  had no test at all — only `cos_beta` did. The mutant negates \(\sin\beta\)
  while leaving \(\cos\beta\) negative, which is precisely the thing
  `conventions.md` §4 forbids ("we do not flip the sign of \(\sin\beta_m\)
  for transmission"). Now pinned by asserting both halves of that sentence in
  one test.
- **`|s| <= 1.0` → `|s| < 1.0`.** At exactly \(|\sin\beta_m| = 1\) — the
  passing-off boundary — the mutant returns NaN. `is_propagating` uses `<=`, so
  an order it counts as propagating would have had no direction to be drawn or
  summed at. The two functions now have to agree at the boundary, with a
  companion showing the boundary is a boundary.

The other two are **equivalent mutants**, verified rather than assumed:

- `<= 1.0` → `<= 2.0` is unobservable, because `arcsin` of anything outside
  \([-1, 1]\) is NaN anyway and `errstate(invalid="ignore")` suppresses the
  warning either way. The guard is doing less than it looks.
- Dropping `dtype=np.float64` from `np.asarray` changes nothing NumPy does not
  already do by promotion.

No test can kill either, and chasing them would mean writing assertions about
nothing. Recorded here so the next reader does not try.

**The lesson, sharpened.** The earlier sweep's conclusion was "a test that only
asserts something changed verifies nothing". This one adds: *a branch with no
test at all is invisible to code coverage when another function exercises the
same line*. `cos_beta(0.5, transmitted=True)` covered the concept; nothing
covered `beta`'s version of it.

---

## Reciprocity is the sharpest available check

Lorentz reciprocity requires $\mathscr{E}_m(\alpha) = \mathscr{E}_m(\beta_m)$.
Scalar theory satisfies it exactly, because the phase
$\cos\alpha + \cos\beta_m$ is symmetric under the exchange.

| | max violation |
|---|---|
| as implemented | 5 × 10⁻¹⁸ |
| with the α↔β symmetry of Φ broken | 4 × 10⁻¹ |

Sixteen orders of magnitude, requiring **no reference data and no closed form**.

This matters because the closed-form tests cannot do it. They compare the solver
against a formula derived the same way the solver computes it, so they validate
the *quadrature*. Reciprocity constrains the *structure* of the phase function.

---

## The cone opens away from the groove-axis vector

Three individually true statements in `conventions.md` §3, whose conjunction is
contradictory under a directed reading:

| Where | Statement |
|---|---|
| frame table | `ĝ = -ẑ` |
| wave vectors | every `k` carries `+\cos\gamma\,\hat{z}` |
| angles | γ is "the *polar* angle measured from the groove axis" |

Measured directly at the reference geometry (α = 25°, γ = 1.5°): the angle from `+ẑ` is
**1.5000°**, and the angle from `ĝ` is **178.5000°**. Every propagating order sits at
exactly 1.5000° from `+ẑ` with an identical `k_z = 0.999657325` — they genuinely lie on
one cone, and **that cone opens along `-ĝ`.**

$$\hat{d} \times \hat{g} = \hat{n} \;\Longrightarrow\; \hat{x} \times \hat{g} = \hat{y}
\;\Longrightarrow\; \hat{g} = -\hat{z}$$

so `ĝ = -ẑ` is forced by the right-handedness relation — which comes from Heilmann et al.
(2024) §2.1 — and is not a free choice. What was loose was the *wording* of the γ bullet:
γ is the acute angle to the groove-axis **line**, i.e. the polar angle about `-ĝ`.

**Nothing was ever wrong.** Every consumer of γ — `facet_graze`, `sin_beta`,
`blaze_wavelength`, `Illumination` — takes it as the scalars `sin γ` / `cos γ`, where the
direction of `ĝ` cannot enter. `Illumination.direction_cosines` gets the sign of its `z`
component right, but nothing asserted it *against* `ĝ`. **No code had ever needed a
directed groove-axis vector.**

Same shape as [the `t̂ = -d̂` finding](#the-profile-parameter-runs-backwards-and-nothing-had-noticed):
invisible until something drew it. The first thing to need a directed `ĝ` is a 3D view of
the diffraction cone, and getting it backwards would open the cone into the grating
instead of away from it — every ray 180° out.

Recorded rather than "fixed" by flipping `ĝ` to `+ẑ`, which would break `d̂ × ĝ = n̂` and
diverge from the reference. Pinned by
`tests/test_illumination.py::TestTheConeOpensAwayFromTheGrooveAxis` — placed where the bug
would have lived, not where it would have been noticed.

---

## The profile parameter runs backwards, and nothing had noticed

`geometry.blaze_direction` returns $\beta_b = 2\delta - \alpha$. That value requires the
active facet's outward normal at azimuth $+\delta$ — the facet *descending* with
increasing $x$. Reflecting the projected incident direction about the other orientation
gives $\sin\beta = -\sin(2\delta + \alpha)$ instead: at the reference geometry
($\delta = 29.5°$, $\alpha = 25°$) that is **−84.0°** where the blaze direction is
**+34.0°**.

But `Blazed.slope()` returns $+\tan\delta$ for $t \leq$ `apex` — the active facet *rises*
with $t$. Rising in $t$ while descending in $x$ means $\hat{t} = -\hat{d}$.

**Nothing was ever wrong.** The scalar solver's blaze peak lands where
`blaze_wavelength` predicts:

| order | predicted $\lambda_b$ | solver peak | ratio |
|---|---|---|---|
| 2 | 4.0498 nm | 3.9970 nm | 0.987 |
| 3 | 2.6999 nm | 2.6863 nm | 0.995 |

(The residual is the sinc envelope peaking between discrete orders, not a sign error.)
`scalar.py` computes $\int e^{+i\Phi_m(t)}e^{-2\pi imt}\,dt$ where §2/§3 imply
$e^{-i\Phi_m}$ under $t \leftrightarrow +x$; the substitution $t \leftrightarrow -x$ maps
one integral onto the **complex conjugate** of the other, and conjugation leaves
$\|G_m\|^2$ bit-identical. The codebase has been self-consistent all along.

**What made it invisible:** no code had ever drawn the groove profile and the diffracted
rays in the same frame. `ProfilePlotPanel` plots height against $t$ left-to-right, which
*implies* $t \leftrightarrow +\hat{d}$ without asserting anything physical, and every
other consumer takes $\|G_m\|^2$ where the handedness cancels. The geometry visualizer is
the first thing that has to place a facet relative to a ray — and getting it backwards
would have drawn the beam striking the anti-blaze facet with the blaze arrow leaving
through the back of the active one.

Recorded in [`conventions.md` §3](conventions.md) rather than repaired by flipping
`scalar.py`'s signs. That flip is equally modulus-preserving and every test would still
pass, but it would move the module's headline formula away from the form transcribed from
the ISSI chapter and Appendix D. **Neither reference states a handedness** — both quote
only $\|G_m\|^2$ — which is the finding rather than an oversight to correct.

---

## An exactly-correct solver produced an unusable report

CI failed on all three matrix jobs while the same commit passed locally.

The reciprocity check tracked its worst violation with a running maximum seeded
at 0.0, updated only on a strictly greater value. A solver that is **bitwise**
reciprocal therefore never recorded which order the measurement came from.
Locally, floating-point noise gave ~5 × 10⁻¹⁸ and the field was set; on CI the
BLAS returned exact zero and the report came back empty.

**The better the solver, the worse the report** — and it only appeared on a
different platform. Fixed by collecting violations and taking the maximum at the
end. The regression test constructs a deliberately exactly-reciprocal solver so
the zero path is exercised everywhere rather than depending on rounding.

---

## Two corrupt files in the legacy collection

- `OGRE/AFM_test.ggp` line 164 has a **missing x value** — a truncated write.
  The reader refuses it rather than silently dropping the point, which would
  shift the profile and change every efficiency computed from it.
- `OGRE/pcgratetest3.txt` has **no efficiency columns at all** — 20,308 rows of
  scan geometry and no data.

Both rejections are correct behaviour, and both are encoded as tests so the
rejections are not mistaken for parser bugs.

### Three `.ggp` header variants exist, and two do not load

| Header | Vendor accepts | Seen in |
|---|---|---|
| `3 0 - Polygonal type` + `Period:` | yes | most files |
| `# 3 0 - Polygonal type` (hashed) | **no** | `*_fix6/7.ggp` |
| `# X(normalized...)` only | **no** | `*_fix{,2,3,4}.ggp` |

`np.savetxt(header=...)` is what produces the hashed variants. We read all
three — they are real files someone needs to load — and write only the first,
flagging the others via `GgpFile.format_valid`.

---

## Thesis and ISSI chapter reconciled

Cross-checking the two primary references: the grating equation, blaze direction
$\beta_b = 2\delta - \alpha$, blaze wavelength, and facet graze angle
$\zeta = \arcsin[\sin\gamma\cos(\delta-\alpha)]$ **all agree exactly**.

Three differences change what the code computes, and four are transcription
errors. All seven are tabulated in [`conventions.md` §10](conventions.md).
The two that matter most:

- The thesis sawtooth phase omits `sin γ` and is inconsistent with its own
  Φ_b twelve lines later. The ISSI form is general and is what we implement.
- The thesis carries an obliquity factor and a Σ-renormalization that ISSI does
  not. Neither is a matter of notation, and neither survived: both are errors,
  and both are gone.
