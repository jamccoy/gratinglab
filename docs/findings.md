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

## The prototype's CXRO reader drops a row

**`CXRO_to_n_k` in `panter1.py` hardcodes `skip_header=3` against a file with
two header lines, so it silently discards the first data row.**

    [0] ' Au Density=19.32'
    [1] ' Energy(eV), Delta, Beta'
    [2] '  200.  0.0112813823  0.00951793324'     <- dropped
    [3] '  200.925018  0.0111386664  0.00957826339'

The lost row is the lowest energy, which is the **longest wavelength** — the
end of the range a grazing-incidence soft X-ray scan is most likely to want.
For the Au table it truncates 6.199 nm to 6.171 nm, a 0.029 nm shortfall, and
the `.ari` files in the corpus carry it because that code wrote them.

Nothing downstream ever noticed, for the same reason nothing noticed `t̂ = −d̂`:
no consumer looked at the endpoint. PCGrate was handed the table and
interpolated inside it.

**Consequence.** `materials.optical.read_cxro` detects the header instead of
counting it — a line is data when its first three whitespace-separated fields
all parse as floats. That is robust to the two- and three-header-line variants
both, and to an export that grows a line, and it cannot drop a row for being at
an end of the table, which is exactly where the range guard is most sensitive.

Pinned by `tests/test_materials.py::TestTheCxroReader::test_the_longest_wavelength_row_survives`,
with a companion measuring the size of the truncation so "a real amount, not
rounding" is a number rather than a claim.

**The general lesson**, and the reason this is written down: *porting faithfully
is not the same as porting correctly*. The plan for this milestone called the
reader "a port, not new code", which was right about the effort and wrong about
the care — a port inherits the bugs unless someone checks the input against the
assumption. The check took one command.

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
| Energy, $\Sigma \leq 1$ | ✗ up to 1.61 | ✓ 1.0000 exactly (Parseval) |
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
| 0.10 | 0.92902 |

The deviation scales with **phase excursion across the groove** — depth relative
to wavelength, hence working order — not with $\lambda/p$ or the number of
propagating orders. That is why $\Sigma$ looked stubbornly constant (~0.9065)
across wildly different geometries at fixed profile.

> **Both numbers were remeasured in M16-C**, after the symmetric flux obliquity
> landed. It is exactly 1 at $m=0$, where all the power sits in the shallow
> limit, so the top row is unchanged — which is the point of the check.

### Appendix D's two extra factors, and what replaced one of them

The thesis carries $\mathcal{E}_n \propto [\cos\beta_n/\cos\alpha]\|G_n\|^2$ and
then renormalises so $\sum_n \mathcal{E}_n = 1$. Both are **removed**, not made
optional — the first because it is asymmetric under $\alpha\leftrightarrow\beta$
and so breaks reciprocity, the second because it destroys the model's own error
signal.

> **The conclusion originally drawn here went one step too far.** This section
> used to continue: "for $N \to \infty$ grooves the efficiency is the
> norm-squared Fourier coefficient and nothing else". That is wrong by up to
> 64% away from Littrow, and the *symmetric* obliquity
> $4\cos\alpha\cos\beta_m/(\cos\alpha+\cos\beta_m)^2$ is what repairs it without
> touching reciprocity. See "The scalar solver was never checked against a
> theory it did not share", below. ISSI §2.1's groove function
> $f(x) = \exp\{ik[n(k)-1]y(x/p)\}$ carries no flux factor either, and is the
> one place this project's normative reference is the one in error.

The renormalisation has no such rehabilitation — forcing $\Sigma$ to unity would
erase exactly the
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

### The full sweep: 95.9%, and where the 4% is

796 mutants across the physics core, 763 killed:

| Module | Mutants | Survivors |
|---|---:|---:|
| `solvers/scalar.py` | 228 | 3 |
| `convergence.py` | 182 | 1 |
| `compare.py` | 136 | 7 |
| `geometry.py` | 128 | 2 (both equivalent) |
| `checks.py` | 96 | **18** |
| `solvers/base.py` | 18 | 2 |
| `result.py` | 8 | 0 |

**Every survivor in `checks.py` is in `check_reciprocity`** — 18 of the 33
total, from a function that is 12% of the mutated surface. Its own score is
81%, against 98–100% everywhere else.

That is worth stating plainly, because this project calls reciprocity *"the
sharpest available check"* and leans on it hard: it is the reason M9 could
choose reciprocity over energy conservation with confidence, and it is the
check that constrains the model rather than the arithmetic. **The check that
validates everything else is the least validated thing here.**

Three kinds of gap, from reading the diffs:

- **The order-selection strategy is unverified.** `np.argsort(strength)` →
  `np.argsort(None)` survives. The function deliberately tests the *strongest*
  orders; nothing asserts that it does, so the selection could be arbitrary and
  every test would still pass — reciprocity holds for whichever orders get
  picked.
- **Boundaries again**, the same shape as the `beta` finding: `>` → `>=` on
  both the order cap and the grazing margin.
- **The no-pairs branch returns a report nobody checks.** `max_violation=0.0`
  → `1.0` survives, so the early return for "nothing could be tested" is
  untested. A report claiming a measurement it never made is exactly what
  `Provenance` exists to prevent elsewhere.

**Closed.** All 18 are killed and `checks.py` is now at 96/96. What it took is
the part worth recording: the strategy could not be tested against the real
solver at all, because the selection happens inside `check_reciprocity` and
leaves no trace on the report. It needed a **recording stand-in** — a solver
reciprocal by construction (efficiency depends on the order index, never on
the incidence azimuth), whose per-order strength the test dictates, and which
remembers every illumination it was asked about. The orders it chose are then
recovered by inverting the reversed azimuth back through `sin_beta`.

That also reached a gap no scalar test could. The reversed illumination must
carry the **original polarization**; dropping it leaves the reverse solve on
`Illumination`'s default of TE, so the check compares two different physical
problems and blames the solver. Scalar neglects polarization, so that mutant
survives against it and always would — catching it needs a backend that
resolves polarization, which a stand-in can be made to do. Precisely the class
of bug the first contributed RCWA will meet.

Two mistakes made while writing these, both worth the warning:

- The first stand-in returned a fixed order window rather than the geometry's
  real propagating set. `check_reciprocity` re-derives that set independently
  and skips anything that disagrees, so every dictated order was silently
  dropped and `pairs_tested` came back 0.
- The first grazing test assumed a near-grazing order existed at a convenient
  geometry. None did — the margin is half a degree, so the geometry has to be
  *solved for*. `sin(beta_1) = sin(89.7°)` gives `alpha = -34.8489°`, which
  puts order +1 inside the margin with four others clear of it. Without that,
  "the skip happens" passes for free.

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

**But it is blind to normalisation**, and that blindness cost something — see
the next two findings. Any factor symmetric in $\alpha \leftrightarrow \beta_m$
passes reciprocity unchanged, so the check cannot tell a correct absolute scale
from a wrong one. And it was never run against a *coated* problem, which is how
a straightforward reciprocity violation sat in the absolute-efficiency path from
M15-D until M16-D.

---

## The scalar solver was never checked against a theory it did not share

Every analytic check in the suite until M16 — sawtooth `sinc²`, the binary
grating, the sinusoid's `J_m²` — descends from the same transmittance-function
picture the solver evaluates. They are strong tests of the quadrature and
structurally incapable of testing the model. Reciprocity constrains the phase
but not the scale (above). Nothing pinned the **absolute normalisation**.

First-order Rayleigh perturbation theory does. It comes from expanding the field
in plane waves and imposing the boundary condition on the corrugated surface —
no Kirchhoff assumption anywhere. Kirchhoff needs gentle slopes, perturbation
theory needs shallow grooves, and a grating that is both must be described
correctly by both. In that overlap, $\mathscr{E}_m = |G_m|^2$ was wrong:

| Mount | $\lambda/p$ | $m$ | $\beta_m$ | scalar / exact |
|---|---|---|---|---|
| in-plane, α=10° | 0.80 | +1 | 38.8° | 1.0137 |
| in-plane, α=10° | 0.80 | −1 | −76.8° | **1.6375** |
| off-plane, α=20°, γ=1.25° | 0.0127 | +1 | 13.9° | 1.0003 |
| off-plane, α=20°, γ=1.25° | 0.0127 | −1 | −67.5° | **1.2156** |

The discrepancy is exactly $(\cos\alpha+\cos\beta_m)^2/(4\cos\alpha\cos\beta_m)$
— to five significant figures, at every geometry tried. It vanishes at Littrow,
which is why forty tests and a thousand runs never saw it, and reaches 22% in
this project's own grazing-incidence regime.

The repair, `geometry.flux_obliquity`, is *symmetric*, so reciprocity survives;
$\leq 1$ by AM–GM, so the energy bound improves rather than degrades (worst-case
sum 1.71 → 1.61); and exactly 1 at $m=0$, so the shallow-limit energy identity
is untouched. Thesis Appendix D's $\cos\beta_m/\cos\alpha$ is the *asymmetric*
form and genuinely does break reciprocity: rejecting it was right, concluding
that no flux factor belonged was not.

Encoded as [`tests/test_perturbation.py`](../tests/test_perturbation.py). Two
traps found while writing it, both worth remembering:

- **`np.allclose` has `atol=1e-8`.** A shallow grating's first-order efficiency
  is ~1e-9, so the default tolerance passed the comparison against *anything*,
  including zero. The first draft of the file passed cleanly against the bug it
  was written to find. `atol=0` is load-bearing.
- **A near-Littrow mount tests nothing.** The first draft also used α=10° at
  λ/p=0.36 and α=25° at λ/p=0.01, where the factor is 0.2% from unity.
  `test_the_mounts_actually_discriminate` now asserts each mount reaches an
  order where the two theories are ≥15% apart.

---

## Naming a coating silently broke reciprocity

`check_reciprocity` was only ever pointed at bare problems. Pointed at a coated
one, on master before M16:

| | max violation |
|---|---|
| no coating, Blazed / Sinusoidal | 1 × 10⁻¹⁶ / 2 × 10⁻¹⁷ |
| `coating="Au"`, Blazed | 1.6 × 10⁻³ |
| `coating="Au"`, Sinusoidal | **3.5 × 10⁻²** |

The cause is plain once looked at: M15-D evaluated $R$ at a graze built from
$\alpha$ alone, so swapping $\alpha$ and $\beta_m$ changed the factor. Every
absolute-efficiency result since M15-D violated the invariant this solver is
built around, and the project's own sharpest check would have caught it on the
first run against a coated problem.

The fix is not to drop the reflectivity but to make it symmetric — weighting
each facet by both the angle it receives at and the angle it emits at,
$\sqrt{r(\zeta_i)\,r(\zeta_{d,m})}$. That restores $1.7\times10^{-16}$.
Pinned as `test_the_facet_model_breaks_reciprocity_which_is_why_it_is_not_default`,
because "the old model is kept for reproducibility" has to include reproducing
what was wrong with it.

---

## Reflectivity is not one number for the whole groove

`docs/theory/scalar.md` justified a single per-wavelength $R$ on the grounds
that "scalar theory has no mechanism by which one order could reflect
differently from its neighbour". There is one. A groove whose reflectivity
varies across the cycle is an **amplitude** grating as well as a phase grating,
and the two Fourier-transform together.

On the project's own reference geometry — Blazed 29.5°/70.5°, γ=1.25°,
α=19.99°, Au — the anti-blaze facet is **16.68% of the period** and its local
graze is negative: it faces away from the beam entirely. The M15 model applied
the active facet's $R(1.2328°)$ across 100% of the period regardless. Resolving
the cycle moves individual orders by **−51% to +12%**, and the movement is
order-dependent where the old factor was flat by construction.

Three details that had to be right:

- **The sign of the facet tilt.** $\tan\delta = +dy/dt$. The opposite sign gives
  0.8119° where the right one gives 1.2328°, and both look like angles. The
  check that settles it is that an ideal sawtooth must return $\delta$ equal to
  its own blaze angle, reproducing `facet_graze` exactly. See also
  "The profile parameter runs backwards", which is what makes the wrong sign
  tempting.
- **Two square roots, not one.** $\sqrt{r_i r_d}$ takes a principal branch of
  the *product* and adds a discontinuity of its own; on Au at 4 nm $\arg r_s$
  reaches −3.124, within 0.017 of the cut. $\sqrt{r_i}\sqrt{r_d}$ leaves only
  the discontinuities $r$ actually has.
- **Brewster is the one that remains.** Where $r_p$ passes through zero its
  phase jumps by π and the geometric mean cannot carry it. Reachable only with
  steep grooves near normal incidence; reported as a provenance warning rather
  than silently returned.

**What it costs.** The visibility mask puts a jump in the integrand at the
shadow boundary, so convergence drops from $O(n^{-2})$ to $O(n^{-1})$ and the
default 2048 points buys ~1e-4 rather than 1e-6 on a blazed profile. A smooth
profile is unaffected. Numbers in `theory/scalar.md` §8.

**Still unvalidated against a rigorous method** — and the corpus cannot fix
that; see below.

---

## The corpus can test the diffraction but not the reflectivity

Run against `OGRE/tastetest_perf_wavescan.txt` (TASTE, period 315.15 nm,
γ=1.25°, α=19.99°), comparing **relative** scalar against PCGrate:

| Profile | Σ, M16 | Σ, pre-M16 | per-order RMS, M16 / pre-M16 |
|---|---|---|---|
| Blazed 29.5°/70.5° | 0.5431 | 0.5514 | 0.11661 / 0.11671 |
| Blazed 29.5° ideal | 0.6301 | 0.6380 | 0.13634 / 0.13704 |
| `AFM_real_echelle.ggp` | 0.6138 | 0.6427 | 0.19639 / 0.19921 |

The flux obliquity improves agreement in every case, and only slightly. That is
consistent rather than disappointing: at this mount the blaze order leaves near
39° against a 19.99° incidence, so $O_m \approx 0.99$ for the orders carrying
the power, and §"the scalar solver was never checked" predicts a small change
exactly here. The residual RMS of ~0.12 is dominated by profile mismatch — which
`.ggp` this run used is still unconfirmed, and `corpus.toml` marks the entry
`confirmed = false`. It is not evidence about the normalisation either way.

**The corpus says nothing at all about M16-D.** Every usable reference run is
**perfect conductivity** — the TASTE table sums to 1.0005, i.e. R ≡ 1, so there
is no reflectivity in it to compare a reflectivity model against. The one
finite-conductivity run in the collection is `panter1_finite`, already marked
`usable = false` for summing to 3.6 at the Rayleigh anomalies.

So the groove-resolved reflectivity model rests on: reciprocity (which it
restores and the alternative breaks), an exact closed-form reduction on a
single-facet profile, and correct limiting behaviour. It has **no external
validation**, and acquiring some needs either a finite-conductivity PCGrate run
or the RCWA backend. Worth stating plainly, because the model changes individual
orders by up to 51% and is now the default.

---

## The measured groove is rounded, not faceted, and the solver reads the depth

The first thing the metrology merge (M17) made testable, and it went the
opposite way to the prediction.

`docs/roadmap.md` has long held that the TASTE residual — scalar's total near
0.55 against a PCGrate total of 1.0005 — is dominated by **profile mismatch**,
because the comparison runs on an idealised `Blazed(29.5°, 70.5°)` while the
`.ggp` the reference run actually used cannot be found. With the metrology
package absorbed, an AFM scan of what is plausibly that grating became
available to the same code (`TASTE_ALS_A205_Ti_Pt_flatten.txt` — A205, Ti/Pt),
so the idealised sawtooth can be replaced with a measured groove and the
hypothesis tested directly.

The scan is **not in the repository** — it is the group's measurement data, held
under `GRATINGLAB_AFM_DIR` (default `~/Documents/afm_scans/`) like the PCGrate
corpus. Everything below is reproducible from it; nothing below is reproducible
from a bare clone, which is the cost of not committing measurements and is
stated here rather than discovered.

Run at the corpus geometry (period 315.15 nm, γ=1.25°, α=19.99°, perfect
conductivity, 8192 quadrature points), mean total over the reference grid:

| Profile | mean Σ | vs PCGrate |
|---|---|---|
| PCGrate (integral method) | 1.0005 | — |
| idealised `Blazed(29.5°, 70.5°)` | 0.5420 | 0.542 |
| measured AFM groove, p = 315.15 nm | 0.3362 | 0.336 |
| measured AFM groove, p = 314.09 nm (its own) | 0.3354 | 0.335 |

**The measured profile is worse, by a lot.** Substituting real geometry for
idealised geometry moved the total *away* from unity. The period makes almost no
difference — 315.15 and 314.09 nm agree to 0.2% — so this is shape, not scale.

### The groove is inconsistent with itself

| Quantity | Value |
|---|---|
| Blaze angle from the **facet fit** (`extract_blaze_angle`, 100 grooves) | 27.91° ± 2.13° |
| Blaze angle implied by the measured **depth**, at a 70.5° anti-blaze | **20.33°** |
| depth/period, measured | 0.3275 |
| depth/period, sharp `Blazed(27.91°, 70.5°)` | 0.4460 |

On a *sharp* two-facet sawtooth these blaze angles are the same number: depth and
facet angle are locked together. Here they differ by **7.6°**. Something is
taking depth out of the groove while leaving the mid-facet slope intact.

> **Do not read that as 3.6σ.** An earlier version of this finding quoted the
> disagreement in units of the fit's own scatter, which sounds decisive and is
> not. The synthetic control (`tests/metrology/fixtures/synthetic_blazed_scan.txt`,
> an ideal 30° sawtooth) recovers 29.90° ± 0.09° fitted against 29.56° implied —
> a 0.33° gap that is pure discretisation, and **3.8σ**, because the scatter on
> ideal data is tiny. Sigma measures the noise, not the disagreement. The
> quantity that separates the two cases is the absolute gap: 0.33° on a sharp
> groove against 7.6° here, more than twenty times larger.

### Three mechanisms could do that. Two are ruled out.

**A flat land.** Real blazed gratings frequently have an unfaceted flat within
the period, and a land of fraction *f* gives
`depth = (1 − f)/(cot δ + cot δ′)` — a depth deficit with no metrology artefact
required. This is the conventional explanation and has to be excluded first.

Excluded. The observed depth needs **f ≈ 25–30%** of the period flat:

| Facet angles assumed | sharp depth | implied land |
|---|---|---|
| 27.91° / 70.5° | 0.4460 | 26.6% |
| 29.5° / 70.5° | 0.4713 | 30.5% |
| 27.91° / 67.45° | 0.4342 | 24.6% |

The profile does not have one. Only **2.9%** of the period lies within 2° of
flat, 7.2% within 5°, and 12.2% within 10° — well short of 25% at any reasonable
threshold. A land of that size would also put a sharp spike at 0° in the slope
distribution, and there is none: the histogram is broad and continuous, spreading
from −50° to 0° on the blaze side and 0° to +70° on the anti-blaze side, with no
mode at zero. A sharp faceted groove with a land would be *tri-modal*. This one
has no modes at all.

**Averaging across grooves.** The pipeline averages 5 grooves whose periods
differ by ±2.7% (314.33 ± 8.50 nm) and whose facet angles differ by ±2.13°, which
would blur an apex all by itself — a defect of the measurement pipeline rather
than of the grating or the tip.

Excluded. The individual grooves are 0.3304, 0.3274, 0.3268, 0.3289, 0.3312 —
mean 0.3289 ± 0.0017. Averaging costs **0.4%** of the depth. Every groove is
individually shallow; the average is not hiding sharp ones.

**Rounding.** What remains, and what the slope histogram positively supports: the
groove has no facets in the sense the fit assumes. The line fit lands on the
steepest part of a continuously curving flank and reports it as "the facet
angle", while the apex and trough are rounded off, removing depth.

### What rounding does *not* tell us

**Whether the grating is rounded or the tip rounded it is not determinable from
this scan.** A tip too blunt to reach the trough produces this signature; so does
a genuinely rounded groove from fabrication. Distinguishing them needs either a
tip characterisation, a scan of the same grating with a sharper tip, or a
cross-section. Nothing here settles it, and the earlier version of this finding
asserted tip convolution without excluding either alternative above.

### Why it matters more than it looks

The scalar phase term is `(2π/λ)·height·sinγ`. It depends on **depth**, not on
facet angle. So this groove diffracts like the shallower blaze its depth implies —
about 20° — regardless of the 27.9° its facets report. Whichever mechanism is
responsible, the facet fit is the number a person would quote for this grating
and the depth is the number the solver uses, and they disagree by 3.6σ.

### What this establishes

- It **falsifies**, for this profile, the claim that swapping in a measured
  groove would close the TASTE gap. It does the reverse.
- It does **not** show the profile-mismatch hypothesis is wrong in general — a
  rounded AFM boundary is not the true grating either, if the rounding is the
  tip's.
- It does **not** confirm this scan is the TASTE reference grating. The evidence
  is a name match, a period agreeing within uncertainty (314.33 ± 8.50 nm
  measured against 315.15 nm known independently from `OGRE/pcgratewavscan.py`),
  and a facet angle consistent with 29.5°. Suggestive, not decisive.
- Scalar theory does not conserve energy by construction (`theory/scalar.md` §5),
  so *part* of the gap to 1.0005 was never attributable to the profile at all.

**Consequence: an absolute efficiency from a raw AFM boundary is not currently
defensible**, and nothing in the pipeline says so. The cheapest available warning
is a diagnostic comparing the depth-implied blaze angle against the fitted one —
they agree on a sharp groove and diverge here, and the divergence is exactly the
quantity that matters. It must be thresholded in **degrees, not sigma**, for the
reason above. Recorded in `theory/metrology.md` §1.

The control now exists as a committed fixture, so the "they agree on a sharp
groove" half is a test rather than an assertion:
`test_depth_and_facet_fit_agree_on_a_sharp_groove`.

---

## Two ways to write a test that cannot fail

Both found while building `test_perturbation.py`, both applicable well beyond it.

**`np.allclose` carries `atol=1e-8`.** A shallow grating's first-order
efficiency is ~1e-9, three orders of magnitude below that floor, so the default
tolerance passed the comparison against *anything* — including a solver
returning zero. The first draft of the file passed cleanly against the very bug
it was written to find. Any comparison of small absolute quantities needs
`atol=0` and a relative tolerance.

**A test geometry can be degenerate without looking it.** The same first draft
used α=10° at λ/p=0.36 and α=25° at λ/p=0.013 — respectable-looking mounts where
the quantity under test sits 0.2% from unity, so agreement was guaranteed
whatever the solver did. The fix is an assertion about the *test setup* rather
than the result: `test_the_mounts_actually_discriminate` requires each mount to
reach an order where the two theories are at least 15% apart.

The pattern behind both: a passing test is evidence only if it could have
failed, and neither of these could. Non-vacuity assertions are cheap and the
suite already uses them elsewhere — these are two more places they were needed.

---

## Névot–Croce was returning an amplitude factor into an intensity

`fresnel.reflectivity` multiplied $|e^{-2k_{iz}k_{tz}\sigma^2}|$ into $|r|^2$.
The factor is the *amplitude* form; the intensity needs its square.

Nothing in the suite could see it. Bounded by 1, monotone in σ, exactly 1 at
σ=0, "Debye–Waller over-damps near θ_c" — every property tested was satisfied by
both forms. The limit that separates them is $n \to 1$, where the interface
disappears, $k_{tz} \to k_{iz}$, and the two models become the same expression:

| graze | NC (before) | DW | NC² |
|---|---|---|---|
| 5° | 0.963213 | 0.927771 | 0.927780 |
| 30° | 0.291214 | 0.084805 | 0.084806 |
| 85° | 0.007467 | 0.000056 | 0.000056 |

Roughened reflectivities were therefore too high — 0.732 where 0.536 was right,
at 15° graze with σ=0.5 nm at λ=2 nm, a 37% error in the factor.

It also falsifies a claim recorded in M15-E. "Debye–Waller over-damps by ~1e-2
near/above θ_c, with no reliable sign far below" was measuring this bug. With
the intensity form the picture is coherent and the sign is reliable in both
regimes: DW over-damps near and above θ_c, under-damps deep below it by ~2e-3,
and the two converge to 1e-5 *far above* θ_c — not below, which is where the old
test looked.

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
