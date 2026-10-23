# Findings

This page records empirical results established while building this project,
each with the evidence that established it. They are recorded so that they are
not rediscovered, and so that a claim can be checked against its evidence.

Findings that constrain how code must behave are also encoded as tests, and
those tests are cross-referenced.

---

## Quadrature error is not monotone in the number of points

Refining a blazed-grating quadrature can make the answer worse, and a
convergence rule that trusts a single agreement will certify a tolerance that
the next refinement violates.

The measurement used the scalar solver at the reference geometry (period
315.15 nm, δ = 29.5°, anti 70.5°, α = 25°, γ = 1.5°) over λ ∈ [1, 5] nm. The
table gives the largest change in efficiency per doubling of
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

Past ~8192 the ratio settles near 4, i.e. second-order convergence, which is the
expected rate for a rule integrating a function with a slope discontinuity.
Below that it swings between 0.4 and 23.

The cause is geometric and is not numerical noise. The blazed profile has a
kink at `t = 1/(1 + tan δ / tan δ') = 0.8331`. That value is not a dyadic
rational, so doubling `n` changes how near the closest node falls to the kink,
and the local error at the kink dominates. A sinusoid, being analytic, shows
none of this behaviour: it is at machine precision (~1.7e-16) from `n = 64`
upward and never moves again. A lamellar profile, which has a jump in place of
a kink, converges cleanly at exactly 4× per doubling with no excursions.

This consequence is encoded in code. `gratinglab.convergence` requires a
*plateau*, defined as `DEFAULT_PLATEAU = 2` consecutive differences below
tolerance, so that three knob values in a row must agree. At a tolerance of
1e-5 the naive one-agreement rule stops at n = 512 on the strength of the
9.6e-6 difference, and the next refinement then moves the answer by 2.4e-5. One
extra solve rejects that stopping point.

For the same reason the harness reports the whole ladder in place of a boolean,
and `converged_at` is the coarsest value of the plateau, not the first value
that happened to agree.

The behaviour is pinned by `tests/test_convergence.py::TestTheMeasurementBehindTheRule`,
which fails if refinement ever becomes monotone here (at which point this entry
and the plateau argument would both need revisiting), and by
`TestThePlateauRule::test_and_plateau_1_would_have_stopped_at_the_dip`, which
shows the naive rule failing on a scripted sequence with the same shape.

---

## The prototype's CXRO reader drops a row

`CXRO_to_n_k` in `panter1.py` hardcodes `skip_header=3` against a file with two
header lines, so it silently discards the first data row.

    [0] ' Au Density=19.32'
    [1] ' Energy(eV), Delta, Beta'
    [2] '  200.  0.0112813823  0.00951793324'     <- dropped
    [3] '  200.925018  0.0111386664  0.00957826339'

The lost row is the lowest energy, which is the longest wavelength. That is the
end of the range a grazing-incidence soft X-ray scan is most likely to need.
For the Au table the defect truncates 6.199 nm to 6.171 nm, a 0.029 nm
shortfall, and the `.ari` files in the corpus carry it because that code wrote
them.

Nothing downstream detected it, for the same reason nothing detected
`t̂ = −d̂`: no consumer examined the endpoint. PCGrate was handed the table and
interpolated inside it.

As a consequence, `materials.optical.read_cxro` detects the header instead of
counting it: a line is treated as data when its first three
whitespace-separated fields all parse as floats. This is robust to both the
two- and three-header-line variants and to an export that grows a line, and it
cannot drop a row for being at an end of the table, which is where the range
guard is most sensitive.

The behaviour is pinned by `tests/test_materials.py::TestTheCxroReader::test_the_longest_wavelength_row_survives`,
with a companion test that measures the size of the truncation, so that the
statement that it is a real amount and not rounding rests on a number.

The general lesson, and the reason this entry was written, is that porting
faithfully is not the same as porting correctly. The plan for this milestone
called the reader "a port, not new code", which was accurate about the effort
but not about the care required: a port inherits the bugs of its source unless
someone checks the input against the assumption. In this case the check took
one command.

---

## Corpus geometry recovered from the data alone

The exported efficiency tables record efficiencies but not the geometry that
produced them, and the vendor project files are binary. Running `strings` on
them yields only a build stamp and a path to an optical-constants file.

The geometry is nonetheless exactly recoverable. Which orders propagate is
fixed by the grating equation, and the constraint factors:

$$m_{\max} = \left\lfloor A/\lambda \right\rfloor, \quad
  m_{\min} = -\left\lfloor B/\lambda \right\rfloor$$

$$A = (1 + \sin\alpha)\,p\sin\gamma \qquad B = (1 - \sin\alpha)\,p\sin\gamma$$

$A$ and $B$ decouple, so each is obtained by a simple interval intersection over
the scan. Across a 500-point wavelength scan this determines both to about
±0.01 nm.

For validation, the period of the TASTE run is independently known: it is
315.15 nm, hard-coded in the group's own plotting script. Using that period
gives

$$\gamma = 1.2500°, \qquad \alpha = +19.99°$$

Both are clean round numbers, which a mis-derived method would not be expected
to produce.

A test asserts that this project's grating equation, driven by the recovered
geometry, regenerates the exact propagating-order set of the reference file at
every sampled wavelength, and the test passes.

The results are stored in [`benchmarks/corpus.toml`](../benchmarks/corpus.toml).
The panter1 period remains unresolved; see [roadmap](roadmap.md#the-panter1-grating-period).

> Sign trap. The first pass assigned $A$ and $B$ the other way round and
> obtained α = −19.99°, which looked plausible. The error was caught by running
> the scalar solver at both signs: the mean summed efficiency was 0.23 at
> −19.99° against 0.55 at +19.99°. A sign error in a recovered quantity gives no
> indication of itself.

---

## A reference file violates energy conservation

`panter1_fix_efftable_finite.txt`, a finite-conductivity run, has summed
efficiency reaching 3.6. This was verified against the raw file, so it is a
defect in the run and not in the parser. The companion perfect-conductivity run
of the same problem conserves energy to four decimal places, so the setup is
sound, and it is the finite-conductivity solve that destabilises.

79 of 561 points (14 %) violate conservation, in five clusters:

| λ range (nm) | points | orders | peak Σℰ |
|---|---|---|---|
| 0.60–0.89 | 30 | 23–33 | 2.94 |
| 0.97–1.03 | 7 | 19–20 | 1.23 |
| 1.11–1.12 | 2 | 18 | 3.60 |
| 1.44–1.49 | 6 | 14 | 2.03 |
| 4.16–4.49 | 34 | 5 | 3.56 |

The number of propagating orders does not explain the violations, since the
last cluster carries only five orders. What predicts them is order passing-off
(Rayleigh anomalies): violating points lie a median 0.020 nm from a wavelength
at which the propagating-order count changes, against 0.125 nm for clean
points, and the last cluster terminates exactly at the 5→4 passing-off
wavelength. The gold optical constants are smooth across all five bands, so no
material edge is involved.

This file should not be used as reference data at those wavelengths. The defect
is encoded as `test_known_defect_finite_conductivity_run_is_unphysical`; if the
run is ever regenerated correctly, that test will fail and should be deleted.

---

## Scalar theory cannot be reciprocal and conserve energy

No scalar formulation satisfies all three of the properties below. The
limitation is structural: reciprocity requires the phase to be symmetric under
$\alpha \leftrightarrow \beta_m$, while Parseval's theorem requires a single
fixed function, which cannot contain $\beta_m$.

| Property | order-dependent $\varphi = kg\sin\gamma(\cos\alpha + \cos\beta_m)$ | $\beta$-free $\varphi = 2kg\sin\gamma\cos\alpha$ |
|---|---|---|
| Energy, $\Sigma \leq 1$ | ✗ up to 1.61 | ✓ 1.0000 exactly (Parseval) |
| Reciprocity $E_m(\alpha) = E_m(\beta_m)$ | ✓ $1\times10^{-17}$ | ✗ 0.44 violation |
| Blaze direction | ✓ exact at all $\alpha$ | ✓ at Littrow only; 7.6° off at 10°, evanescent past ~15° |

Physical optics is reciprocal but does not conserve energy, and the pure
transmittance picture has the opposite properties. The solver keeps
reciprocity. The energy-conserving variant is not offered as an option, because
it would silently exchange a well-understood energy defect for a reciprocity
violation, which a user could easily select without realising the cost.

Referencing the phase to the blaze direction $\beta_b = 2\delta - \alpha$ appears
to offer a third option but does not. $\beta_b$ is an output of the model (it is
where the $\operatorname{sinc}^2$ envelope peaks), so feeding it back in is
circular, and the "zero error at every $\alpha$" it produces is a tautology.

### Check of the implementation in the shallow-groove limit

The implementation is not at fault. In the shallow-groove limit the sum is
exactly unity, and it degrades smoothly with depth:

| depth / period | $\Sigma_m E_m$ |
|---|---|
| 0.0001 | 1.00000 |
| 0.10 | 0.92902 |

The deviation scales with the phase excursion across the groove, that is, with
depth relative to wavelength and hence with working order. It does not scale
with $\lambda/p$ or with the number of propagating orders. This is why $\Sigma$
remained nearly constant (~0.9065) across very different geometries at fixed
profile.

> Both numbers were remeasured in M16-C, after the symmetric flux obliquity
> landed. That factor is exactly 1 at $m=0$, where all the power sits in the
> shallow limit, so the top row is unchanged, as the check requires.

### Appendix D's two extra factors, and what replaced one of them

The thesis carries $\mathcal{E}_n \propto [\cos\beta_n/\cos\alpha]\|G_n\|^2$ and
then renormalises so that $\sum_n \mathcal{E}_n = 1$. Both factors are removed,
and neither is retained as an option. The first is asymmetric under
$\alpha\leftrightarrow\beta$ and so breaks reciprocity; the second destroys the
model's own error signal.

> The conclusion originally drawn here went one step too far. This section
> used to continue: "for $N \to \infty$ grooves the efficiency is the
> norm-squared Fourier coefficient and nothing else". That statement is wrong
> by up to 64% away from Littrow, and the symmetric obliquity
> $4\cos\alpha\cos\beta_m/(\cos\alpha+\cos\beta_m)^2$ repairs it without
> affecting reciprocity. See "The scalar solver was never checked against a
> theory it did not share", below. ISSI §2.1's groove function
> $f(x) = \exp\{ik[n(k)-1]y(x/p)\}$ carries no flux factor either, and it is the
> one place where this project's normative reference is the one in error.

The renormalisation has no comparable rehabilitation: forcing $\Sigma$ to unity
would erase the validity signal described above.

The full derivation is in [theory/scalar.md §5](theory/scalar.md). The
deviation is reported as a provenance warning and is never rescaled.

---

## Mutation testing found a real gap

Mutation testing here means deliberately corrupting the physics and checking
whether the suite detects the change. The first sweep applied seven mutations,
of which six were caught and one survived.

### The obliquity factor was unverified

Inverting $\cos\beta_m/\cos\alpha \to \cos\alpha/\cos\beta_m$ passed the entire
suite. The only test asserted that the result changed when the option was
enabled; none asserted that it changed correctly.

This is the thesis-versus-ISSI discrepancy recorded in `conventions.md` §5,
where an error was plausible and not merely hypothetical. The factor was pinned
to the exact ratio, with a second test showing that the two ratios are
distinguishable for the geometry used, so that the assertion could not become
vacuous.

> Since superseded. The obliquity option was removed entirely once the factor
> was established to be an error and not a convention; see "Scalar theory
> cannot be reciprocal and conserve energy" above. The lesson still applies,
> and it is the reason that section could be written with confidence: a test
> that only asserts that something changed verifies nothing.

### In-plane tests are blind to γ errors

Replacing $\sin\gamma \to \sin^2\gamma$ left 34 of 37 scalar tests passing.
In-plane cases have $\sin\gamma = 1$, which makes every power of it identical,
so they cannot detect an error in the γ handling at all. Only the off-plane
cases failed, and the off-plane mount is the primary application.

Coverage at the time was 9 off-plane tests against 13 in-plane. Every
closed-form check now runs in both mounts.

After both fixes, all 12 mutations are caught.

### A sweep with an automated mutation tool found two more

The hand-written sweep consisted of seven mutations chosen by someone who knew
where to look. `mutmut` generated 19 for `geometry.beta` alone, and four
survived. Two of them were real gaps:

- **`np.pi - b` → `np.pi + b` for transmitted orders.** `beta(..., transmitted=True)`
  had no test at all; only `cos_beta` did. The mutant negates \(\sin\beta\)
  while leaving \(\cos\beta\) negative, which is the case `conventions.md` §4
  forbids ("we do not flip the sign of \(\sin\beta_m\) for transmission"). It is
  now pinned by one test that asserts both halves of that sentence.
- **`|s| <= 1.0` → `|s| < 1.0`.** At exactly \(|\sin\beta_m| = 1\), the
  passing-off boundary, the mutant returns NaN. `is_propagating` uses `<=`, so
  an order it counts as propagating would have had no direction at which to be
  drawn or summed. The two functions are now required to agree at the
  boundary, with a companion test showing that the boundary is a real one.

The other two are *equivalent mutants*, that is, changes with no observable
effect; this was verified and not assumed:

- `<= 1.0` → `<= 2.0` is unobservable, because `arcsin` of any value outside
  \([-1, 1]\) is NaN in any case and `errstate(invalid="ignore")` suppresses the
  warning either way. The guard therefore does less than it appears to.
- Dropping `dtype=np.float64` from `np.asarray` changes nothing that NumPy does
  not already do by type promotion.

No test can kill either mutant, and attempting to would mean writing assertions
about nothing. They are recorded here so that the next reader does not try.

This sweep extends the conclusion of the earlier one, that a test which only
asserts that something changed verifies nothing. It adds that a branch with no
test at all is invisible to code coverage when another function exercises the
same line. `cos_beta(0.5, transmitted=True)` covered the concept, but nothing
covered the version of it in `beta`.

### The full sweep: 95.9%, and where the 4% is

The sweep produced 796 mutants across the physics core, of which 763 were
killed:

| Module | Mutants | Survivors |
|---|---:|---:|
| `solvers/scalar.py` | 228 | 3 |
| `convergence.py` | 182 | 1 |
| `compare.py` | 136 | 7 |
| `geometry.py` | 128 | 2 (both equivalent) |
| `checks.py` | 96 | **18** |
| `solvers/base.py` | 18 | 2 |
| `result.py` | 8 | 0 |

Every survivor in `checks.py` is in `check_reciprocity`. These account for 18 of
the 33 survivors in total, from a function that makes up 12% of the mutated
surface. Its own score is 81%, against 98–100% everywhere else.

This matters because the project calls reciprocity "the sharpest available
check" and relies on it heavily. It is the reason M9 could choose reciprocity
over energy conservation with confidence, and it is the check that constrains
the model as opposed to the arithmetic. The results show that the check which
validates everything else was the least validated part of the code.

Reading the diffs showed three kinds of gap:

- **The order-selection strategy is unverified.** `np.argsort(strength)` →
  `np.argsort(None)` survives. The function deliberately tests the strongest
  orders, but nothing asserts that it does, so the selection could be arbitrary
  and every test would still pass, since reciprocity holds for whichever orders
  are picked.
- **Boundaries again**, with the same shape as the `beta` finding: `>` → `>=`
  on both the order cap and the grazing margin.
- **The no-pairs branch returns a report that nothing checks.**
  `max_violation=0.0` → `1.0` survives, so the early return for the case where
  nothing could be tested is itself untested. A report claiming a measurement it
  never made is what `Provenance` exists to prevent elsewhere.

These gaps are now closed. All 18 mutants are killed and `checks.py` is at
96/96. The strategy could not be tested against the real solver at all, because
the selection happens inside `check_reciprocity` and leaves no trace on the
report. Testing it required a *recording stand-in*: a solver that is reciprocal
by construction (efficiency depends on the order index and never on the
incidence azimuth), whose per-order strength is dictated by the test, and which
records every illumination it is asked about. The orders it chose are then
recovered by inverting the reversed azimuth back through `sin_beta`.

The stand-in also reached a gap that no scalar test could. The reversed
illumination must carry the original polarization. If it is dropped, the
reverse solve falls back to the default TE of `Illumination`, so the check
compares two different physical problems and attributes the difference to the
solver. Scalar theory neglects polarization, so that mutant survives against
the scalar solver and always would; catching it requires a backend that
resolves polarization, which a stand-in can be made to do. This is the class of
bug the first contributed RCWA backend is expected to encounter.

Two mistakes were made while writing these tests, and both are recorded as
warnings:

- The first stand-in returned a fixed order window in place of the geometry's
  real propagating set. `check_reciprocity` re-derives that set independently
  and skips any order that disagrees, so every dictated order was silently
  dropped and `pairs_tested` came back 0.
- The first grazing test assumed that a near-grazing order existed at a
  convenient geometry. None did: the margin is half a degree, so the geometry
  has to be solved for. `sin(beta_1) = sin(89.7°)` gives `alpha = -34.8489°`,
  which places order +1 inside the margin with four others clear of it. Without
  that geometry, the test that the skip happens would pass without testing
  anything.

---

## Reciprocity is the sharpest available check

Lorentz reciprocity requires $\mathscr{E}_m(\alpha) = \mathscr{E}_m(\beta_m)$.
Scalar theory satisfies it exactly, because the phase
$\cos\alpha + \cos\beta_m$ is symmetric under the exchange.

| | max violation |
|---|---|
| as implemented | 5 × 10⁻¹⁸ |
| with the α↔β symmetry of Φ broken | 4 × 10⁻¹ |

The two cases differ by sixteen orders of magnitude, and the check requires
neither reference data nor a closed form.

The closed-form tests cannot provide this. They compare the solver against a
formula derived in the same way as the solver computes it, so they validate the
quadrature. Reciprocity constrains the structure of the phase function.

Reciprocity is, however, blind to normalisation, and that limitation had a
cost, described in the next two findings. Any factor symmetric in
$\alpha \leftrightarrow \beta_m$ passes reciprocity unchanged, so the check
cannot distinguish a correct absolute scale from a wrong one. It was also never
run against a coated problem, which is how a straightforward reciprocity
violation remained in the absolute-efficiency path from M15-D until M16-D.

---

## The scalar solver was never checked against a theory it did not share

Every analytic check in the suite until M16 (the sawtooth `sinc²`, the binary
grating, and the sinusoid's `J_m²`) descends from the same transmittance-function
picture that the solver evaluates. They are strong tests of the quadrature and
are structurally incapable of testing the model. Reciprocity constrains the
phase but not the scale, as described above. Nothing constrained the absolute
normalisation.

First-order Rayleigh perturbation theory does constrain it. It is obtained by
expanding the field in plane waves and imposing the boundary condition on the
corrugated surface, with no Kirchhoff assumption at any point. Kirchhoff theory
requires gentle slopes and perturbation theory requires shallow grooves, so a
grating satisfying both conditions must be described correctly by both. In
that overlap, $\mathscr{E}_m = |G_m|^2$ was found to be wrong:

| Mount | $\lambda/p$ | $m$ | $\beta_m$ | scalar / exact |
|---|---|---|---|---|
| in-plane, α=10° | 0.80 | +1 | 38.8° | 1.0137 |
| in-plane, α=10° | 0.80 | −1 | −76.8° | **1.6375** |
| off-plane, α=20°, γ=1.25° | 0.0127 | +1 | 13.9° | 1.0003 |
| off-plane, α=20°, γ=1.25° | 0.0127 | −1 | −67.5° | **1.2156** |

The discrepancy is exactly $(\cos\alpha+\cos\beta_m)^2/(4\cos\alpha\cos\beta_m)$,
to five significant figures, at every geometry tried. It vanishes at Littrow,
which explains why forty tests and a thousand runs did not reveal it, and it
reaches 22% in this project's own grazing-incidence regime.

The repair, `geometry.flux_obliquity`, has three properties. It is symmetric, so
reciprocity survives. It is $\leq 1$ by the AM–GM inequality, so the energy
bound improves (worst-case sum 1.71 → 1.61). It is exactly 1 at $m=0$, so the
shallow-limit energy identity is untouched. The $\cos\beta_m/\cos\alpha$ of
thesis Appendix D is the asymmetric form and does break reciprocity: rejecting
it was correct, but concluding that no flux factor belonged was not.

The result is encoded as [`tests/test_perturbation.py`](../tests/test_perturbation.py).
Two traps were found while writing it:

- **`np.allclose` has `atol=1e-8`.** A shallow grating's first-order efficiency
  is ~1e-9, so the default tolerance passed the comparison against any value,
  including zero. The first draft of the file passed cleanly against the bug it
  was written to find. The setting `atol=0` is therefore essential.
- **A near-Littrow mount tests nothing.** The first draft also used α=10° at
  λ/p=0.36 and α=25° at λ/p=0.01, where the factor is 0.2% from unity.
  `test_the_mounts_actually_discriminate` now asserts that each mount reaches an
  order where the two theories are ≥15% apart.

---

## A coating broke reciprocity in the absolute-efficiency path

`check_reciprocity` had only ever been applied to bare problems. Applied to a
coated one, on master before M16, it gave:

| | max violation |
|---|---|
| no coating, Blazed / Sinusoidal | 1 × 10⁻¹⁶ / 2 × 10⁻¹⁷ |
| `coating="Au"`, Blazed | 1.6 × 10⁻³ |
| `coating="Au"`, Sinusoidal | **3.5 × 10⁻²** |

The cause is clear on inspection. M15-D evaluated $R$ at a graze constructed
from $\alpha$ alone, so swapping $\alpha$ and $\beta_m$ changed the factor.
Every absolute-efficiency result since M15-D therefore violated the invariant
around which this solver is built, and the project's own sharpest check would
have detected it on the first run against a coated problem.

The fix keeps the reflectivity and makes it symmetric, weighting each facet by
both the angle at which it receives and the angle at which it emits,
$\sqrt{r(\zeta_i)\,r(\zeta_{d,m})}$. This restores $1.7\times10^{-16}$. The old
behaviour is pinned as `test_the_facet_model_breaks_reciprocity_which_is_why_it_is_not_default`,
because keeping the old model for reproducibility includes reproducing what was
wrong with it.

---

## Reflectivity varies across the groove

`docs/theory/scalar.md` justified a single per-wavelength $R$ on the grounds
that "scalar theory has no mechanism by which one order could reflect
differently from its neighbour". Such a mechanism exists. A groove whose
reflectivity varies across the cycle is an amplitude grating as well as a phase
grating, and the two are Fourier-transformed together.

On the project's own reference geometry (Blazed 29.5°/70.5°, γ=1.25°,
α=19.99°, Au), the anti-blaze facet occupies 16.68% of the period and its local
graze is negative, meaning that it faces entirely away from the beam. The M15
model nevertheless applied the active facet's $R(1.2328°)$ across 100% of the
period. Resolving the cycle moves individual orders by −51% to +12%, and the
change is order-dependent, whereas the old factor was flat by construction.

Three details had to be correct:

- **The sign of the facet tilt.** $\tan\delta = +dy/dt$. The opposite sign gives
  0.8119° where the correct one gives 1.2328°, and both values are plausible
  angles. The deciding check is that an ideal sawtooth must return $\delta$
  equal to its own blaze angle, reproducing `facet_graze` exactly. See also
  "The profile parameter runs backwards", which explains why the wrong sign is
  easy to adopt.
- **Two square roots in place of one.** $\sqrt{r_i r_d}$ takes a principal branch
  of the product and introduces a discontinuity of its own; on Au at 4 nm,
  $\arg r_s$ reaches −3.124, within 0.017 of the branch cut.
  $\sqrt{r_i}\sqrt{r_d}$ leaves only the discontinuities that $r$ actually has.
- **Brewster's angle is the one discontinuity that remains.** Where $r_p$ passes
  through zero its phase jumps by π, and the geometric mean cannot represent the
  jump. This case is reachable only with steep grooves near normal incidence,
  and it is reported as a provenance warning instead of being returned
  silently.

The model has a cost. The visibility mask puts a jump in the integrand at the
shadow boundary, so convergence drops from $O(n^{-2})$ to $O(n^{-1})$, and the
default 2048 points achieve ~1e-4 instead of 1e-6 on a blazed profile. A smooth
profile is unaffected. The numbers are in `theory/scalar.md` §8.

The model has not yet been validated against a rigorous method, and the corpus
cannot provide that validation, as the next section explains.

---

## The corpus can test the diffraction but not the reflectivity

The comparison was run against `OGRE/tastetest_perf_wavescan.txt` (TASTE,
period 315.15 nm, γ=1.25°, α=19.99°), comparing relative scalar efficiencies
against a perfect-conductivity reference table from an independent
integral-method code. Comparisons against the commercial code are not published
here, because its licence terms have not been checked, so those figures remain
with the private corpus. This project's own totals are:

| Profile | Σ, M16 | Σ, pre-M16 |
|---|---|---|
| Blazed 29.5°/70.5° | 0.5431 | 0.5514 |
| Blazed 29.5° ideal | 0.6301 | 0.6380 |
| `AFM_real_echelle.ggp` | 0.6138 | 0.6427 |

The flux obliquity improves agreement in every case, but only slightly. This is
consistent with expectation: at this mount the blaze order leaves near 39°
against a 19.99° incidence, so $O_m \approx 0.99$ for the orders carrying the
power, and §"the scalar solver was never checked" predicts a small change in
exactly this case. The residual per-order disagreement is dominated by profile
mismatch. Which `.ggp` this run used is still unconfirmed, and `corpus.toml`
marks the entry `confirmed = false`. The comparison is therefore not evidence
about the normalisation in either direction.

The corpus provides no evidence about M16-D. Every usable reference run is for
perfect conductivity: the TASTE table sums to unity, i.e. R ≡ 1, so it contains
no reflectivity against which a reflectivity model could be compared. The one
finite-conductivity run in the collection is `panter1_finite`, already marked
`usable = false` for summing to 3.6 at the Rayleigh anomalies.

The groove-resolved reflectivity model therefore rests on three grounds:
reciprocity (which it restores and the alternative breaks), an exact
closed-form reduction on a single-facet profile, and correct limiting
behaviour. It has no external validation, and obtaining some requires either a
finite-conductivity reference run or the RCWA backend. This limitation is
significant because the model changes individual orders by up to 51% and is
now the default.

---

## The measured groove is rounded rather than faceted, and the solver uses its depth

This was the first result the metrology merge (M17) made testable, and it went
in the opposite direction to the prediction.

`docs/roadmap.md` has long held that the TASTE residual (a scalar total near
0.55 against the unit total of the perfect-conductivity reference) is dominated
by profile mismatch, because the comparison uses an idealised
`Blazed(29.5°, 70.5°)` while the `.ggp` the reference run actually used cannot
be found. With the metrology package absorbed, an AFM scan of what is plausibly
that grating became available to the same code
(`TASTE_ALS_A205_Ti_Pt_flatten.txt`, A205, Ti/Pt). The idealised sawtooth could
therefore be replaced with a measured groove and the hypothesis tested
directly.

The scan is not in the repository. It is the group's measurement data, held
under `GRATINGLAB_AFM_DIR` (default `~/Documents/afm_scans/`) in the same way as
the reference corpus. Everything below is reproducible from the scan, and
nothing below is reproducible from a bare clone; that is the cost of not
committing measurements.

The run used the corpus geometry (period 315.15 nm, γ=1.25°, α=19.99°, perfect
conductivity, 8192 quadrature points). The table gives the mean total over the
reference grid, on which a perfect-conductivity reference sums to unity:

| Profile | mean Σ |
|---|---|
| idealised `Blazed(29.5°, 70.5°)` | 0.5420 |
| measured AFM groove, p = 315.15 nm | 0.3362 |
| measured AFM groove, p = 314.09 nm (its own) | 0.3354 |

The measured profile performs considerably worse. Substituting real geometry for
idealised geometry moved the total away from unity. The period makes almost no
difference (315.15 and 314.09 nm agree to 0.2%), so the effect is one of shape
and not of scale.

### Depth and facet angle disagree

| Quantity | Value |
|---|---|
| Blaze angle from the **facet fit** (`extract_blaze_angle`, 100 grooves) | 27.91° ± 2.13° |
| Blaze angle implied by the measured **depth**, at a 70.5° anti-blaze | **20.33°** |
| depth/period, measured | 0.3275 |
| depth/period, sharp `Blazed(27.91°, 70.5°)` | 0.4460 |

On a sharp two-facet sawtooth these two blaze angles are the same number,
because depth and facet angle are locked together. Here they differ by 7.6°.
Some feature of the groove removes depth while leaving the mid-facet slope
intact.

> This disagreement should not be read as 3.6σ. An earlier version of this
> finding quoted it in units of the fit's own scatter, which appears decisive
> but is not. The synthetic control
> (`tests/metrology/fixtures/synthetic_blazed_scan.txt`, an ideal 30° sawtooth)
> recovers 29.90° ± 0.09° fitted against 29.56° implied. That 0.33° gap is
> entirely discretisation, yet it is 3.8σ, because the scatter on ideal data is
> very small. Sigma measures the noise and not the disagreement. The quantity
> that separates the two cases is the absolute gap: 0.33° on a sharp groove
> against 7.6° here, more than twenty times larger.

### Three candidate mechanisms, two of them excluded

The first candidate is a flat land. Real blazed gratings frequently have an
unfaceted flat within the period, and a land of fraction *f* gives
`depth = (1 − f)/(cot δ + cot δ′)`, a depth deficit that requires no metrology
artefact. This is the conventional explanation and has to be excluded first.

It is excluded. The observed depth would require f ≈ 25–30% of the period to be
flat:

| Facet angles assumed | sharp depth | implied land |
|---|---|---|
| 27.91° / 70.5° | 0.4460 | 26.6% |
| 29.5° / 70.5° | 0.4713 | 30.5% |
| 27.91° / 67.45° | 0.4342 | 24.6% |

The profile has no such land. Only 2.9% of the period lies within 2° of flat,
7.2% within 5°, and 12.2% within 10°, which is well short of 25% at any
reasonable threshold. A land of that size would also produce a sharp spike at
0° in the slope distribution, and there is none. The histogram is broad and
continuous, spreading from −50° to 0° on the blaze side and from 0° to +70° on
the anti-blaze side, with no mode at zero. A sharp faceted groove with a land
would be *tri-modal*, whereas this histogram has no modes at all.

The second candidate is averaging across grooves. The pipeline averages 5
grooves whose periods differ by ±2.7% (314.33 ± 8.50 nm) and whose facet angles
differ by ±2.13°, which could by itself blur an apex. That would be a defect of
the measurement pipeline, as opposed to one of the grating or the tip.

It is also excluded. The individual grooves have depth/period values of 0.3304,
0.3274, 0.3268, 0.3289 and 0.3312, with a mean of 0.3289 ± 0.0017. Averaging
costs 0.4% of the depth. Every groove is individually shallow, so the average
is not concealing sharp ones.

The remaining candidate is rounding, which the slope histogram positively
supports. In this explanation the groove has no facets in the sense the fit
assumes. The line fit lands on the steepest part of a continuously curving
flank and reports it as the facet angle, while the apex and trough are rounded
off, which removes depth.

### What the rounding does not establish

This scan cannot determine whether the grating itself is rounded or the tip
rounded it. A tip too blunt to reach the trough produces this signature, and so
does a groove that was genuinely rounded during fabrication. Distinguishing
them requires a tip characterisation, a scan of the same grating with a sharper
tip, or a cross-section. Nothing here settles the question, and the earlier
version of this finding asserted tip convolution without excluding either
alternative above.

*Narrowed since (M22): a tip in nominal condition is now excluded. Erosion with
the specified probe recovers no depth from this scan, and reproducing the
deficit in the forward direction requires an apex worn to ~80 nm. See "A
nominal tip does not explain the rounded groove" below.*

### Consequence for the solver

The scalar phase term is `(2π/λ)·height·sinγ`, which depends on depth and not on
facet angle. This groove therefore diffracts like the shallower blaze its depth
implies, about 20°, irrespective of the 27.9° its facets report. Whatever the
mechanism, the facet fit gives the number a person would quote for this
grating, the depth gives the number the solver uses, and the two disagree by
3.6σ.

### What this establishes

- For this profile, it falsifies the claim that substituting a measured groove
  would close the TASTE gap. The substitution has the opposite effect.
- It does not show that the profile-mismatch hypothesis is wrong in general,
  because a rounded AFM boundary is not the true grating either if the rounding
  comes from the tip.
- It does not confirm that this scan is of the TASTE reference grating. The
  evidence is a name match, a period agreeing within uncertainty
  (314.33 ± 8.50 nm measured against 315.15 nm known independently from
  `OGRE/pcgratewavscan.py`), and a facet angle consistent with 29.5°. This
  evidence is suggestive but not decisive.
- Scalar theory does not conserve energy by construction (`theory/scalar.md` §5),
  so part of the gap to unity was never attributable to the profile at all.

It follows that an absolute efficiency computed from a raw AFM boundary is not
currently defensible, and nothing in the pipeline reports this. The simplest
available warning is a diagnostic comparing the depth-implied blaze angle with
the fitted one. The two agree on a sharp groove and diverge here, and the
divergence is the quantity that matters. For the reason given above, the
diagnostic must be thresholded in degrees and not in sigma. This is recorded in
`theory/metrology.md` §1.

The control now exists as a committed fixture, so the statement that the two
angles agree on a sharp groove is covered by a test:
`test_depth_and_facet_fit_agree_on_a_sharp_groove`.

---

## Two forms of test that cannot fail

Both were found while building `test_perturbation.py`, and both apply well
beyond it.

The first is the default tolerance of `np.allclose`, which carries `atol=1e-8`.
A shallow grating's first-order efficiency is ~1e-9, three orders of magnitude
below that floor, so the default tolerance passed the comparison against any
value, including that of a solver returning zero. The first draft of the file
passed cleanly against the bug it was written to find. Any comparison of small
absolute quantities needs `atol=0` and a relative tolerance.

The second is a test geometry that is degenerate without appearing so. The same
first draft used α=10° at λ/p=0.36 and α=25° at λ/p=0.013. These are
reasonable-looking mounts, but at both the quantity under test lies 0.2% from
unity, so agreement was guaranteed whatever the solver did. The fix is an
assertion about the test setup in addition to the result:
`test_the_mounts_actually_discriminate` requires each mount to reach an order
where the two theories are at least 15% apart.

In both cases the test could not have failed, and a passing test is evidence
only if it could have failed. Non-vacuity assertions are inexpensive, and the
suite already uses them elsewhere; these were two further places where they
were needed.

---

## Névot–Croce factor applied in amplitude form to an intensity

`fresnel.reflectivity` multiplied $|e^{-2k_{iz}k_{tz}\sigma^2}|$ into $|r|^2$.
That factor is the amplitude form, and the intensity requires its square.

Nothing in the suite could detect the error. The factor is bounded by 1,
monotone in σ, and exactly 1 at σ=0, and "Debye–Waller over-damps near θ_c";
every property tested was satisfied by both forms. The limit that separates
them is $n \to 1$, in which the interface disappears, $k_{tz} \to k_{iz}$, and
the two models become the same expression:

| graze | NC (before) | DW | NC² |
|---|---|---|---|
| 5° | 0.963213 | 0.927771 | 0.927780 |
| 30° | 0.291214 | 0.084805 | 0.084806 |
| 85° | 0.007467 | 0.000056 | 0.000056 |

Roughened reflectivities were therefore too high: 0.732 where 0.536 was
correct, at 15° graze with σ=0.5 nm at λ=2 nm, a 37% error in the factor.

The error also falsifies a claim recorded in M15-E. The statement that
"Debye–Waller over-damps by ~1e-2 near/above θ_c, with no reliable sign far
below" was a measurement of this bug. With the intensity form the picture is
coherent and the sign is reliable in both regimes: Debye–Waller over-damps near
and above θ_c, under-damps deep below it by ~2e-3, and the two converge to 1e-5
far above θ_c. The old test looked below θ_c, where they do not converge.

---

## The cone opens away from the groove-axis vector

`conventions.md` §3 contains three statements that are individually true but
contradictory when taken together under a directed reading:

| Where | Statement |
|---|---|
| frame table | `ĝ = -ẑ` |
| wave vectors | every `k` carries `+\cos\gamma\,\hat{z}` |
| angles | γ is "the *polar* angle measured from the groove axis" |

Measured directly at the reference geometry (α = 25°, γ = 1.5°), the angle from
`+ẑ` is 1.5000° and the angle from `ĝ` is 178.5000°. Every propagating order
lies at exactly 1.5000° from `+ẑ` with an identical `k_z = 0.999657325`, so the
orders do lie on one cone, and that cone opens along `-ĝ`.

$$\hat{d} \times \hat{g} = \hat{n} \;\Longrightarrow\; \hat{x} \times \hat{g} = \hat{y}
\;\Longrightarrow\; \hat{g} = -\hat{z}$$

so `ĝ = -ẑ` is forced by the right-handedness relation, which comes from
Heilmann et al. (2024) §2.1, and is not a free choice. What was imprecise was
the wording of the γ bullet: γ is the acute angle to the groove-axis line, i.e.
the polar angle about `-ĝ`.

No computed result was ever wrong. Every consumer of γ (`facet_graze`,
`sin_beta`, `blaze_wavelength` and `Illumination`) takes it as the scalars
`sin γ` / `cos γ`, in which the direction of `ĝ` cannot appear.
`Illumination.direction_cosines` gets the sign of its `z` component right, but
nothing asserted it against `ĝ`. No code had ever needed a directed
groove-axis vector.

The situation parallels [the `t̂ = -d̂` finding](#the-profile-parameter-runs-backwards-and-nothing-had-noticed):
the inconsistency was invisible until something drew it. The first component to
need a directed `ĝ` is a 3D view of the diffraction cone, and a reversed vector
would open the cone into the grating instead of away from it, putting every ray
180° out.

The inconsistency is recorded and was not corrected by flipping `ĝ` to `+ẑ`,
which would break `d̂ × ĝ = n̂` and diverge from the reference. It is pinned by
`tests/test_illumination.py::TestTheConeOpensAwayFromTheGrooveAxis`, which is
placed where the bug would have occurred, as opposed to where it would have
been noticed.

---

## The profile parameter runs backwards, and nothing had noticed

`geometry.blaze_direction` returns $\beta_b = 2\delta - \alpha$. That value
requires the active facet's outward normal to lie at azimuth $+\delta$, so that
the facet descends with increasing $x$. Reflecting the projected incident
direction about the other orientation gives $\sin\beta = -\sin(2\delta + \alpha)$
instead: at the reference geometry ($\delta = 29.5°$, $\alpha = 25°$) that is
−84.0° where the blaze direction is +34.0°.

However, `Blazed.slope()` returns $+\tan\delta$ for $t \leq$ `apex`, so the
active facet rises with $t$. A facet that rises in $t$ while descending in $x$
implies $\hat{t} = -\hat{d}$.

No computed result was ever wrong. The scalar solver's blaze peak lands where
`blaze_wavelength` predicts:

| order | predicted $\lambda_b$ | solver peak | ratio |
|---|---|---|---|
| 2 | 4.0498 nm | 3.9970 nm | 0.987 |
| 3 | 2.6999 nm | 2.6863 nm | 0.995 |

(The residual arises because the sinc envelope peaks between discrete orders; it
is not a sign error.) `scalar.py` computes $\int e^{+i\Phi_m(t)}e^{-2\pi imt}\,dt$,
whereas §2/§3 imply $e^{-i\Phi_m}$ under $t \leftrightarrow +x$. The
substitution $t \leftrightarrow -x$ maps one integral onto the complex conjugate
of the other, and conjugation leaves $\|G_m\|^2$ bit-identical. The codebase has
therefore been self-consistent throughout.

The inconsistency remained invisible because no code had ever drawn the groove
profile and the diffracted rays in the same frame. `ProfilePlotPanel` plots
height against $t$ from left to right, which implies $t \leftrightarrow +\hat{d}$
without asserting anything physical, and every other consumer takes
$\|G_m\|^2$, in which the handedness cancels. The geometry visualizer is the
first component that has to place a facet relative to a ray, and with the
handedness reversed it would have drawn the beam striking the anti-blaze facet
and the blaze arrow leaving through the back of the active one.

The convention is recorded in [`conventions.md` §3](conventions.md); it was not
repaired by flipping the signs in `scalar.py`. That flip would equally preserve
the modulus, and every test would still pass, but it would move the module's
headline formula away from the form transcribed from the ISSI chapter and
Appendix D. Neither reference states a handedness (both quote only
$\|G_m\|^2$), and that absence is itself the finding; it is not an oversight to
be corrected.

---

## An exactly reciprocal solver produced an empty reciprocity report

CI failed on all three matrix jobs while the same commit passed locally.

The reciprocity check tracked its worst violation with a running maximum seeded
at 0.0 and updated only on a strictly greater value. A solver that is bitwise
reciprocal therefore never recorded the order from which the measurement came.
Locally, floating-point noise gave ~5 × 10⁻¹⁸ and the field was set; on CI the
BLAS returned exact zero and the report came back empty.

The more accurate the solver, the less complete the report, and the defect
appeared only on a different platform. It was fixed by collecting the
violations and taking the maximum at the end. The regression test constructs a
deliberately exactly-reciprocal solver, so that the zero path is exercised on
every platform and does not depend on rounding.

---

## Two corrupt files in the legacy collection

- `OGRE/AFM_test.ggp` line 164 has a missing x value, the result of a truncated
  write. The reader refuses the file instead of silently dropping the point,
  which would shift the profile and change every efficiency computed from it.
- `OGRE/pcgratetest3.txt` has no efficiency columns at all: it contains 20,308
  rows of scan geometry and no data.

Both rejections are correct behaviour, and both are encoded as tests so that
they are not mistaken for parser bugs.

### Three `.ggp` header variants, two of which the vendor does not load

| Header | Vendor accepts | Seen in |
|---|---|---|
| `3 0 - Polygonal type` + `Period:` | yes | most files |
| `# 3 0 - Polygonal type` (hashed) | **no** | `*_fix6/7.ggp` |
| `# X(normalized...)` only | **no** | `*_fix{,2,3,4}.ggp` |

`np.savetxt(header=...)` is what produces the hashed variants. This project
reads all three, since they are real files that need to be loaded, and writes
only the first, flagging the others via `GgpFile.format_valid`.

---

## Thesis and ISSI chapter reconciled

Cross-checking the two primary references shows that the grating equation, the
blaze direction $\beta_b = 2\delta - \alpha$, the blaze wavelength, and the
facet graze angle $\zeta = \arcsin[\sin\gamma\cos(\delta-\alpha)]$ all agree
exactly.

Three differences change what the code computes, and four are transcription
errors. All seven are tabulated in [`conventions.md` §10](conventions.md). The
two that matter most are:

- The thesis sawtooth phase omits `sin γ` and is inconsistent with the thesis's
  own Φ_b twelve lines later. The ISSI form is general and is the one
  implemented here.
- The thesis carries an obliquity factor and a Σ-renormalization that ISSI does
  not. Neither is a matter of notation, and neither was retained: both are
  errors, and both have been removed.

---

## The integral solver matches an independent integral-method code, and identifies the TASTE profile

This is the first comparison the project has produced between two rigorous
methods. The new perfectly-conducting integral solver was compared against the
TASTE wavescan in the reference corpus (`OGRE/tastetest_perf_wavescan.txt`, 541
wavelengths, 0.6–6.0 nm, γ = 1.25°, α = +19.99°), computed by an independent,
commercial implementation of the same method. Orderwise, across every
propagating order of the full scan, exactly one candidate profile agrees within
the reference's own energy-balance error, and every alternative, including the
design sawtooth and the mirrored candidates, is two to three orders of
magnitude further off. Comparisons against the commercial code are not
published here, because its licence terms have not been checked, so the
figures remain with the private corpus.

The comparison has three results:

1. **The unidentified `.ggp` is identified.** The roadmap's question of which
   `.ggp` goes with which run is answered for TASTE: the run was made from
   `AFM_test.ggp` and not from the design sawtooth. A margin of that size
   amounts to identification and cannot be attributed to tolerance slack. It is
   encoded as a positive test and a negative control in
   `tests/test_integral_corpus.py`.
2. **The solver stack is validated end to end**, covering the conical
   reduction, the kernels, the quadrature and the amplitude extraction, against
   an independent implementation of the same method on a measured profile. The
   project's own energy balance on that scan is 1.1e-8.
3. **The polarization mapping is pinned.** PCGrate's exported `Eff.TE` column
   corresponds to this project's groove-referenced TE at this mount. Each
   PCGrate table carries exactly one polarization, which the earlier scalar-era
   comparisons never had to take into account.

Two caveats apply. `AFM_test.ggp` still fails `read_ggp` (the truncated-write
defect above), so the corpus test parses around the damaged tail, with the
defect documented in place. In addition, the scalar-versus-reference deficit
reported earlier in this file ("the corpus can test the diffraction but not the
reflectivity") was computed on the idealised sawtooth, which this comparison now
shows the reference run never used.

---

## Scalar validity in a conical mount follows the reduced wavelength

The λ/period guard was testing the wrong ratio, and the flagship off-plane
regime passed a check that it should have failed by a factor of 46.

For a perfect conductor the conical problem decouples exactly into an in-plane
problem at the reduced wavelength λ/sin γ, which is the fact on which the
integral solver is built. The condition that the structure be much larger than
the wavelength, which Kirchhoff theory requires, must therefore be judged on
λ/(p·sin γ). The old guard used λ/p, which at the reference geometry (period
315.15 nm, γ = 1.25°, λ = 2.226 nm) reads a comfortable 0.007 while the
effective ratio is 0.32.

The following was measured for scalar against the integral solver, both
relative, TE, `Sinusoidal(0.05)`, with λ/p held fixed at 0.001 while the cone
closes:

| γ | λ/(p·sin γ) | max abs scalar-vs-integral |
|---:|---:|---:|
| 10° | 0.0058 | 1.6e-3 |
| 5° | 0.0115 | 2.4e-3 |
| 2.5° | 0.0229 | 3.4e-3 |
| 1.25° | 0.0458 | 4.5e-3 |

The disagreement grows monotonically with 1/sin γ at fixed λ/p, while the naive
ratio does not change. At the reference geometry itself (Blazed 29.5°/70.5°,
reduced ratio 0.32) the orderwise disagreement is severe: the integral solver
puts 0.42 in m=+2 and 0.21 in m=+3, where bare scalar puts 0.69 in m=+3.

This has two consequences. First, the guard now tests the reduced ratio (the
in-plane case is unchanged), and it fires on the project's own flagship regime
for λ ≳ 0.7 nm. That is the correct outcome and not collateral damage. Second,
the long-standing per-order disagreement between scalar and reference ("the
corpus can test the diffraction but not the reflectivity", above) is no longer
unexplained and becomes a statement about regime. Scalar theory remains useful
there for order positions and the blaze envelope, but per-order numbers require
the integral method.

Alongside the γ ladder, the committed cross-method suite
(`tests/test_cross_method.py`, the first scalar-versus-integral tests that need
no external corpus) measures the other two failure axes on `Sinusoidal`
profiles. Depth_fraction {0.01, 0.02, 0.05, 0.2} at λ/p = 0.04 gives maximum
absolute differences of {3.6e-4, 8.4e-4, 4.4e-3, 1.0e-1}, so depth is the
dominant axis, with a ×285 increase over a ×20 depth range. λ/p {0.05, 0.1, 0.4}
at depth 0.02 gives {8.4e-4, 2.2e-3, 4.1e-3}. In the overlap regime (shallow
grooves and a small reduced ratio) the two solvers agree to a few 1e-4 absolute
in both mounts, and the integral solver's own energy-balance deviation is below
1e-7, so none of the gap is mesh error.

---

## The facet-normal test misses cast shadows, almost entirely on the exit side

The visibility masks tested only whether a point's own facet faces the ray. The
shadow the apex casts across the trough was invisible to them. That shadow
covers a fraction of a percent of the period toward the incident beam, but
10–50% of it toward individual diffracted directions, and it changes the blaze
order by −32%.

The geometry is as follows. For a groove-invariant surface, 3D occlusion reduces
exactly to a 2D test in the transverse plane, and γ drops out: the ray's
groove-parallel component slides along the invariant direction and can neither
create nor remove an occlusion. The scan therefore runs at the azimuth alone,
α for incidence and β_m for each order. On an ideal sawtooth (blaze slope s_b,
anti-blaze slope s_a, apex at t_a) the cast shadow past the trough has the
closed form

    Δ = (1 − t_a)(s_a − s_r)/(s_b + s_r),   s_r = cot θ

which is zero until the ray dips below the anti-blaze slope (θ > 19.5° for the
29.5°/70.5° reference groove). The table gives the extra fraction of the period
masked beyond the facet-normal test, measured on the reference geometry:

| direction | θ | extra masked |
|---|---:|---:|
| incident | 19.99° | 0.38% |
| exit m=+3, λ=2.226 | 39.0° | 14.7% |
| exit m=+2, λ=3 | 32.1° | 9.5% |
| exit m=+2, λ=4 | 55.3° | 28.3% |

On the incident side the shadow is almost entirely captured already, since
α = 19.99° lies just past the 19.5° threshold. This is presumably why the
omission went unnoticed: the `shadowed_fraction` a person would check barely
changes. The exit directions are much steeper, and there the omitted shadow is
large. In the full pipeline (`visibility="horizon"`, Au, local model) the blaze
order at its blaze wavelength moves from 0.514 → 0.348 (−32%), a change of the
same order of magnitude as the M16 groove-resolution change.

Reciprocity survives by construction, because occlusion along a straight ray is
the same from either end (measured at 1e-16 with and without coating). The
horizon test also subsumes the orientation test entirely: a back-facing point
has a falling ray-adapted height u = g + p·cot θ·t, which places it under the
running horizon.

Two sign traps, each of which fails silently by reporting no cast shadow
instead of raising an error, are pinned by tests
(`tests/test_geometry.py::TestHorizonVisible` derives the sliver trigonometry
independently). The first is that the profile parameter runs against the
dispersion direction ("The profile parameter runs backwards", above), so a ray
with α > 0 travels toward +t and its occluders lie at smaller t. The second is
that the scan must use `problem.height_nm` and not the period-normalised
`profile.height`, because the ray's run-to-drop ratio compares height against
period directly.

Cast shadowing moved the blaze order from 0.69 toward the rigorous 0.21 of the
(out-of-regime) perfect-conductor comparison. The direction is correct, but this
is suggestive and does not constitute validation. As with the geometric-mean
weight, the horizon model ships as opt-in until a finite-conductivity A/B
comparison exists.

---

## s and p are degenerate at the grazes this project runs

Every scalar reflectivity model averages facet-local s and p 50/50, because the
groove-TE/TM → facet-s/p rotation is deliberately left unowned. The
measurements show that this simplification costs under half a percent in the
working regime, so the mapping remains unimplemented, and the solver now
quantifies the cost of the simplification instead of assuming it is small.

The fractional intensity split |R_s − R_p|/R_s on the vendored Au table is:

| λ (nm) | ζ = 0.5° | 1.25° | 2° | 5° |
|---:|---:|---:|---:|---:|
| 0.62 | 1.2e-4 | 3.9e-4 | 2.2e-3 | 2.8e-2 |
| 1.0 | 3.3e-4 | 8.9e-4 | 1.7e-3 | 2.4e-2 |
| 3.0 | 1.9e-3 | 4.9e-3 | 7.9e-3 | 2.4e-2 |
| 6.0 | 2.2e-3 | 5.5e-3 | 8.9e-3 | 2.5e-2 |

At amplitude level, |r_s − r_p|/|r_s| ≤ 9e-3 and arg(r_p/r_s) ≤ 8e-3 rad at
ζ ≤ 2°. Since the unpolarized mean is already computed, any rotation of the
polarization basis can move the result by at most half the split, which is
< 0.5% at every graze the reference groove presents (the worst local graze is
1.94°). The mapping begins to matter only above ~5° facet graze, that is, for
steep grooves near normal incidence, where the Brewster warning already fires.

Implementing the rotation would therefore add false precision here. What was
missing was an estimate of the cost. The provenance now evaluates the split at
the guarded graze and at the scan's wavelength extremes, and warns above 2%, so
that a steep-groove near-normal run (measured split: 40%) reports that its
unpolarized average is a modelling decision with a measurable cost.

---

## Sub-cell shadow boundaries: three schemes compared

Binary visibility masks limit the masked models to O(1/n) convergence, with an
error of ~5e-4 at the default 2048 points on the coated blazed reference. The
scheme that resolved this places every boundary geometrically and reaches 2e-7
at the same 2048 points. The two schemes that failed are recorded because each
appears to be a reasonable design until it is measured.

1. **Interpolating the mask function across the straddling pair** (the obvious
   design) mishandles discontinuities. At a profile corner the orientation
   function `sin ζ(t)` jumps instead of crossing zero, so the linear estimate is
   systematically biased, and on the reference sawtooth it assigned 0.8 of a
   cell of unit-amplitude weight to a back-facing sample. As measured, it
   performs worse than binary masks on blazed profiles (1.5e-3 vs 5.0e-4 at
   n = 2048).
2. **Extrapolating the horizon deficit from the shadow side, at sample level**,
   fixes the reappearance boundary but not the entry boundary. A shadow's level
   was the occluding apex's height at a sample, and the true corner between
   samples lies up to one cell higher, which shifts both boundaries of the
   stretch by O(h). The measured lit-fraction errors of 2e-4–9e-4 plateau across
   n before dropping, which is the signature of a boundary mislocated by a fixed
   physical length.
3. **The surviving scheme** (`geometry.horizon_weights`) recovers the entry
   corner by intersecting the secants through the two samples on each side of
   it. The lit-side pair must be (i−2, i−1), because the last discretely-lit
   sample can already lie past the corner. Every shadow stretch is occluded by
   its own entry crest (a lit point is at its running maximum), so the exit is
   where u climbs back to that corner's height, one period-drop lower if the
   stretch wraps t = 0. On a polygon (Blazed, Lamellar, and every measured
   groove) the lit fraction is then exact to machine precision at any n
   (asserted at 1e-13 in `tests/test_geometry.py::TestHorizonWeights`).

One composition rule mattered as much as the boundaries: incident and exit
masks are combined by taking the minimum, not the product. Every direction's
shadow starts at the same apex corner, so the two masks share boundary cells,
and multiplying their sub-cell weights counts the shared boundary twice. This
was measured as the entire remaining error, concentrated on the blaze order,
whose integrand is nearly constant, so that its efficiency is the square of the
lit overlap.

The final measured ladder, for the coated blazed reference with
`visibility="horizon"` against a 131072-point run, is 4.4e-6 (n=1024), 2.0e-7
(2048), 7.8e-7 (4096), 1.0e-5 (8192) and 1.8e-7 (16384). The residuals are
small and erratic, and every rung lies below the binary scheme's n=2048 error
by 30× or more. The `facet-normal` mode keeps binary masks deliberately: its
boundary function jumps at corners where no crossing estimate is sharp, and
bit-for-bit reproducibility is the purpose of that mode. Orders away from the
blaze converge spectrally (1e-8 at n=2048), which shows that the shadow
boundaries were the only obstacle between the masked models and the
convergence of the pure phase integral.

---

## Finite conductivity, with signs fixed by Fresnel checks

The M20 coupled conical system (Goray & Schmidt 2010) was derived independently
in project conventions instead of being transcribed. Their normal points into
the metal and their potentials carry a factor 2, which is the kind of mapping in
which a copied sign error survives every syntax check. In place of trust in the
transcription, the signs were fixed by the following checks, with the measured
outcomes:

| Check | Measured |
|---|---|
| Flat interface, in-plane: TE amplitude vs `fresnel.amplitude` r_s, TM vs r_p, phase included | agrees to <1e-4, all of gold/X-ray/glass indices, two angles |
| Flat interface, conical (cos γ = 0.8): full 2×2 (E_z,B_z) reflection matrix vs the rotated-Fresnel closed form | 6.2e-5, cross-coupling terms included |
| R + A = 1 with A from the independent boundary integral (their eq. 26, re-derived; reduces exactly to 1−|r|² on a flat interface, verified analytically before any code ran) | 1e-5 flat; ~1e-4 smooth sinusoid; 1.2e-4 for Au at 2 nm, 1.5° graze, P = 256 |
| R + T = 1, lossless substrate; boundary integral equals the transmitted Rayleigh sum | few 1e-4; A−T ≤ 9e-5 |
| Perfect-conductor limit, n = 100i, in-plane and conical | 1.3e-3 (TE) / 7.4e-4 (TM) worst order; the residual is the ~1 nm field penetration, i.e. physics, not error; absorption at numerical zero |
| G&S Table 3 (dielectric sinusoid, conical, pure E_z; cross-method against Li's CM) | within 6e-5 of every tabulated order, reflected and transmitted, at half their collocation count |
| Lorentz reciprocity, lossy conical Au at X-ray graze | ~5e-5, through `check_reciprocity` unchanged |

Two facts about the regime were established:

1. **Below the critical angle the metal side sets the mesh.** The field decays
   into the metal on the λ/√(2δ) scale. For Au at 2 nm and 1.5° graze that
   scale is ~21 nm, against a 76 nm vacuum reduced wavelength, so the metal-side
   floor (~210 points on the 600 nm test sinusoid) is the binding one. The
   solver enforces it with the same 6-nodes-per-wavelength refusal.
2. **The absorbed fraction is large where it matters.** For Au at 2 nm and
   1.5° graze on the test sinusoid, R = 0.70 and A = 0.30. The absorbed
   fraction is the quantity the perfect-conductivity mode could not produce; it
   is now recorded per wavelength on the scan and consumed by
   `check_energy_balance` automatically.

---

## Goray & Schmidt's Table 4 is a publisher's duplication of Table 3

This was noticed while transcribing validation data. In the published PDF,
Table 4 ("dielectric sine grating for E_z = 0") prints exactly the same
efficiency and polarization-angle values as Table 3 ("for B_z = 0") in all
twelve rows and all four digit columns, with only the footnoted polarization
state differing (δ = 8.499°, ψ = 180° vs δ = 81.501°, ψ = 0). This was
confirmed against the rendered page as well as against the extracted text.

The two incident states cannot physically give identical efficiencies, and the
solver resolves the question: the true E_z = 0 case differs decisively
(specular 1.07% against the printed 10.33%; T₊₁ 16.7% against 9.9%) while
conserving R + T = 1 to 2e-5. A test in `tests/test_integral_finite.py` pins the
disagreement so that the finding survives refactoring, and the computed E_z = 0
values in the test file serve as the project's record of that case.

A note on decoding, for anyone returning to their Tables 1–2: the footnote
angles follow Li's convention, with tan δ = |E_s/E_p| and ψ the relative phase.
This is recoverable because the Table 3/4 footnotes bracket the pure-z states
(81.501° + 8.499° = 90°, and the same split is reproduced by this project's
incident basis rotation). Table 5's δ = 0 and Table 6's δ = 90° at θ = 0 are
exactly the pure E_z and pure B_z states.

---

## A nominal tip does not explain the rounded groove

Eroding the TASTE scan with the nominal probe changes nothing it measures, and
burying an ideal facet to the measured depth requires a tip forty times blunter
than nominal. The rounded-groove finding above left open whether the grating or
the tip is responsible, and the tip-deconvolution milestone (M22) narrows the
question: a tip in anything like its specified condition is excluded. What
remains is either a severely worn tip or a grating that is genuinely rounded.

The scan is the same as in the finding above
(`TASTE_ALS_A205_Ti_Pt_flatten.txt`, 512 px over 2 µm, i.e. 3.91 nm/px), and so
is the pipeline, now with `tip_correction='erosion'` and the probe modelled as a
cone of 18° half-angle capped by a spherical apex (`core/tip.py`, Villarrubia
1997). The nominal probe is R ≈ 1 nm ("2 nm wide"), and wear only increases R.

In the reconstruction direction, the measured image is eroded. This can only
deepen a groove that the tip failed to reach into:

| apex radius R | depth/period | implied blaze | pixels certain |
|---|---|---|---|
| uncorrected | 0.3275 | 20.33° | n/a |
| 1 nm (nominal) | 0.3275 | 20.33° | 97.7% |
| 2 nm | 0.3275 | 20.33° | 95.9% |
| 5 nm | 0.3275 | 20.33° | 81.8% |
| 10 nm | 0.3274 | 20.32° | 74.8% |

The correction recovers no depth. The measured surface is already reachable by
all of these tips, so there is no depth hidden by the tip for erosion to
return. A surface whose troughs had been rounded by the tip would come back
deeper, and this one does not change at the reported precision.

In the forward direction, an ideal sharp `Blazed(27.91°, 70.5°)` groove (the
fitted facet angle, at the measured 314.1 nm period, sampled at the scan's own
pitch) is dilated with the same tip family, and the depth shown by the
resulting image is examined:

| apex radius R | image depth/period | implied blaze |
|---|---|---|
| sharp truth | 0.4460 | 27.91° |
| 1 nm (nominal) | 0.4439 | 27.77° |
| 2 nm | 0.4433 | 27.74° |
| 10 nm | 0.4316 | 27.00° |
| 20 nm | 0.4176 | 26.11° |
| 40 nm | 0.3909 | 24.40° |
| **80 nm** | **0.3362** | **20.89°** |

The nominal tip costs 0.14° of implied blaze, against an observed gap of 7.6°.
Reproducing the measured depth requires R ≈ 80 nm, a probe worn to forty to
eighty times its specified apex. The 18° flank is not responsible either: at
72° from the surface it is steeper than the 70.5° anti-blaze facet, so the facet
remains reachable at any radius, and the entire deficit comes from the apex.

The remaining question could be settled by the probe's actual condition at scan
time (an 80 nm apex is far outside specification and is often visible in the
vendor's tip check), by a re-scan with a fresh tip, or by a cross-section. Until
one of these exists, the rounding should be treated as plausibly real, and the
consequence for absolute efficiency stated in the finding above stands
unchanged, because erosion with any defensible tip model does not sharpen this
groove.

Note added 2026-10-22: the probe is a pyramid, and one of the reasons given
above does not hold along the cantilever. Bruker quotes the ScanAsyst-Air at
R = 2 nm nominal and 12 nm maximum (a radius and not a width, so "R ≈ 1 nm"
above understates it; the 2 nm rows already cover the nominal case), with a 15°
front face, a 25° back face and 17.5° sides, each ±2.5°. The tables above were
made with an 18° cone, which matches the sides, so they describe a scan across
the cantilever. The statement "The 18° flank is not the culprit ... the facet
remains reachable at any radius" does not hold along the cantilever: the back
face rises at 65°, below the 70.5° facet, and buries it. The forward test was
repeated with that pyramid, using the same ideal groove and pitch:

| tip, flank facing the anti-blaze facet | image depth/period | implied blaze |
|---|---|---|
| 17.5° sides, R = 2 nm | 0.4435 | 27.75° |
| 25° back face, R = 2 nm | 0.4224 | 26.41° |
| 25° back face, R = 12 nm | 0.4097 | 25.61° |
| 27.5° back face (spec limit), R = 12 nm | 0.4014 | 25.08° |
| 27.5° back face, R = 40 nm | 0.3694 | 23.03° |

The worst tip within the specification, in its worst orientation, therefore
accounts for 2.8° of the 7.6° gap (1.8° of it from the orientation alone, at
12 nm), and the conclusion stands: no tip within Bruker's specification
reproduces the measured 20.3°. It is still unknown how the TASTE scan was
oriented relative to the cantilever, which determines whether the orientation
contribution applies at all. Erosion of the real scan with the along-cantilever
presets has not been run, because the scan is not in this checkout.

The behaviour is pinned by `tests/metrology/test_tip.py::TestTheRealScan`, which
skips without the group's data and otherwise fails if nominal-tip erosion ever
starts to change the measured depth, which would mean either that the scan had
changed or that the erosion had started to invent surface.
