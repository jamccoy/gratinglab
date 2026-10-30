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
7.6°.

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

### Five `.ggp` control-string variants on disk, and what PCGrate documents

| Control strings | Seen in |
|---|---|
| none, points only | Aditi's files, and the default output |
| `3 0 - 20 harmonics - Polygonal type` + `Period:` | PCGrate's own 35 polygonal exports |
| `3 0 - Polygonal type` + `Period:` | legacy files written by this project |
| `# 3 0 - Polygonal type` (hashed) | `*_fix6/7.ggp` |
| `# X(normalized...)` only | `*_fix{,2,3,4}.ggp` |

`np.savetxt(header=...)` is what produces the hashed variants.

**Answered from the vendor manual, 2026-10-30.** The *User's Guide of
PCGrate-S(X) Series, version 6.1* authors a polygonal `.ggp` step by step in
10.2.2 (pp. 62–63) and lists the thirteen format errors the program rejects in
10.2.3 (p. 63). The first control string is `3 0 - 20 harmonics - Polygonal
type`, verbatim, and all 35 polygonal files the program ships carry exactly
that line. This project wrote `3 0 - Polygonal type`, which omits the
`- 20 harmonics -` segment, so `write_ggp(header=True)` produced a first
control string a direct load rejects; and `_classify` compared against that
short form, so every one of the vendor's own files read back as
`commented-canonical` with `format_valid` False, despite containing no `#` at
all. Both are fixed: the documented pair is what is written and what counts as
`canonical`, the short form is labelled `short-type-line` and flagged, and the
geometry rules of 10.2.3 are now asserted against the vendor's files and
against this writer's output in `tests/test_ggp.py`.

The guide also settles why the headerless practice works: a file with no
control strings, or with control strings in another format, is explicitly legal
input — through the Border Profile Editor, which exports it into the program,
rather than a direct load. Aditi loads these through the border-profile import
and reports it does not want a header, so the points are still written alone by
default.

Nothing numerical moved. No solver and no comparison gates on `format_valid` or
`header_variant`; only the tests read them, and the point data parsed correctly
under both labels.

Two things the guide does not get quite right about its own program: 24 of the
35 shipped files separate the columns with a tab, where 10.2.2 says to use the
space key, and all 35 end in a line break, where a trailing one was reported
(2026-10-27) to be rejected. The reader accepts tabs; the writer still omits
the trailing break.

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

## The integral solver solved every blazed grating as its mirror image

Found 2026-10-24 by the gratinglab app, whose JavaScript port of the scalar
solver was being checked against this package's integral method on a title
figure. At `Blazed(9)`, p = 1667 nm, 550 nm, α = ±40°, unpolarized, perfect
conductor, 400 graded points:

| α | integral (old) | integral (fixed) | scalar |
|---:|---|---|---|
| +40° | m=−1 0.56, m=0 0.17 | m=+1 0.80, m=0 0.17 | m=+1 0.89, m=0 0.10 |
| −40° | m=−1 0.80, m=0 0.17 | m=+1 0.56, m=0 0.17 | m=+1 0.27, m=0 0.10 |

`blaze_direction(9, 40)` is −22°, which is where m = +1 goes, so the scalar
solver and the documented convention agreed and the integral solver was the
odd one out. Its old answer at +α is its new answer at −α with the orders
negated, digit for digit: a reflection, not an error of accuracy.

**The cause** is one line. `physical_boundary` placed the curve at
`x = t P`, while [`t̂ = −d̂`](#the-profile-parameter-runs-backwards-and-nothing-had-noticed)
puts it at `x = −t P`. That entry says no numerical result depends on the
direction, which is true of `|G_m|²` and was the reason it went unnoticed, but
the integral method does not take `|G_m|²`: it draws the boundary. It now
places it at `x = (1 − t) P` with the nodes reversed and `nx` negated, which
keeps them advancing in `+x` as every operator assumes.
The corner-graded mesh's polygon is reflected with them, and its vertex list
has to start at the image of the closing vertex, not merely be reversed: the
wall that closed the profile's list now opens the curve, and left at the end it
falls outside the three periods the corner subtraction integrates, which put a
0.07 offset on every node. Two tests of ungraded TM corners had been tuned
against the mirrored groove; their incidence angles change sign so they still
measure the case they were written for, since the new solver at −α is the old
one at +α exactly.

**Nothing caught it because every cross-method case was a sinusoid**, and a
mirrored sinusoid is the same sinusoid half a period along. The guard is
`tests/test_cross_method.py::TestAgreementOnAnAsymmetricGroove`: a 3° facet at
λ/p = 0.05, in-plane and off-plane, where integral and scalar agree to 0.017
and both put about 0.96 in the order `blaze_direction` predicts. On the old
boundary the same test fails by 0.97.

**The PCGrate corpus pinned the other half.** The TASTE match fed the `.ggp`'s
`x` straight in as `t`, on the old boundary, and matched orderwise to 7.8e-4 at
α = +19.99° with the mirrored candidates far off. Since α and the order labels
were recovered in this project's frame, that match says the `.ggp` runs along
`+d̂`. So the adapters reflect too (`FromProfileData.along_dispersion`), and the
fixed solver on the reflected points reproduces the old solver on the raw
points to 2.2e-14, measured on an asymmetric measured-style groove in the
TASTE mount (TE, equal-arc and graded) and in-plane (TM). The corpus test
could not be run when this was fixed, because the corpus was on no machine to
hand; that identity is what carries its result across, and it should still be
run. Its negative control, the design sawtooth, is now solved the other way
round from when its margin was measured, so its 0.05 bound is unverified.

The scalar solver, by contrast, **had** been solving every `.ggp` and every
AFM profile mirrored relative to PCGrate, since it honoured `t̂ = −d̂` and the
files never passed through a reflection. That is the half of this fix that
changes results on measured data.

**Two earlier conclusions in this file rest on the mirrored solver.** At the
reference geometry (`Blazed(29.5, 70.5)`, p = 315.15 nm, γ = 1.25°,
α = 19.99°, λ = 2.226 nm, TE):

| | m=+1 | m=+2 | m=+3 |
|---|---:|---:|---:|
| integral, old | 0.05 | 0.42 | 0.21 |
| integral, fixed | 0.13 | 0.06 | 0.76 |
| scalar | | | 0.69 |

So the "severe" orderwise disagreement quoted under
[scalar validity in a conical mount](#scalar-validity-in-a-conical-mount-follows-the-reduced-wavelength)
was mostly the reflection, and the cast-shadowing note that the horizon model
"moved the blaze order from 0.69 toward the rigorous 0.21" now reads the other
way: the rigorous value is 0.76, above bare scalar. The shadowing comparison has
been re-run, in gold as well, and its result is recorded where the claim was made,
under [the facet-normal test](#the-facet-normal-test-misses-cast-shadows-almost-entirely-on-the-exit-side).
The reduced-ratio guard itself rests on the sinusoid ladder, which a reflection
cannot touch, and stands.

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

**That reading was wrong, and so was its direction.** The 0.21 came from the
integral solver drawing the groove mirrored (see
[the entry on that](#the-integral-solver-solved-every-blazed-grating-as-its-mirror-image)).
Re-run on 2026-10-25 with the fixed solver, and with the finite-conductivity
A/B the paragraph above was waiting for, at the same geometry (γ = 1.25°,
α = 19.99°). The integral is graded to order 4 and agrees with itself between
512 and 1024 points to 2e-4 or better throughout:

| λ, blaze order | bare scalar | horizon | integral |
|---|---:|---:|---:|
| Au, unpolarized, 2.226 nm, m=+3 | 0.514 | 0.348 | 0.574 |
| Au, 3.0 nm, m=+2 | 0.415 | 0.343 | 0.637 |
| Au, 4.0 nm, m=+2 | 0.253 | 0.165 | 0.253 |
| Au, 5.0 nm, m=+1 | 0.334 | 0.325 | 0.748 |
| perfect, TE, 2.226 nm, m=+3 | 0.688 | 0.466 | 0.763 |
| perfect, TE, 3.0 nm, m=+2 | 0.573 | 0.452 | 0.851 |

The rigorous blaze efficiency is at or above bare scalar at every wavelength
here, so the horizon model never moves toward it. At the blaze wavelength in
gold it widens the gap from 0.06 to 0.23. This does not fault the shadow
geometry, which the closed-form tests above pin; it faults removing the
shadowed power from a scalar result that is already short, its orders summing
to 0.35 to 0.54 against the integral's 0.74 to 0.76. The geometry is also
outside scalar's regime (reduced ratio 0.32), so it says nothing about the
horizon model inside it. Opt-in remains right, and the A/B it was waiting for
now argues against it rather than for it.

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

---

## A graded corner mesh, and three errors only a solve could find

Measured 2026-10-23 while adding `corner_grading` to the integral solver. The
mesh itself works: on three cornered profiles the TM energy-balance deviation at
P = 512 fell from 1.3e-3, 2.9e-2 and 1.4e-2 to 1.8e-6, 6.6e-6 and 4.2e-6
(`theory/integral.md` §5). What is worth recording is how it got there, because
each operator passed a direct accuracy test while being wrong.

**The principal-value diagonal was wrong at the nodes next to a corner.** K, L
and D_t V take their diagonal as the average of the kernel at the two adjacent
nodes. A node next to a corner has one neighbour across it, where the kernel is
of order 1/distance, so the average is not the limit along the node's own facet.
Measured on a lamellar at P = 256: those diagonal entries read 0.04 to 0.64
against an identity term of 1/2, where every other node reads zero. Taking the
limit one-sided, from the neighbour on the same side, brought the graded TM
results in the table above down by a further factor of 2 to 10. **The
equal-arc mesh has the same flaw, more weakly (0.02 to 0.1), and it is not fixed
there**, so that the default results stay exactly what they were. It is a
candidate contributor to the first-order corner behaviour that this mesh was
built to remove.

**A comparison function applied to the wrong factor.** The tangential layer
removes its Cauchy singularity with a Hilbert comparison applied exactly in
Fourier space. Written in τ and applied to the density σ, it was exact on every
smooth test density, and `D_t V σ` matched the derivative of `V σ` to 4e-4
mid-facet. In the finite-conductivity solve it sent the absorption to 5e5, with
the reflected sum still a sane 0.46. A corner density is singular in arc length
and, through the grading, unbounded in τ, and the spectral transform of an
unbounded function is garbage. Applied to σ s′, which the grading makes smooth,
the same solve gives R + A within 4e-4 of 1. The lesson is the one
`test_integral_operators.py` already states for `D_t`: an operator verified on
smooth data says nothing about the densities a solve actually produces.

**The energy balance passes on a TE result that is wrong.** On the ideal
sawtooth the equal-arc TE balance reads 2.6e-10 at P = 512 while order −1 sits
6e-5 from the graded answer and moves by 2e-4 between rungs. The theorem holds
for any density solving the discrete system, so it is blind to an error common
to every order. The graded mesh settles that order to 5e-6. Reciprocity, which
is not the same theorem, did see the TM corner error: 2.3e-2 on the equal-arc
sawtooth at P = 256, 1.3e-4 graded.

**What is not resolved.** With finite conductivity the graded mesh helps TE
everywhere and the conical X-ray blaze greatly (8.3e-2 to 4.1e-4), but TM on a
visible-light metal (n = 0.2 + 3.5i) stalls near 4e-3 at order 4 and diverges at
order 6, and the conical X-ray lamellar is erratic on both meshes. The solver
therefore refuses `corner_grading` with `conductivity="tabulated"`. Two
candidate causes, neither tested: Kress's graded Nyström proof needs the
contributions of the nodes nearest the corner modified, which is not done; and
a negative permittivity gives a corner oscillatory edge singularities that no
grading resolves.

Pinned by `tests/test_integral_graded.py`.

## Finite conductivity on the graded mesh: two errors that lived next to the corners

Measured 2026-10-24, finishing the milestone above. Both causes the previous
entry named were wrong, and the two real ones are the same mistake: an error
confined to the few nodes nearest a corner is harmless in the perfectly
conducting TM equation, which is second kind, and is not harmless in the
finite-conductivity system, which is first kind in its densities.

**The negative permittivity was not the problem.** The TM corner exponent at
ε = n² ≈ −12.2 + 1.4i (n = 0.2 + 3.5i) is ν = 0.600 − 0.009i at a 90° or 270°
corner and 0.564 − 0.010i at the blaze's 80° apex, against 2/3 and 0.643 for a
perfect conductor, and real to two decimals. The singularity is a little stronger than the
perfect conductor's and resolved the same way.

**K⁺ near a corner.** The trace block is `(I/2 + K⁺)V⁻`. For a node a distance
a from a corner, the double-layer kernel over the facet across it varies on the
scale a, and the graded nodes on that facet are spaced on the same scale, so the
rectangle rule's error is a fixed fraction of the integrand at the first few
nodes, whatever P is. Subtracting the density's value at the target and
supplying the Laplace part of `∫ K` exactly, as subtended angles of the polygon's
segments, removes it. The visible-metal TM lamellar went from a stall near 3e-3
to 3.4e-6 at P = 512, and the X-ray lamellar on the graded mesh from 1.2 at
P = 256 (worse than equal arc) to convergence.

**The prototype of that fix had its own corner error, and only the blaze showed
it.** The target's own segment contributes zero, and the first version found it
by testing collinearity. On the lamellar every node coordinate is exact, because
the facets are axis-aligned. On the blaze a node 1e-9 nm from a vertex carries
rounding of 1e-13 nm, the collinearity test missed, the own segment counted as
a whole 1/2, and the visible blaze read 5.9. The segment is now the one nearest
the node. `test_corner_subtraction_completes_the_rule` checks the subtraction
against an independent adaptive quadrature on all three profiles, and it fails
when the correction is negated or removed (checked by doing both).

**D_t V⁻ near a corner.** After the first fix every case converged, at third
order, except the conical X-ray lamellar, which wandered at p = 6 below P = 320
and stalled near 7e-3 at p = 5 all the way to P = 768. The same lamellar in-plane
converged to 1e-6 at every order, which pointed at the one operator only the
conical system uses. On solved densities its weighted disagreement with the
spectral derivative of the trace was 0.22 against a total of 11 at p = 5, and
largest at wall nodes 1e-3 to 1e-1 nm from a corner. Replacing it by that
composition gives under 1e-6 at P = 768 at every order from 4 to 6. **That reverses a rule this
code states**: `D_t` composed with an assembled `V` was measured useless on the
equal-arc mesh, which is why the single operator exists. On the graded mesh the
trace is smooth in τ, so the composition is the accurate one, while the single
operator's Cauchy remainder is not continuous across a corner. Each mesh now
uses the operator that is right on it.

**Short sides make it worse, and the reason is the allocation.** Nodes go to
sides in proportion to length, so by the allocation's arithmetic the first
node on a 15.8 nm wall sits about 10^(p−1) times farther from its corner than the first node on the 157 nm floor
beside it. An allocation proportional to `L^(1/p)`, which matches the two scales,
made the X-ray lamellar converge smoothly even with the old operator, and cost
accuracy on the 1000 nm blaze (8.2e-4 at p = 6, P = 512). It was a diagnosis,
not a change, and the mesh is as it was.

**Re-measured with both fixes in, the allocation should stay proportional to
length.** Shares proportional to `L^(1/p)` and to `L^(1/2)` were tried against
`L` on all four finite-conductivity cases and the three perfectly conducting
profiles of §5, at p = 4 and 6. At P = 512 and p = 4 the alternatives lose on
the X-ray blaze (2.9e-6 against 2.2e-5 and 1.5e-5), the X-ray lamellar (2.0e-6
against 5.8e-6 and 3.9e-6), the visible blaze (2.4e-5 against 3.7e-5 and
3.1e-5), the perfectly conducting lamellar (1.8e-6 against 8.1e-6 and 1.9e-5)
and the ideal sawtooth (6.6e-6 against 1.5e-4 and 1.4e-4). They gain only on
the visible lamellar, by 15%. The short-wall trouble was the old `D_t V⁻`, and
with it gone the scale mismatch at a corner costs nothing measurable.

**Two things the energy balance could not see.** The single operator balances
the X-ray lamellar to 1.1e-6 at P = 1024 and the composition to 8e-7 at P = 768,
and at P = 1280 their efficiencies differ by 2.9e-5. Only the ladder of efficiencies settles it: the old
path moves by 1.9e-5 between P = 1024 and 1280, and the new one by 1.1e-6. And
at p = 6 past P = 512 on the 1000 nm blaze the solve degrades to 6.8e-2, which
looked like a third failure and is the regime the separation floor already
refuses: the perfectly conducting solve is refused there for the same mesh.

**Not changed.** The equal-arc mesh is untouched in both boundary conditions,
and the conical cornered case is still poor there (0.08 to 0.16 at P = 512). The
`K` subtraction was tried on the perfectly conducting TM equation and left off:
it helps the lamellar and hurts the ideal sawtooth by up to 25 times.

Pinned by `tests/test_integral_graded.py`; the ladders are in
`theory/integral.md` §5.

---

## What a cross-section shows that a scan cannot

Three cleaved cross-sections of the 160 nm OGRE laminar gratings say that an
etched trench narrows with depth by a third of its width, that a resist line
before transfer does not, and that both the crest and the floor meet the
sidewall on a radius. None of those three facts is available from an AFM scan
of the same surface, and the first two are not representable at all by
`Lamellar`, which was until now the only binary groove this package offered.

The micrographs are Zeiss SmartSEM TIFFs read by
[`metrology/io/sem.py`](../src/gratinglab/metrology/io/sem.py), which takes the
scale from the instrument's own record rather than from a scale bar. Each row
below is one image, measured by the same procedure: the period from the
strongest spatial frequency along the rows, the depth from the extent of the
rows that carry it, and the trench width at each depth from the dark fraction
of that row.

| image | nm/pixel | period | depth | trench width, mouth → floor | sidewall |
|---|---|---|---|---|---|
| `ogre_xsec/xsec_050` | 2.233 | 163.3 nm | 69.2 nm | 0.70 → 0.42 | 72° |
| `FESEM_picks/xsec_058` | 5.582 | 158.8 nm | 67.0 nm | 0.63 → 0.32 | 70° |
| `110_2ionmillKOH_xsec/xsec_04` | 2.484 | 169.5 nm | 223.5 nm | 0.54 → 0.48 | 89° |

The third is the resist line grating used as an ion-mill mask, imaged before
transfer: it holds its width over three times the depth of the other two. The
taper is therefore introduced by the transfer etch and is not a property of the
lithography, which is the sort of attribution a single measurement technique
cannot make.

The widths are read from a brightness threshold, so they carry the edge
brightness of the micrograph and the plane of the cleave with them, and the
sidewall angles are good to a degree or two rather than better. That is ample
for the conclusion drawn, which is that the angle is not 90°.

### What this changed

[`Trapezoidal`](../src/gratinglab/profiles.py) was added: a binary groove with a
sidewall angle and a corner radius, which reduces to `Lamellar` exactly when the
sidewall is vertical and the corners are sharp. `Lamellar` itself was left
alone, because it is the idealisation the published benchmark tables are
computed for and the validation rests on it being exactly that.

### What it does not establish

These are not the grating whose efficiency the ALS campaign measured. They come
from the same programme and have the same period, but no cross-section in the
collection is of OGRE laminar I itself, and nothing here licenses using these
numbers as that grating's geometry.
The entry is evidence that the shape to assume for an etched grating is a
tapered one, not a measurement of any particular grating in the corpus.

Nor does it close "The measured groove is rounded rather than faceted". That
finding names a cross-section as one of the three things that could settle
whether the rounding belongs to the grating or to the probe, and a cross-section
of *that* grating would indeed settle it. The collection holds none: its
cross-sections are of laminar gratings and date from 2021, and the TASTE scans
run from 2022. The two have never been taken of one sample, which is the
experiment that would close the question.

### A pixel size that is wrong by a factor of two

SmartSEM records two pixel sizes and the full-precision one in the tag's leading
block is the *display* pixel, not the image pixel. They agree at a store
resolution of 1024 × 768 and differ by exactly two at 2048 × 1536, which is the
resolution of the 2023 background scans. Taking the block as written doubles
every length measured from those images while leaving the micrograph looking
entirely correct — the same shape of trap as the `Sens. ZsensSens` key in the
Nanoscope reader, and caught the same way, by cross-checking one recorded
quantity against another. The reader anchors the scale on `AP_IMAGE_PIXEL_SIZE`
and takes the block's extra digits only after confirming the whole-number ratio
between them.

Reproduced by `tests/metrology/test_sem_reader.py`, which builds its own
SmartSEM files and needs none of the group's data. The micrographs themselves
are the group's and are not in this repository.

---

## The measured comparison ran, and the rigorous solver does not reproduce the ALS campaign

The comparison the roadmap deferred on 2026-10-24 for want of data has now
run. It does not validate the solvers, and the reason is not the solvers.

The beamline tree was found on a third machine (neither MC1 nor MAC1), at
`$ALS_DATA_ROOT`, together with the PCGrate corpus. `als_eff_meas` reduces
`ogre_laminar_I_oct21` to 251 absolute per-order efficiencies with
uncertainties, across orders −2 to 7 and 39 wavelengths from 1.0332 to
2.8178 nm (1200 down to 440 eV). `tests/test_measured.py::TestOgreLaminarCampaign`
passes for the first time anywhere, including its 251-point baseline, so the
`to_measured()` conversion written out in that fixture is confirmed against
real data.

**Two of the roadmap's open questions are closed in passing.**

1. *The alpha sign convention.* `als_eff_meas`'s recovered `alpha` is already
   in `Illumination.alpha_deg`'s convention, with no flip. At α = +27.22°
   every one of the 251 measured points propagates under gratinglab's grating
   equation; at α = −27.22°, 107 of them are evanescent. The two-candidate
   search in the fixture's `_oriented()` is therefore no longer a guess, and
   the sign is pinned by a test.
2. *Convergence at this geometry.* All three configurations plateau, against a
   tolerance of 1e-5 chosen from the measurement rather than from habit: the
   campaign's median σ is 4e-5, so 1e-5 is below the error bar of every point
   the model is compared against.

| configuration | knob | ladder | last two differences | converged at |
|---|---|---|---|---:|
| integral, `conductivity="tabulated"` | `boundary_points` | 128–2048 | 9.1e-6, 1.3e-6 | **512** |
| integral, `conductivity="perfect"` | `boundary_points` | 128–4096 | 5.5e-6, 8.1e-7 | **1024** |
| scalar | `quadrature_points` | 1024–1048576 | 1.5e-5, 7.7e-6, 3.9e-6 | **262144** |

The scalar ladder's differences halve exactly on each doubling, which is clean
first-order convergence and the reason it needs 2.6e5 points to reach a
tolerance the integral method reaches at 512 boundary nodes.

### The comparison cannot be a prediction, because the groove is not recorded

The campaign records the period (158.0 ± 0.5 nm, and the reduction re-derives
157.8 ± 2.7 nm from the measured dispersion, which agrees) and recovers the
angles (γ = 1.857° ± 0.025, α = 27.215° ± 0.487). It does not record the
groove depth or the duty cycle, and the original analysis never modelled the
groove at all: its legacy script plots the measurement against bare Au
reflectivity at 1.7°.

So the best available comparison is a fit, and it is reported as one. Depth was
scanned from 6 to 44 nm in 2 nm steps against duty cycle 0.30 to 0.93, 280
rigorous sweeps of 39 wavelengths each, with `conductivity="tabulated"` and
`coating="Au"` so that the efficiencies are absolute. Two free parameters
against 251 points is heavily over-constrained, and the order-to-order ratios
at fixed wavelength are what carry the constraint.

**The surface is shallow and multi-modal, and no groove fits.** Over the whole
grid the RMS residual moves only between 0.0706 and 0.1841. The best is
16 nm / 0.70 at 0.0706; a second, essentially equal basin sits at 34 nm / 0.60
at 0.0712. Duty cycles above 0.72 are monotonically worse, so 0.70 is a real
ridge and not a truncated search.

### Residuals at the best fit

Each solver at its own converged setting, 16 nm deep, duty 0.70:

| method | RMS | max abs | where | mean | χ² / 251 pts |
|---|---:|---:|---|---:|---:|
| integral, tabulated | 0.0706 | 0.2909 | order 0, 2.583 nm | −0.0148 | 1.36e9 |
| integral, perfect | 0.1313 | 0.5162 | order 0, 1.107 nm | −0.0703 | 3.74e9 |
| scalar | 0.0852 | 0.4142 | order 0, 1.127 nm | −0.0071 | 1.05e9 |

The χ² column is recorded only to be dismissed. Every one of the 251 points
carries a σ, so the number has its full denominator, and a reduced χ² of 5.4
million is not a goodness-of-fit: the campaign's uncertainties are
counting-statistics propagated through the reduction, with no term for the
systematics that dominate here. They bound the measurement's repeatability,
not its accuracy.

Per order, integral with tabulated conductivity:

| order | n | RMS | mean | ⟨measured⟩ | ⟨model⟩ |
|---:|---:|---:|---:|---:|---:|
| −2 | 14 | 0.0158 | −0.0137 | 0.0189 | 0.0326 |
| −1 | 37 | 0.0405 | −0.0023 | 0.0947 | 0.0970 |
| 0 | 39 | 0.1417 | −0.0504 | 0.2036 | 0.2540 |
| 1 | 39 | 0.0850 | −0.0017 | 0.1263 | 0.1279 |
| 2 | 39 | 0.0422 | −0.0147 | 0.0488 | 0.0635 |
| 3 | 34 | 0.0297 | −0.0018 | 0.0303 | 0.0321 |
| 4 | 25 | 0.0242 | −0.0148 | 0.0249 | 0.0398 |
| 5 | 15 | 0.0228 | −0.0197 | 0.0217 | 0.0413 |
| 6 | 7 | 0.0146 | −0.0142 | 0.0149 | 0.0291 |
| 7 | 2 | 0.0050 | −0.0047 | 0.0062 | 0.0109 |

Per wavelength, binned, with the ratio of the measured to the modelled sum over
the orders that wavelength measured:

| eV band | λ (nm) | pts | RMS | measured/model total |
|---|---|---:|---:|---:|
| 1120–1200 | 1.033–1.107 | 47 | 0.0518 | 0.911 |
| 1020–1100 | 1.127–1.216 | 42 | 0.0641 | 0.904 |
| 920–1000 | 1.240–1.348 | 39 | 0.0516 | 0.896 |
| 820–900 | 1.378–1.512 | 30 | 0.0619 | 0.906 |
| 720–800 | 1.550–1.722 | 30 | 0.0784 | 0.815 |
| 620–700 | 1.771–2.000 | 25 | 0.0606 | 0.760 |
| 520–600 | 2.066–2.384 | 24 | 0.0865 | 0.758 |
| 440–500 | 2.480–2.818 | 14 | 0.1412 | 0.829 |

That ratio is a genuine energy statement and not a partial sum: at this
geometry the campaign measured **every order the model says propagates**, at
all 39 wavelengths. The real grating diffracts 15% less than the model does,
averaging 0.848 ± 0.071 and reaching 0.70 at 620 eV.

### What dominates the residuals

**Zero order, and then order 1.** Order 0 carries RMS 0.1417 against a mean
measured value of 0.2036, and holds the largest single residual at every
method's maximum. Orders 0 and 1 together are 85% of the total squared
residual while being 31% of the points.

**Not the recorded inputs.** Each input uncertainty was propagated through the
converged solver at the best-fit groove, and none of them is large enough to
matter:

| perturbation | RMS | max shift in the model |
|---|---:|---:|
| baseline | 0.07064 | — |
| γ ∓ 0.025° | 0.07112 / 0.07044 | 0.0170 |
| α ∓ 0.487° | 0.07069 / 0.07065 | 0.0090 |
| polarization TE / TM | 0.07056 / 0.07071 | 0.0009 |
| period ∓ 0.5 nm | 0.07066 / 0.07061 | 0.0022 |

The unrecorded polarization, which the roadmap flagged as something the
comparison would have to settle, turns out not to matter at all: TE and TM
differ by at most 9.1e-4 across all 251 points at this graze, which is
consistent with the s/p degeneracy already recorded in this file and makes the
fixture's `"unpolarized"` default harmless here. α's ±0.487° looks alarming and
is not: it moves the model by 0.009, an eighth of the residual.

**Not the absolute level alone.** The obvious reading of a 15% energy deficit is
that the diffraction is right and the reflectivity is wrong, which is the
division this file has drawn before under "the corpus can test the diffraction
but not the reflectivity". That reading is wrong here. Normalising both
measurement and model at each wavelength by their own sum over the measured
orders leaves an RMS of 0.1200 in normalised units, against 0.0706/0.6 ≈ 0.118
for the raw residual carried into the same units. Removing the level removes
none of the disagreement. The pattern itself disagrees.

**The pattern disagreement alternates with the parity of the order.** In
normalised fractions, measured minus modelled runs −0.019, +0.024, −0.020,
+0.034, −0.008, +0.006 for orders −2 through 3: every even order over-modelled,
every odd order under-modelled. For a lamellar groove that is the duty-cycle
signature — the even orders vanish at duty 0.5 and grow as the duty leaves it —
so the data wants a duty closer to 0.5 than 0.70, while the envelope and the
absolute level want 0.70. No single rectangle satisfies both, which is why the
RMS surface is shallow in the first place.

**Not sidewall slope.** A real etched groove has finite sidewalls, and the
over-modelling of orders 4 to 6 by 40–50% is the signature of corners sharper
than the real ones. `FromProfileData` carries a trapezoid, so this is testable,
and it fails as an explanation. Sweeping the sidewall run at the best fit:

| sidewall run | RMS, all orders | RMS, orders ≥ 4 | mean, orders ≥ 4 |
|---:|---:|---:|---:|
| 0 nm | 0.0706 | 0.0222 | −0.0158 |
| 8 nm | 0.0728 | 0.0247 | −0.0212 |
| 18 nm | 0.0901 | 0.0192 | −0.0152 |
| 25 nm | 0.1047 | 0.0096 | −0.0007 |

A 25 nm sidewall fixes the high orders almost exactly, cutting their RMS by 2.3
times and their bias to nothing, and in doing so pushes the overall RMS up by
half. The high orders and the low orders make incompatible demands, in this
profile family as in the rectangular one.

**Not the solver's own consistency.** At the same settings the finite-conductivity
energy theorem holds to 1.4e-6 (R + A = 1 at `boundary_points` 1024; 7.0e-6 at
512), four to five orders of magnitude tighter than the disagreement with the
measurement. Worth noting in passing: `check_energy_balance` returns FAIL on
this scan, because 16 of the 39 points exceed unity — by 7e-6. The report is
correct and the tolerance is doing its job, but a reader seeing FAIL should
read the magnitude before concluding anything.

**Scalar fails differently, and as predicted.** Its RMS (0.0852) is not much
worse than the rigorous solver's, which flatters it. Per order it is
structurally wrong: it under-models order 3 by seventeen times (0.0018 against
0.0303 measured) and orders 4 to 6 by an order of magnitude, while
over-modelling orders −2 and −1. This is the regime its own guard already
refuses — λ/(p·sin γ) runs 0.20 at 1.033 nm to 0.55 at 2.818 nm, against the
0.32 this file already records as the point where scalar theory departs.
Perfect conductivity is worse still (0.1313) and for a plainer reason: it
returns efficiency relative to a perfect reflector, so comparing it to absolute
measured data charges it the whole 35–39% the model absorbs. It is kept in the
table as a diagnostic of the diffraction pattern, not as a candidate.

### The verdict on 0.3, from this campaign

*Superseded on 2026-10-29 by "The verdict on 0.3, from the TASTE campaign" at
the end of this file, which met the criterion. What follows is why this campaign
could not, and it is still the reason: the absence of a groove. Nothing in it
was wrong.*

**The criterion is not met.** The roadmap's gate for the solvers to lose their
experimental status is "use on a real problem, with the integral solver
compared against measured efficiencies ... and not only against published
benchmarks and its own physics checks". The comparison has now been performed,
which satisfies the sentence read literally, and that is the weaker half of
what the gate was for. The gate exists to establish that the solver is right
about a real grating, and this comparison does not establish it.

What stands in the way is not the solver and not the measurement, but the
absence of the grating's geometry. The best fit over a two-parameter family
lands at an RMS of 0.0706, which is 18% of the largest measured efficiency in
the campaign, 1700 times the median σ, and four orders of magnitude above the
solver's own numerical self-consistency. The disagreement survives every
recorded uncertainty, both profile families tried, and the removal of the
absolute level. The weight of this project's other evidence — orderwise
agreement with an independent commercial integral code, Goray & Schmidt's
Table 3 to 6e-5, two-sided energy theorems, reciprocity — says the solver is
unlikely to be the error, but this comparison cannot say so, because a fitted
profile can absorb a solver error and a solver can be blamed for a profile
error, and nothing here separates them.

This is the same wall the TASTE run hit, recorded above under "the measured
groove is rounded rather than faceted": an absolute efficiency computed from an
assumed groove is not defensible, and the project has now met that conclusion
twice from opposite directions, once with a measured profile and no reference
and once with a reference and no profile.

What would close it is metrology on this grating, not more solving: an AFM scan
of OGRE laminar I, fed through `FromProfileData` with the tip correction the
metrology side already has, turning the fit into a prediction. Failing that, a
campaign that records its groove. Until one of those exists the solvers stay
experimental, and the honest statement of their status is unchanged by this
work: validated against another implementation of the same method and against
published tables, not yet against a measured grating.

The second of those two routes is the one that opened: the December 2018 TASTE
campaign both records its groove and has an AFM scan of the grating, and the
comparison against it is at the end of this file. The AFM scan of OGRE laminar I
is still worth taking, because it would make *this* campaign's 251 points a
prediction too.

Reproduced by `tests/test_measured.py::TestOgreLaminarCampaign`, which skips
without `als_eff_meas` and `$ALS_DATA_ROOT`. The measurement itself is not in
this repository and must not be committed.

---

## The TASTE campaign: a measured groove and a measured efficiency, at last

The gate under "The measured comparison ran" closed on the absence of a
groove: the OGRE campaign recorded none, so the comparison could only be a
two-parameter fit, and a fit cannot separate a solver error from a profile
error. The December 2018 TASTE campaign is the case that absence pointed at.
It has both halves of the pair — 296 absolute per-order efficiencies and an
AFM scan of the grating they were measured on — and it was published, as
McCoy, McEntaffer & Miles, *ApJ* **891**, 114 (2020), erratum *ApJ* **926**,
105 (2022), with its mount in a table and its own model for comparison.

This entry runs: what the reduction had to recover before the data could be
used at all, a convention error the data found, what the metrology says the
groove is, the comparison and where its residual lives, the paper's own
suspicion about the apex tested as a labelled hypothesis, two defects the
comparison turned up in the solvers themselves, and the verdict on 0.3.

### The published reduction is not the surviving code

`als_eff_meas` gained two campaigns, `taste_dec18_round3a` (440–800 eV) and
`taste_dec18_round3b` (80–420 eV). They are two and not one because each half
has its own arc fit: the low-energy half was taken with the beamline's
triple-mirror order sorter inserted, which moves the beam on the grating. The
paper allows for this and calls the effect small — "the effect is small and not
apparent in the measured absolute efficiency data" — and on γ it is not small.
Round 3a's fit is the published Table 1 to the digit (γ = 1.7304° against
1.73 ± 0.03, α = 23.3925° against 23.4 ± 0.6, and six more); round 3b's is
γ = 2.2125°, sixteen published sigmas away.

Both halves reproduce their saved `results/` to better than 5 × 10⁻⁹, which is
exact at the eight decimals those files carry. Getting there meant discovering
that **the driver script on disk is not the record of the published run**, three
times over:

| what the surviving code says | what the saved output needed |
|---|---|
| `als_eff_3a.py`: thirty `#delete` blocks | twenty-two points excluded |
| `als_eff_3b.py`: six `#delete` blocks | none of them applied |
| `december_2018_taste3b.py`: `R_correct = 1.05 # 0.85` | 0.825, and 1.05 cannot complete the reduction |
| `als.py`: `throw = ufloat(235.0, 2.0)`; Table 1: 232.0 ± 1.4 mm | 233.0 mm |

For OGRE the driver and the saved output agreed, and
`als_eff_meas/docs/VERIFICATION.md` warned only that the `#delete` comments
disagreed with the `#delete` code. Here the whole script post-dates the run.
`results/` is byte-identical to `results_3a_rev23`, the last of that round's
twenty-one saved revisions and the state the paper was submitted from a week
later, so the saved output is the record and the script is a later draft.

Two of those four are worth their own line.

**`R_correct` is not cosmetic.** It scales the fitted arc radius between the
arc fit and the order prediction, so it decides how many orders are counted as
propagating. The original's per-energy noise floors are hand-written lists with
one entry per order, so a wrong value does not give wrong numbers — it raises a
broadcast error. Those lists are an accidental and strict check on it, and they
reject the 1.05 that the code actually sets.

**The throw is nobody's published number.** 232.0 ± 1.4 mm admits a tenth order
at 440–500 eV against the nine the noise lists are written for there; 235.0
drops a fifteenth at 740–800 eV. 233.0 reproduces all eight of round 3a's geometry
numbers and all 296 efficiencies, and sits inside Table 1's uncertainty. Round
3a's efficiencies are insensitive to it anywhere in 232.1–233.9; round 3b's are
not, and pick out 233.0, so one value serves the whole beamtime.

One thing could not be recovered. Round 3b's module defines no direct-beam
block, and the `dir_x` its saved α and yaw imply, −0.183 ± 0.005, matches no
scan surviving in that round: the direct beam sits at −0.161 mm at 80 eV,
drifts to +0.088 at 180 eV where the gain and the order sorter changed, and
ends at −0.008 at 420 eV. Round 3a's holds −0.057 to −0.064 across its whole
sweep. The campaign uses round 3b's own first direct-beam pair instead, which
gives α = 19.24° against the saved 19.4 and yaw 0.729° against 0.730 — close,
and still not it. It costs nothing, because `direct_x` feeds only α and yaw and
neither enters the efficiency reduction, but the α this campaign reports is not
the mount. Table 1 is the mount.

**Two measured points should not exist.** At Table 1's mount the grating
equation makes order −4 at 680 eV and order −2 at 340 eV evanescent, and the
campaign has a measured point for each. Both are the most negative order their
energy recorded, at the edge of the arc, and order −2 at 340 eV is precisely
the one point of round 3b whose reduction is sensitive to `R_correct`. The model
is not wrong to call them evanescent and the reduction is not wrong to have
found something there, but one of the two is mistaken.
`tests/test_measured.py::TestTasteCampaign` names them, so their number cannot
grow quietly.

### The reduction reproduces the paper's prose, not only its files

Agreement with `results/` proves the port, not the physics. Three statements the
paper makes in words, which no saved file carries, come out of the reduction
independently:

- **"Peak-order efficiency ranges from about 75% down to 25%."** The strongest
  non-zero order runs 0.741 at 160 eV down to 0.157 at 760 eV, and the order
  carrying it walks +1 → +7 as the energy rises, which is the paper's orders 0
  through 7.
- **"The sharp cutoff in the measured n = 2 curve at 160 eV."** Order 2 appears
  at 160 eV and not below, in the reduction as in Figure 13.
- **"Total relative diffraction efficiency ranges from about 95% to 88% as
  photon energy increases from 240 to 800 eV."** Summing every measured order
  and dividing by the Fresnel reflectivity of gold at ζ = 1.73° gives 0.967 at
  240 eV falling to 0.882 at 800 eV.

The third is worth more than a cross-check. Overlaid on it, the Névot–Croce
factor for σ = 1.5 nm rms — the paper's Figure 15 argument — tracks it to
within 0.013 across the whole band, and to within 0.002 below 500 eV:

| eV | measured total / Fresnel Au | Névot–Croce, σ = 1.5 nm |
|---:|---:|---:|
| 240 | 0.967 | 0.971 |
| 320 | 0.952 | 0.951 |
| 400 | 0.938 | 0.937 |
| 480 | 0.925 | 0.926 |
| 560 | 0.909 | 0.916 |
| 640 | 0.895 | 0.907 |
| 720 | 0.892 | 0.901 |
| 800 | 0.882 | 0.895 |

So **the absolute level of this measurement is accounted for to about 1%** by
gold's reflectivity and a 1.5 nm rms roughness, with nothing structured left
over. That is a stronger starting point than the OGRE campaign gave: there the
model diffracted 15% more than the grating did, and normalising the level away
removed none of the disagreement. Here the level is explained before the
comparison begins, so what the comparison below finds is almost entirely a
disagreement about the *distribution among orders*.

One consequence for the solvers. `Problem.roughness` is consumed by the scalar
solver alone; the three rigorous backends solve the smooth boundary and ignore
it. The paper applies its roughness factor outside PCGrate for the same reason,
so the comparison below is run both ways: the rigorous model as the solver
returns it, and the same model multiplied by the Névot–Croce factor above, which
is the paper's procedure.

### The erratum does not touch any of this

The erratum (*ApJ* **926**, 105) corrects three things: a missing space in
Equation 1, a swap of the φ and ϕ glyphs in the online version, and a
clarification about why the coated bottom plateau looks widened (non-conformal
EBPVD rather than surface composition). It also replaces Equation 14 with
φ = arcsin(Δx₀/r) − α, noting that the published form held only for small roll
and that Table 1's value was derived with the corrected one.

`als_eff_meas.geometry.arc_angles` computes `roll = arcsin(delta_x_zero / R) -
alpha`, which is the erratum's form, so the port carries the correction without
having been told about it: the reduction's 1.1744° reproduces the saved −1.17
and sits inside the published 1.14 ± 0.04. Roll enters nothing this comparison
uses, but it is worth knowing that the number the port reproduces is the
corrected one and not the printed equation's.

That same function carries a known grouping quirk — `eta` is computed as
`(|Δy|/L)·cos(roll)` and not `|Δy|/(L·cos(roll))`, preserved because it is what
produced the published angles. TASTE's roll is twice OGRE's, so it is worth
restating the size: the two forms differ by cos²(roll) inside the arcsin, which
here is 4 parts in 10⁴, or 0.0007° of η against a quoted ±0.04°. Still
immaterial, now at twice the roll.

### A measured groove has no handedness of its own

Before any of the numbers below meant anything, a convention error had to be
found, and it was found by the data rather than by reading.

An AFM scan fixes a shape, not an orientation. Nothing in a scan records which
way the sample faced the beam, or the sign of the fast scan axis against the
dispersion direction. The convention does fix it: a `.ggp` and the metrology
pipeline's `x_norm` run along `+d̂`, `t` runs against it, and
`FromProfileData.along_dispersion` is the one place the two meet
(`docs/conventions.md`, "Direction of the profile parameter"). Building
`FromProfileData(t=..., y=...)` straight from a `.ggp` skips that handover and
solves the groove as its mirror image.

The failure is not loud. At 440 eV, where the campaign's blaze order is n = +4,
the measured groove through `along_dispersion` puts 0.437 into +4 and 0.092
into all negative orders together; mirrored, +4 keeps 0.076 and order −2 alone
grows to 0.224. Both look like a blazed grating. One of them blazes where the
grating does not, and compared against measured data it reads as a profile
error rather than as a convention error — the same trap as "The integral solver
solved every blazed grating as its mirror image", one layer further out, in the
caller instead of the solver.

`tests/test_measured.py::TestTheMeasuredGroovesHandedness` pins it, needs
neither the measurement nor an optical-constant table, and so always runs.

### What the AFM says the groove is

Three scans of this grating survive, all at 400 nm period and 2 µm scan width:
`TASTE_post-test.0_0000{0,1}` from 2018-12-18, ten days after the beamtime, and
`bareTASTE_ALS.0_00000` from 2022-09-01, a re-scan of the same sample at
512 rows rather than 128. The mould before coating is there too
(`final_TASTE.0_00000`, 2018-11-16, peak-to-valley 130 nm, which is the 130 nm
PMMA the paper names); the coated grating measures 10 nm shallower, exactly as
the paper reports.

Through `metrology.boundary` with tip correction at the probe's nominal
condition (Bruker ScanAsyst-Air, R = 2 nm, 17.5° sides):

| scan | rows | grooves | period | depth / period | groove residual | pixels certain |
|---|---:|---:|---:|---:|---:|---:|
| post-test 2018, scan 0 | 128 | 4 | 396.2 nm | 0.2933 | 3.68 nm RMS | 96.7% |
| post-test 2018, scan 1 | 128 | 4 | 412.3 nm | 0.2880 | 11.46 nm RMS | 96.1% |
| bare 2022 | 512 | 3 | 408.8 nm | 0.2847 | 1.79 nm RMS | 96.3% |

The groove residual is the scatter of the individual grooves about their
average, and it is what the 2022 re-scan bought: 1.79 nm against 3.7 and
11.5 nm, from four times the rows. Tip correction changes the depth of none of
them — 0.2933 to 0.2930 for the worst probe within Bruker's specification,
0.2933 to 0.2933 for the nominal one — which is the same null result recorded
under "A nominal tip does not explain the rounded groove", now on a second
grating: the measured surface is already reachable by any defensible tip, so
erosion has no hidden depth to return.

Swath statistics on the same scans, which re-measure the same grooves in bands
of rows and so describe variation rather than independent samples:

| scan | swaths | depth / period | period | shape spread between swaths |
|---|---:|---:|---:|---:|
| bare 2022 | 16 | 0.2854 ± 0.0023 | 407.3 ± 5.1 nm | 1.07 nm RMS |
| post-test 2018 | 8 | 0.2882 ± 0.0037 | 398.5 ± 8.3 nm | 1.93 nm RMS |

Three things follow, and one does not.

**The period is 400 nm.** The AFM's own estimate runs 398.5 ± 8.3 to
407.3 ± 5.1 nm depending on the scan, and the reduction's re-derivation from the
measured dispersion gives 401 ± 8 nm in round 3a and 382 ± 14 in round 3b.
Three of those four cover the lithographic 400 nm within their own scatter; the
fourth, round 3b's 382 ± 14, falls 1.3σ short of it and is the noisier of the
two dispersion estimates. So the model takes 400 nm and the AFM is used for
shape only — and that 382 returns, uncomfortably, under "One number explains the
residual" below. The alternative, trusting the AFM laterally and rescaling the
depth to the 400 nm period, moves depth/period from 0.2847 to 0.2910, and its
effect on the comparison is in the budget below.

**The three scans agree on the shape.** Registered to their best lateral shift
and compared point for point, the three normalised grooves sit 3.0 to 4.4 nm
RMS apart, against a depth of 114 nm — and the 2022 scan with and without tip
correction differ by 0.32 nm. Four years and a different session change the
measured groove by about 3% of its depth.

**The apex is flat and the floor is flat.** The measured groove is a trapezoid:
a bottom land of 93–95 nm, a steep side rising at 66–67°, a top land of
68–70 nm, and a facet descending at 27.9° ± 0.2 (the fitted blaze angle, from
60 grooves across 20 row groups of the 2022 scan; the 2018 scans fit 29.4 ± 0.2
and 28.5 ± 0.1). The paper's own fitted trapezoid, from PCGrate and never
compared against this scan, is a flat top of 77 nm and a bottom of 85 nm at
blaze 27° and depth 120 nm. The agreement is close enough to be worth saying
plainly: **the AFM independently confirms the groove the paper fitted to the
efficiency data**, to about 10 nm on each land and 6 nm in depth.

What does not follow is the steep side. The measured 66–67° is where the
probe's own flanks sit — 65° for the 25° back face, 72.5° for the 17.5° sides —
so the AFM cannot distinguish that wall from the near-vertical 88.8° the
paper's trapezoid needs, and the measurement is consistent with either. This is
the same tip limit recorded under "A nominal tip does not explain the rounded
groove", on the other facet.

### Convergence, and one solver that does not reach the tolerance

The tolerance is taken from the measurement rather than from habit: the
campaign's median σ is 3.0 × 10⁻⁵, so 10⁻⁵ sits below the error bar of every
point the model is compared against. Three probe wavelengths spanning the band
— 1.5498, 3.0996 and 15.4980 nm, i.e. 800, 400 and 80 eV — with `plateau` at its
default of two, so a pass needs three knob values in a row that agree.

| configuration | knob | ladder | last two differences | converged at |
|---|---|---|---|---:|
| C-method, `conductivity="tabulated"` | `truncation` | 32–128 | 6.9e-6, 5.3e-6 | **32** |
| integral, `conductivity="perfect"` | `boundary_points` | 256–1024 | 3.1e-6, 2.0e-6 | **256** |
| scalar | `quadrature_points` | 65536–1048576 | 2.5e-6, 3.1e-6 | **65536** |
| integral, `conductivity="tabulated"` | `boundary_points` | 512–2048 | 3.4e-4, 4.9e-5 | **not reached** |
| RCWA | `truncation` | 20–80 | 3.0e-3, 9.7e-4 | **not reached** |

The two failures are different in kind.

**RCWA is in the wrong regime and says so by not plateauing.** Its differences
fall by about a factor of three per doubling, so 10⁻⁵ is four or five more
doublings away, at a truncation near 10³ — an eigenproblem of order 4000 per
slab, with the slab count tied to the truncation. The reason is the aspect
ratio: this groove is 114 nm deep at wavelengths down to 1.55 nm, i.e. 74
wavelengths deep, staircased into slabs. RCWA's answer at truncation 80 is
nevertheless good, as the cross-method table below shows; what it cannot do at
this geometry is *demonstrate* that it is.

**The integral solver's tabulated mode converges, slowly, and is accurate well
before it gets there.** Its differences run 3.4 × 10⁻⁴ then 4.9 × 10⁻⁵ at 2048
nodes, which is a clean factor of seven per doubling and puts the tolerance one
more doubling out, at 4096. That rung was not run: at 2048 one wavelength
already costs over a minute, and the point is already made by the cross-method
comparison, where the integral solver at 1024 nodes agrees with the *converged*
C-method to 2 × 10⁻⁵ RMS. So `converged=False` here is a statement about the
plateau test and not about the answer, and this file's own rule applies —
"`converged=False` is itself a result and does not indicate a missing one".

A reader deciding what to run at this geometry should note the cost, which is
the opposite way round from the OGRE case: the C-method reaches the tolerance at
truncation 32 and sweeps all 37 wavelengths in 24 seconds, where the integral
method needs 1024 nodes and 50 seconds to get within 2 × 10⁻⁵ of it.

### The comparison

The model is a prediction, not a fit. Every input is measured: the period from
the lithography and confirmed by two dispersion re-derivations, the mount from
the paper's Table 1, the groove from the AFM scan through
`metrology.boundary`, the optical constants from a CXRO export, the roughness
from the same AFM. Nothing was adjusted to improve the agreement.

One practical note on the optical constants. The table this repository carries
for the test suite starts at 200 eV, and this campaign starts at 80, so an
absolute model across the whole band is impossible with it — `lookup` refuses
rather than extrapolate, which is correct. The runs below use
`Au_CXRO.dat` from the beamline tree, a CXRO export in the same format covering
80–1300 eV, supplied through `GRATINGLAB_MATERIALS` as `materials/SOURCES.md`
describes. That is the table the paper itself used, so the model and the paper's
own model see identical n(λ). Nothing was added to this repository.

| model | n | RMS | max abs | where | mean |
|---|---:|---:|---:|---|---:|
| C-method, tabulated, truncation 64 | 296 | 0.0383 | 0.1079 | order +5, 520 eV | −0.0058 |
| integral, tabulated, 1024 nodes | 296 | 0.0383 | 0.1079 | order +5, 520 eV | −0.0058 |
| RCWA, truncation 80 | 296 | 0.0383 | 0.1080 | order +5, 520 eV | −0.0057 |
| scalar, 65536 points | 296 | 0.0767 | 0.4932 | order 0, 80 eV | +0.0257 |

**The three rigorous methods agree with each other two hundred times more
closely than any of them agrees with the measurement.**

| pair | RMS difference | max difference |
|---|---:|---:|
| C-method vs integral | 2.0e-5 | 2.4e-4 |
| C-method vs RCWA | 2.0e-4 | 7.3e-4 |
| integral vs RCWA | 2.0e-4 | 7.3e-4 |

This is the thing the OGRE comparison could not do. Three independent
formulations — an integral equation on the boundary, a coordinate
transformation, and a Fourier-space eigenmode expansion through a staircase —
land on the same 296 numbers to within 7 × 10⁻⁴, so the 0.038 that separates all
of them from the grating is not the solver. It is the geometry or the
measurement, and those are now the only candidates left.

### Where the residual is

Per order, C-method at its converged setting:

| order | n | RMS | mean | ⟨measured⟩ | ⟨model⟩ |
|---:|---:|---:|---:|---:|---:|
| −4 | 7 | 0.0142 | −0.0061 | 0.0156 | 0.0217 |
| −3 | 15 | 0.0161 | −0.0062 | 0.0207 | 0.0269 |
| −2 | 24 | 0.0127 | +0.0020 | 0.0263 | 0.0243 |
| −1 | 32 | 0.0172 | −0.0013 | 0.0273 | 0.0286 |
| 0 | 37 | 0.0223 | −0.0047 | 0.0879 | 0.0926 |
| 1 | 37 | 0.0361 | −0.0064 | 0.1458 | 0.1522 |
| 2 | 33 | 0.0360 | −0.0062 | 0.1112 | 0.1174 |
| 3 | 29 | 0.0400 | −0.0020 | 0.1021 | 0.1040 |
| 4 | 26 | 0.0505 | +0.0002 | 0.0975 | 0.0973 |
| 5 | 22 | 0.0579 | −0.0055 | 0.0973 | 0.1028 |
| 6 | 16 | 0.0665 | −0.0058 | 0.1125 | 0.1183 |
| 7 | 11 | 0.0601 | −0.0463 | 0.0778 | 0.1241 |
| 8 | 7 | 0.0353 | −0.0269 | 0.0150 | 0.0419 |

Three features, and none of them is the OGRE pattern.

**It is not the zero order.** Order 0 carries RMS 0.0223 against a mean measured
0.0879 — the *best* of the seven orders that carry real power. In the OGRE
comparison order 0 held RMS 0.1417 and the single largest residual of every
method, and orders 0 and 1 together were 85% of the total squared residual on
31% of the points. Here order 0 is modelled better than every order above it.

**There is no parity alternation.** The per-order means run −0.006, −0.006,
+0.002, −0.001, −0.005, −0.006, −0.006, −0.002, +0.000, −0.006, −0.006 for
orders −4 through 6: all small, all the same sign but one, with no even/odd
structure. The duty-cycle signature that made the OGRE fit self-contradictory is
simply absent, which it should be — this groove is not a rectangle and has no
duty cycle.

**The residual grows with order number.** RMS climbs almost monotonically from
0.013 at order −2 to 0.067 at order 6 — orders 1 and 2 are the only inversion,
and they differ by 1 × 10⁻⁴ — and orders
7 and 8 are over-modelled by 0.046 and 0.027 against measured means of 0.078 and
0.015 — that is, the model puts 60% and 180% more into them than the grating
does. The paper saw the same thing and said so: "measured peak-order efficiency
becomes increasingly diminished relative to the model as order number increases
beyond n = 4". It is the high orders, the ones most sensitive to the fine
structure of the groove, that the model gets wrong — which is where a groove
measured on three periods with a 1.8 nm residual would be expected to fail.

Per energy band, with the ratio of the measured to the modelled sum over the
orders that band measured:

| eV band | n | RMS | measured/modelled total |
|---|---:|---:|---:|
| 80–160 | 11 | 0.0433 | 0.975 |
| 180–260 | 22 | 0.0476 | 0.968 |
| 280–360 | 31 | 0.0326 | 0.958 |
| 380–460 | 40 | 0.0399 | 0.938 |
| 480–560 | 47 | 0.0430 | 0.916 |
| 580–660 | 54 | 0.0392 | 0.900 |
| 680–800 | 91 | 0.0326 | 0.887 |

That ratio is a genuine energy statement: at every one of the 37 wavelengths the
campaign measured every order the model says propagates, bar the two evanescent
points named above. It falls monotonically from 0.975 to 0.887, and it is the
only part of this comparison that is *not* a disagreement about shape — it is
the roughness the rigorous solvers do not model.

### The roughness, put back the way the paper puts it

`Problem.roughness` reaches the scalar solver only; the three rigorous backends
solve the smooth boundary. The paper multiplies its own PCGrate output by the
Névot–Croce factor for σ = 1.5 nm rms at ζ = 1.73°, and doing the same here:

| model | RMS | max abs | mean | 680–800 eV total ratio |
|---|---:|---:|---:|---:|
| rigorous, smooth boundary | 0.0383 | 0.1079 | −0.0058 | 0.887 |
| × Névot–Croce, σ = 1.5 nm | 0.0363 | 0.1136 | −0.0007 | 0.987 |

The factor is not fitted — σ = 1.5 nm is the AFM's number and ζ = 1.73° is
Table 1's — and it takes the mean residual from −0.0058 to −0.0007 and the
energy ratio from 0.887 to 0.987. **The band-wide level discrepancy is
roughness, to within 1.3%, and nothing else.** What it does not do is reduce the
scatter much (0.0383 to 0.0363) or the maximum (which grows), because the
factor is smooth in wavelength and the residual is not: it lives in the
distribution among orders.

This is a cleaner division than the OGRE comparison could draw. There,
normalising the level away left all of the disagreement behind. Here the level
is explained by an independently measured roughness and the remainder is
entirely orderwise.

### What the recorded uncertainties do

Each input perturbed by its own recorded uncertainty, through the converged
C-method, against the measured groove:

| perturbation | RMS | max shift in the model |
|---|---:|---:|
| baseline | 0.03827 | — |
| γ − 0.03° | **0.02754** | 0.0410 |
| γ + 0.03° | 0.04901 | 0.0399 |
| α − 0.6° | 0.03839 | 0.0259 |
| α + 0.6° | 0.03807 | 0.0171 |
| period ∓ 0.5 nm | 0.03747 / 0.03908 | 0.0028 |
| polarization TM | 0.03845 | 0.0068 |
| groove depth ∓ 1 swath σ (0.0023) | 0.03448 / 0.04211 | 0.0127 |
| AFM lateral scale: depth × 1.022 | 0.04876 | 0.0347 |
| roughness 0 or 3 nm | 0.03827 | 0.0000 |
| fitted blaze angle ∓ 0.23° | 0.03827 | 0.0000 |

**γ's published ±0.03° moves the model by as much as the whole residual.** This
is the opposite of the OGRE finding, where no recorded uncertainty came within
an eighth of the disagreement. At γ = 1.70°, one published sigma low, the RMS
drops to 0.0275; at 1.76° it rises to 0.049. So the comparison does not
distinguish "the model is right and γ is 1.70" from "γ is 1.73 and 0.038 of the
residual is real", and the data cannot be made to: γ enters through
`sin γ`, which sets the whole order spacing, and a 1.7% change in it is a 1.7%
change in every diffracted angle.

Two entries are zero by construction and are listed so that no reader mistakes
their absence for a null result. `Problem.roughness` does not reach a rigorous
solver, as above. `FromProfileData.blaze_angle` reaches only the scalar solver's
facet reflectivity; a rigorous solver takes the boundary itself and never
consults the fitted facet angle, which is why ±0.23° on a fit from 60 grooves
changes nothing here.

The polarization, which the campaign did not record, again does not matter: TE
and TM differ by at most 6.8 × 10⁻³ across all 296 points, consistent with the
s/p degeneracy this file already records for these grazes.

### The apex: the paper's suspicion, tested as a hypothesis

The paper offers two accounts of its own groove. Its model is an acute trapezoid
with a *flat* apex, 77 nm wide, which it attributes to the 100 nm top staircase
step that the GEBL pattern leaves unexposed and the reflow therefore barely
touches. Its discussion then suspects something else — that the apex is
"slightly rounded as a result of the thermal reflow process, but this is
difficult to verify through AFM because the shape of the SCANASYST-AIR tip is
convolved with the true grating topography". Flat and rounded are different
shapes, and the measurement can tell them apart.

Everything below is the integral solver, because the C-method cannot do it: see
the capability note after the table. All at the recorded mount, all at the
paper's facet angle of 27°, and all at its depth of 120 nm except the two rows
that say otherwise. The first four rows hold the total land at the paper's
77 + 85 = 162 nm and only move it between the floor and the apex, so they are a
one-parameter family; the last two change the depth and the bottom land
respectively, and are labelled with what they change. **None of these is the
main result and none of them is fitted into it.**

| groove | depth (nm) | steep side | RMS | max abs |
|---|---:|---:|---:|---:|
| sharp apex, the land all at the floor | 119.9 | 88.8° | 0.0481 | 0.1695 |
| flat apex 40 nm | 120.0 | 88.8° | 0.0460 | 0.1777 |
| **flat apex 77 nm — the paper's own fit** | 120.0 | 88.8° | **0.0345** | 0.1573 |
| flat apex 120 nm | 120.0 | 88.8° | 0.0335 | 0.1093 |
| sharp apex, reflowed over 20 nm | 115.0 | 88.8° | 0.0446 | 0.1375 |
| sharp apex, reflowed over 40 nm | 109.9 | 88.8° | 0.0418 | 0.1302 |
| sharp apex, reflowed over 80 nm | 99.8 | 88.8° | 0.0456 | 0.1567 |
| sharp apex, reflowed over 160 nm | 79.7 | 88.8° | 0.0930 | 0.3279 |
| the paper's fit at the measured depth, 114 nm | 113.9 | 82.8° | 0.0284 | 0.1341 |
| the paper's fit with the bottom land cut to 37 nm, giving the measured 67° steep side | 120.0 | 67.1° | 0.0288 | 0.1554 |
| **the measured groove, bare 2022** | 113.9 | 67.1° | **0.0383** | 0.1079 |
| the measured groove, post-test 2018 | 117.3 | 66.4° | 0.0631 | 0.1928 |

"Reflowed over w" is the sharp trapezoid smoothed by a periodic boxcar average
of width w, which is the simplest model of a thermal reflow: it rounds every
convex corner, fills every concave one, conserves the enclosed area, is monotone
in w, and needs no tangency algebra where a 27° facet meets a near-vertical
wall. It costs depth, as rounding must, and it rounds the floor as well as the
apex, which reflow also does.

**An unpointed apex is required, and the data prefer it flat.** Going from a
sharp apex to the paper's 77 nm land improves the RMS by 28%, from 0.0481 to
0.0345, and widening the land to 120 nm improves it a little further. Rounding
the same sharp apex does less: 40 nm of reflow gets to 0.0418 and more reflow
then makes it worse, because rounding takes depth away as it goes while a land
does not. So the paper's *model* is supported and the paper's *suspicion* is
not the better of the two accounts it offered. On this measurement the apex is
flat-topped rather than rounded, which is what the mechanism the paper names —
an unexposed top step that the reflow does not reach — actually predicts.

**The AFM groove is not the best of these, and that is informative.** The
measured groove gives 0.0383, which the paper's trapezoid beats when its depth
is brought to the measured 114 nm (0.0284) or its bottom land cut until the
steep side matches the measured 67° (0.0288). Each of those is a one-parameter
nudge within the AFM's own uncertainty, so neither is a rival measurement; both
are evidence that the measured groove is not yet pinned down to the level the
comparison can feel. The measured groove's
defence is that it is a measurement rather than a three-parameter shape, and
its weakness is visible in this file already: it is an average over three
periods of a 2 µm scan, and the second scan of the same grating gives 0.0631.
The honest reading is that these four numbers — 0.028, 0.029, 0.038, 0.063 —
bracket what is knowable about this groove, and that no member of the family
does better than 0.028.

**The C-method's guard against a vertical facet is exact, and 88.8° is not
exact.** The flat-apex trapezoids have a near-vertical anti-blaze side, and the
C-method does not solve them: it returns RMS between 1.6 and 3.3 with individual
efficiencies up to 31, against the integral solver's 0.033 to 0.048 on the same
four grooves. The README says that backend "refuses vertical facets and
undercuts", and it does — `cmethod/_series.py` raises when a profile has two
heights at one position (`steps == 0.0`), and for `Lamellar` and
`Blazed(antiblaze_angle=90)` by type. A wall at 88.8° has a finite step, so it
passes the guard, and the solver then answers confidently and wrongly on a
boundary whose slope is 48. `check_energy_balance` catches it immediately, and
that is currently the only thing that does. Recorded in `docs/roadmap.md` as an
open question; the guard needs a slope ceiling rather than an equality test.

### One number explains the residual, and nothing can supply it

Every feature of the residual at the recorded geometry points at a single
scalar. The grating equation reads
`sin β_m = m λ / (period · sin γ) − sin α`
(`geometry.sin_beta`), so `period · sin γ` sets every diffracted angle and hence
every order's blaze wavelength. Table 1 and the lithography give
400 nm × sin 1.73° = 12.076 nm. The measurement wants 4.6% less.

The sharpest form of the evidence is not the RMS but the position of the blaze
order, which is a yes/no at each of the 37 energies:

| model | RMS | RMS × Névot–Croce | peak order matches |
|---|---:|---:|---:|
| measured groove, Table 1 geometry | 0.0383 | 0.0363 | **31 / 37** |
| the same, γ = 1.65° (period 400 nm) | 0.0157 | 0.0111 | **37 / 37** |
| the same, period 382 nm (γ = 1.73°) | 0.0142 | 0.0107 | **37 / 37** |

At the recorded geometry the model's peak order runs one 20 eV step ahead of the
measurement at 420, 520, 620, 640, 720 and 740 eV — always ahead, never behind.
Shrink `period · sin γ` by 4.6% and all 37 agree. The paper reports that "the
measured peak-order positions roughly match those predicted by the model", which
is what a one-step offset looks like on a plot of eight curves.

**No candidate for the missing 4.6% survives its own measurement.**

- *The period is not 382 nm.* The AFM reads 396 to 412 nm across three scans,
  with swath scatter of ±5 to ±8 nm, and round 3a's dispersion re-derives
  401 ± 8 nm. Round 3b's re-derivation gives 382 ± 14 nm, which is the value the
  comparison wants, but it is the noisier of the two and it disagrees with the
  AFM at nearly 2σ. A 4.5% error in an e-beam-written period is also not a
  plausible failure mode.
- *γ is not 1.65°.* It would need the fitted arc radius to be 6.71 mm against
  7.036 ± 0.124 (2.6σ low), or the throw to be 244 mm against the published
  232.0 ± 1.4 (8.8σ high) — and the throw is already the quantity this
  reduction had to move, to 233.0 mm. Nor is 1.65° some other recorded angle:
  the graze is 1.56° ± 0.04 and the small-angle combination of graze and yaw is
  1.706°.
- *It is not α.* Its published ±0.6° moves the model by 0.026 at most and
  changes the RMS by 0.0003.
- *It is partly, but not mostly, the facet angle.* The peak order depends on
  2δ − α as well as on `period · sin γ`. Holding the mount and the steep side at
  their measured values and scanning the facet gives an optimum at 26.5–27.0°
  with RMS 0.0283 and 35–36 of 37 peak orders — better than the AFM's own
  27.9° (0.0371, 33/37), and notably the paper's 27°. But all three AFM scans
  fit 27.9° or more, and the AFM's lateral scale reads 2% *high*, which biases a
  measured facet angle shallow, not steep. So the facet cannot be 1.4° below
  every scan, and in any case it gets 36/37 where the dispersion scale gets
  37/37.
- *It is not the beamline's energy scale.* λ enters only as
  `λ / (period · sin γ)`, so a 4.6% calibration error would be
  indistinguishable. The scan files carry a fourth column that
  `als_eff_meas.io.read_scan` documents as "the actual monochromator energy"
  and reads for nothing; it is a constant 500.3 in every file of both rounds,
  at every nominal energy from 80 to 800 eV, so it is not an energy readback
  and cannot test this. The nominal setpoints are all there is.

So the comparison has produced a sharp, falsifiable statement that the recorded
geometry cannot accommodate: the dispersion scale of this measurement is 4.6%
smaller than `400 nm × sin 1.73°`. Either one of the two published numbers is
wrong by more than its stated uncertainty, or something in the mount is not
captured by the pair `(period, γ)`. That is a measurement question, and it is
the one this campaign now hands back.

### What the solver's own checks say, and the one thing they caught

At the converged C-method setting, `check_energy_balance` returns
`R + A = 1.00000` at all 37 wavelengths, and RCWA at truncation 80 does the
same. Both compute the absorbed fraction from the solution rather than as
`1 − R`, so this is a real theorem and not an identity.

**The integral solver fails it, and the failure is a finding.** With
`conductivity="tabulated"` at 1024 nodes, `check_energy_balance` returns FAIL
with the total in [1.0001, 1.0827] at 37 of 37 points. The reflected orders are
not the problem — they match the C-method to 2 × 10⁻⁵ RMS over the whole
campaign — the independently computed absorption is:

| eV | \|n\| for Au | Σ R, integral − C-method | A, integral | A, C-method | R + A, integral |
|---:|---:|---:|---:|---:|---:|
| 800 | 0.9971 | −5.7e-6 | 0.34891 | 0.34883 | 1.00007 |
| 440 | 0.9946 | −3.3e-7 | 0.32667 | 0.32593 | 1.00073 |
| 200 | 0.9888 | −3.1e-7 | 0.22232 | 0.22012 | 1.00220 |
| 160 | 0.9762 | −1.4e-5 | 0.09029 | 0.08661 | 1.00367 |
| 120 | 0.9445 | −3.2e-6 | 0.06810 | 0.05627 | 1.01182 |
| 80 | 0.8880 | +3.4e-4 | 0.16107 | 0.07875 | 1.08266 |

The error grows as \|n\| leaves unity — not quite monotonically; there is
structure of a few parts in 10⁴ around the 340 eV region — and reaches 8% of the
incident power at 80 eV, where the integral solver's absorption is twice the
C-method's. **It also grows with refinement** — 3.9% at 512 nodes, 8.3% at 1024,
8.1% at 2048 — so it is not a discretisation error that more mesh would cure but
something in the absorption integral itself. It is also what spoils that
solver's convergence plateau: dropping the 80 eV point from the probe moves the
last two differences from 2.7e-4, 3.4e-4 to 1.3e-4, 1.9e-5.

Nothing in this entry's comparison depends on it, because the comparison uses
reflected orders and those are sound. But an absorption that is wrong by 8% is
wrong, it is invisible to a caller who does not run the energy check, and it is
recorded in `docs/roadmap.md` as an open question with the diagnosis above.

### The verdict on 0.3, from the TASTE campaign

**The comparison the criterion asks for is met. The solvers keep their
experimental status until the two defects this comparison found are fixed**:
the integral solver's absorption with tabulated metals
([#6](https://github.com/jamccoy/gratinglab/issues/6)) and the C-method on
near-vertical walls ([#7](https://github.com/jamccoy/gratinglab/issues/7)), both recorded below.

The gate asks for "use on a real problem, with the integral solver compared
against measured efficiencies ... and not only against published benchmarks and
its own physics checks", and the roadmap's gloss on it, written when the OGRE
comparison failed, says what the words were for: the gate "exists to establish
that the solver is right about a real grating", and "needs a measured groove for
the grating that was measured", because "a fit cannot separate a profile error
from a solver error".

That separation has now been made, and not by one solver's word:

- **The comparison is a prediction.** Period from the lithography and confirmed
  by dispersion, mount from the paper's Table 1, groove from the grating's own
  AFM scan through `metrology.boundary`, optical constants from a CXRO export,
  roughness from the same scan. No parameter was adjusted. 296 absolute
  per-order efficiencies over 80–800 eV and orders −4 to 8 come out at RMS
  0.0383, and 0.0363 with the roughness factor the paper applies.
- **The residual is demonstrably not the solver.** Three independent rigorous
  formulations agree with each other to 2 × 10⁻⁴ RMS and 7 × 10⁻⁴ at worst while
  all three sit 0.038 from the grating. The OGRE comparison could only assert
  that the solver was probably not the error, on the strength of evidence
  gathered elsewhere. This one measures it.
- **The absolute level is accounted for.** The sum over all measured orders
  divided by Fresnel gold tracks the Névot–Croce factor for the AFM's own
  1.5 nm rms to within 1.3% across the whole band, and applying that factor
  takes the mean residual from −0.0058 to −0.0007. What is left of the level is
  at the 1% mark and carries no structure.
- **What is left is one scalar, and it is a measurement question.** The whole
  residual is consistent with `period · sin γ` being 4.6% smaller than
  400 nm × sin 1.73°: supply that and the RMS falls to 0.014 and all 37
  blaze-order positions agree instead of 31. Neither published number can move
  that far within its own uncertainty, so the question goes back to the mount,
  not to the solver.

For comparison, the OGRE attempt reached RMS 0.0706 on a *fitted* two-parameter
groove, with order 0 holding the largest residual of every method, a duty-cycle
parity signature the fit could not satisfy, and no recorded uncertainty within
an eighth of the disagreement. This campaign reaches half that residual with
nothing fitted, models the zero order better than any other order, shows no
structural signature of a wrong profile family, and has its whole disagreement
reproduced by one input — the dispersion scale `period · sin γ` — moving 4.6%,
which is 2.7σ on the published γ or 1.3σ on the noisier of the two dispersion
re-derivations of the period.

**Three things this does not claim.**

It does not claim the solvers are validated to the measurement's precision. The
campaign's median σ is 3.0 × 10⁻⁵ and the residual is 1200 times that. What is
established is agreement at the level the grating's *geometry* is known, which
is as far as a comparison against a real grating can be taken without better
metrology on it, and the honest statement of the remaining uncertainty is the
0.028 to 0.063 bracket that the groove family spans.

It does not rest on the integral solver's absorption, which is wrong by up to 8%
at this geometry (above). The gate names the integral solver, and the integral
solver's *reflected orders* are what the comparison used and are sound to
2 × 10⁻⁵ against a converged independent method. Its absorption is a separate
defect and a separate open question.

It does not close the rounded-groove question for the 315 nm TASTE grating. That
entry — "The measured groove is rounded rather than faceted" — concerns a
different grating, and nothing here settles whether *its* rounding is the
grating or the probe. What this campaign adds is a second grating from the same
programme where the same question arises, and where the efficiency data answer
it: the apex is unpointed, and flat rather than rounded.

The remaining work for 0.3 is therefore not more solving. It is the release
mechanics, the two defects above, and the 4.6% in `period · sin γ`. The README's
"experimental" wording and the solver table that repeats it still say what they
said; they change at that release and not in this entry, which only records that
the criterion they point at has been satisfied.

**Released as 0.3 (2026-10-29).** Both defects are dealt with: #6 by the
opt-in `profile_smoothing="auto"`, which the energy-balance warning names when
it fires on a measured polygon, and #7 by a slope floor read from
the profile before solving, and the three rigorous solvers dropped
"experimental". The 4.6% in `period · sin γ` remains open: it is a question
about the measurement's geometry rather than a solver defect, and it is
recorded with the campaign.

Reproduced by `tests/test_measured.py::TestTasteCampaign`, which needs
`als_eff_meas` with both TASTE campaigns registered and its private beamline
tree at `$ALS_DATA_ROOT` and skips without either, and by
`tests/test_measured.py::TestTheMeasuredGroovesHandedness`, which needs neither
and always runs. The measurement is not in this repository and must not be
committed; nor are the AFM scans, which live in the group's own tree, so the
profile the committed tests use is the paper's published trapezoid rather than
the measured groove.

---

## Trapezoidal through the four solvers

`Trapezoidal` arrived with all four profile representations and no solver
having been run on it. Mapping it found one solver giving a wrong answer, one
profile defect, and one known fault in a new place. All three are recorded
here with the map itself.

The case is OGRE-like: period 160 nm, depth 68 nm (0.425 of the period), a
crest 0.30 of the period wide, CXRO gold (the vendored `Au.txt`), extreme
off-plane at graze 1.5° and azimuth 0, at 1.5, 3.0 and 4.5 nm. The sidewall
ran from 90° to 58°; below about 50.5° that groove has no floor left. A second, 24 nm deep groove (0.15 of the period)
reached 30°, which is 60° off vertical. Corner radius ran from 0 to the
largest that fits. Each solver climbed its own ladder: the integral solver at
400, 800 and 1600 boundary points (default `corner_grading=4`), the C-method at
truncation 16 to 128, and RCWA at truncation 10, 20 and 40, with its staircase
slices following at twice the truncation. 44 cases and 470 solves, measured
2026-10-29 on MC1 with `OPENBLAS_NUM_THREADS=1`. A perfect-conductor variant
ran wherever it applies, and RCWA refuses it, as it must.

### The C-method solved a sharp trapezoid wrongly

Before the fix, `cmethod/_series.py` sent only `Blazed` and `FromProfileData`
down its exact polygon path. A `Trapezoidal` took the FFT path meant for smooth
profiles. The FFT of a slope that jumps at four corners converges at first
order in the sample count: the slope coefficients were off by 8.5e-4 at 4096
samples and 1.9e-4 at 16384. So the answer sat on a floor that truncation could
not move. Against the same trapezoid given as `FromProfileData`:

| sidewall, depth | N = 16 | N = 32 | N = 64 |
|---|---:|---:|---:|
| 72°, 0.425 | 3.60e-4 | 3.81e-4 | 3.86e-4 |
| 45°, 0.20 | 1.99e-5 | 1.96e-5 | 1.94e-5 |

The integral solver sided with the polygon path. At 45° it was 3.5e-7 from it
at 1600 points against 1.9e-5 from the sampled answer, and at 72° 1.8e-5
against 3.4e-4 at 800 points. The corner warning never fired, because it read
polygons only, and the sampled path's energy balance stalled at 1.3e-3. A
reader had nothing to go on.

**Fixed.** A sharp trapezoid is now the polygon of its corners, with the
answer identical to `FromProfileData`'s and the warning given. A rounded one is
a new kind of series: the straight runs are integrated in closed form and each
fillet by Gauss-Legendre in the angle round its centre, where the integrands
are entire however steep the wall. A rule in position instead was still moving
by 7e-5 at four times its nodes on an 89° wall, because the slope has a
square-root singularity just past the arc's end there. In angle, quadrupling
the nodes changes the coefficients by under 3e-14 from 10° to 89.9°. A
2^22-sample FFT agrees to its own second-order error, 1e-10 at 72° with a
0.005-period fillet. `tests/test_cmethod.py::TestTrapezoidalSeries`.

### The map

Each entry is the largest |Δ efficiency| over orders and wavelengths, against
the integral solver at 800 points. That solver's own 400-to-800 change is the
first column, and on sharp corners it is 7e-7 to 9e-5. Its 1600-point rungs,
run on the gold sharp cases, move it a further 5e-6 to 1e-5.

| case (gold, TE) | int 400→800 | C 64 | C 128 | RCWA 10 | RCWA 20 | RCWA 40 | RCWA 40 − C 128 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 90°, sharp | 9.0e-5 | refused | refused | 6.4e-3 | 8.8e-4 | 1.2e-4 | |
| 89°, sharp | 8.8e-5 | 2.1e-1 | 6.8e-2 | 6.3e-3 | 8.8e-4 | 1.2e-4 | 6.8e-2 |
| 87°, sharp | 8.3e-5 | 1.9e-2 | 6.8e-4 | 5.9e-3 | 8.7e-4 | 1.3e-4 | 6.7e-4 |
| 85°, sharp | 8.0e-5 | 1.3e-3 | 1.7e-4 | 5.7e-3 | 9.0e-4 | 1.3e-4 | 1.7e-4 |
| 80°, sharp | 7.1e-5 | 1.5e-4 | 2.5e-5 | 6.2e-3 | 9.8e-4 | 1.5e-4 | 1.5e-4 |
| 76°, sharp | 6.3e-5 | 7.0e-5 | 1.4e-5 | 6.4e-3 | 1.0e-3 | 1.6e-4 | 1.5e-4 |
| 72°, sharp | 6.3e-5 | 5.0e-5 | 9.5e-6 | 7.2e-3 | 1.3e-3 | 2.4e-4 | 2.4e-4 |
| 65°, sharp | 4.6e-5 | 1.8e-5 | 5.6e-6 | 7.4e-3 | 1.3e-3 | 2.2e-4 | 2.1e-4 |
| 58°, sharp | 4.2e-5 | 5.4e-6 | 5.4e-6 | 7.2e-3 | 1.4e-3 | 2.7e-4 | 2.6e-4 |
| 24 nm deep, 45°, sharp | 8.2e-6 | 3.5e-6 | 2.9e-6 | 1.4e-3 | 2.3e-4 | 3.5e-5 | 3.3e-5 |
| 24 nm deep, 30°, sharp | 7.8e-6 | 5.0e-7 | 7.9e-7 | 1.4e-3 | 2.5e-4 | 4.9e-5 | 5.0e-5 |

TM, run at 90°, 85° and 72°, is within 10% of TE in every column. Each solver costs, for
three wavelengths: integral 800 points about 16 s, RCWA truncation 40 about
20 s, C-method truncation 128 about 3 s. The integral solver at 1600 points
took 30 to 55 minutes under load, which is why the ladder stops at 800.

**Where all three converge, they agree.** On every sharp case from 58° to 85°
the integral solver is converged to about 1e-5 at 800 points, the C-method at
128 is within 2e-4 of it, and RCWA closes on it at about 5.5 times per
doubling of truncation. RCWA 40's gap to the C-method equals RCWA's own
remaining error, so nothing is left over that would be a disagreement about
the physics.

**RCWA does not care about the sidewall.** Its error at a given truncation is
the same within a factor of two from 58° to 90°: the staircase steps are
vertical whatever the true wall does. What it pays for is the staircase
itself, which converges roughly at first order. A vertical, sharp
`Trapezoidal` is also solved as a staircase, at about twelve times the cost of
`Lamellar`'s single exact slab. Every slab is identical, so the answer still
agrees to rounding.

**The C-method's practical limit on this groove is about 85°.** There it needs
truncation 128 for 1.7e-4. At 87° it reaches only 6.8e-4 at 128, and at 89°
it is wrong by 6.8e-2 at 128 and still moving by 0.15 between 64 and 128.
Its energy balance shows each of those failures: 4e-4 at 85°, 2e-3 at 87° and
4e-2 at 89°, all at truncation 128. Below 80° it is the fastest of the three
by an order of magnitude and the best converged. Branch `cmethod-steep-walls`
(issue #7) derives a truncation floor that refuses the hopeless cases; it
measures the same knee independently on the TASTE groove.

**Perfect conductor:** the integral solver at 800 points against the C-method
at 128 gives 5.5e-5 at 58° and 1.3e-4 at 72°. The C-method converges more
slowly without absorption, at 1.4e-4 for 72° in TM and 4.8e-4 at 85°. The
integral solver balances to 1e-11 in TE.

### Rounded corners, where the integral solver's absorption fails

| case (gold, TE, 72°) | radius | int bal 800 | C 64→128 | C bal 128 | RCWA 40 − C 128 | int 800 − C 128 |
|---|---:|---:|---:|---:|---:|---:|
| sharp | 0 | 9e-7 | 5.0e-5 | 1e-5 | 2.4e-4 | 9.5e-6 |
| 0.3 nm | 0.002 | 2 | | 9e-6 | 2.7e-4 | 1.0e-2 |
| 0.8 nm | 0.005 | 0.2 | | 5e-6 | 4.0e-4 | 3.5e-3 |
| 1.6 nm | 0.01 | 0.04 | 1.6e-5 | 4e-6 | 4.7e-4 | 1.4e-3 |
| 3.2 nm | 0.02 | 0.01 | | 2e-6 | 6.0e-4 | 6.0e-4 |
| 8 nm | 0.05 | 2e-3 | 3.4e-6 | 6e-7 | 9.3e-4 | 1.2e-4 |

With the corners rounded the C-method becomes the reference. The slope is
continuous, so it converges quickly and balances to a few 1e-6. RCWA agrees
with it to within RCWA's own staircase error, which is 1.5e-3 between
truncation 20 and 40 at 1.6 nm.

The integral solver's R + A − 1 reads 2 at 0.3 nm and 0.2 at 0.8 nm. It falls
with the mesh far too slowly to be convergence, and its reflected orders are
off by up to 1e-2. A fillet reaches it as a polyline of one vertex per degree,
every vertex turning by less than `CORNER_TURN_DEGREES`. No corner is detected
and the mesh is equal-arc and ungraded, while each vertex is still a real,
weak corner of the metal interface. That is exactly the fault diagnosed on
measured AFM grooves (issue #6, branch `integral-absorption`), met here on an
analytic profile. A perfect conductor has no absorption to get wrong, and on
the same 1.6 nm fillet the integral solver and the C-method agree to 2.8e-5.
The solver does warn that the summed efficiency strays from unity. Issue #6
needs a design decision, so it was not fixed here.
`tests/test_trapezoidal_solvers.py::TestRoundedCorners` pins that the failure
is reported; it is written to fail once the fault is fixed.

### The largest fillet had a phantom corner

At `max_corner_radius` on the 72° groove, the crest's straight run has zero
width. `Trapezoidal.boundary` appended its end point anyway, a second copy of
the previous arc's end to within rounding. The segment between them, about
1e-17 long and pointing anywhere, read as a corner of more than 10°. The graded
mesh then clustered two nodes onto it, and the integral solver refused every
`boundary_points` with "puts two nodes 0 nm apart". Fixed by skipping the
zero-width run. The integral solver now solves that groove, though its
absorption still fails as in the section above, and its reflected orders are
within 6e-4 of the C-method. The C-method itself loses digits on that profile
from truncation 32 (condition number 3e12, 1.5e20 at 128) and says so.
`tests/test_trapezoidal_solvers.py::TestTheLargestFillet`.

### The reduction to Lamellar

With `sidewall_angle=90` and `corner_radius=0`, against `Lamellar` of the same
depth and duty:

| solver | gold | perfect conductor |
|---|---|---|
| RCWA, truncation 10 / 20 | 2.0e-13 / 1.2e-12 (TE), 6.1e-13 / 9.5e-13 (TM) | refused |
| integral, 400 / 800 points | 6.0e-6 / 1.1e-7, against Lamellar's own 400-to-800 change of 9.0e-5 | 4.4e-16 (TE), 2.5e-13 (TM) |
| C-method | refused: a vertical wall has no finite slope | refused |

The integral solver's gold residual is mesh error, not a disagreement. The two
polygons are the same shape translated by half a period, so their graded meshes
meet the period boundary in different places. The perfect-conductor case shows
the same geometry to rounding.
`tests/test_trapezoidal_solvers.py::TestVerticalAndSharpIsLamellar`.
