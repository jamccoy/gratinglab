# Findings

Empirical results established while building this, each with the evidence that
established it. Recorded so they are not rediscovered, and so a claim can be
checked rather than taken on trust.

Findings that constrain how code must behave are also encoded as tests; those
are cross-referenced.

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

## Scalar theory does not conserve energy

With the order-dependent phase that both primary references specify, summed
efficiency **exceeds unity** by up to ~12 % across mounts (~7 % off-plane).
A passive grating cannot return more power than it receives, so these values are
unphysical.

This is a property of the formulation, not of the implementation. For a fixed
phase, $\sum_m \operatorname{sinc}^2(x - m) = 1$ is an exact identity; making
$x$ depend on the order breaks it. Running the identical machinery with
`phase_reference="specular"` restores a genuine Parseval pair and satisfies
$\sum_m|G_m|^2 = 1$ to $10^{-10}$ — which is what proves the quadrature and
normalisation are correct.

Full derivation: [theory/scalar.md §5](theory/scalar.md#5-energy-is-not-conserved-and-that-is-the-formulations-fault).
The excess is reported as a provenance warning and never rescaled.

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
hypothetical. Now pinned to the exact ratio, with a second test proving the two
ratios are actually distinguishable for the geometry used, so the assertion
cannot go vacuous.

### In-plane tests are blind to γ errors

Breaking $\sin\gamma \to \sin^2\gamma$ left **34 of 37** scalar tests passing.
In-plane cases have $\sin\gamma = 1$, which makes any power of it identical, so
they cannot detect an error in the γ handling at all. Only the off-plane cases
failed — and off-plane is the primary application.

Coverage at the time was 9 off-plane tests against 13 in-plane. Every
closed-form check now runs in **both** mounts.

After both fixes, all 12 mutations are caught.

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
- The obliquity factor differs between them and is not a matter of notation —
  the two redistribute efficiency between orders differently.
