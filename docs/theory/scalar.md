# Scalar (Kirchhoff) diffraction

Implemented in [`solvers/scalar.py`](../../src/gratinglab/solvers/scalar.py).
Symbols and sign conventions follow [`conventions.md`](../conventions.md); every
equation below is written in those.

---

## 1. What the model assumes

Scalar theory replaces the six coupled vector wave equations with one scalar
equation. That is only defensible under a specific set of assumptions, and the
whole point of this page is that they are stated rather than implied.

| Assumption | Consequence when it fails |
|---|---|
| The field is scalar — polarization is neglected | TE and TM come out identical. Real gratings separate them, sometimes by a lot. |
| Structure is large compared to the **reduced** wavelength, λ/sin γ ≪ p | Sub-wavelength features are misrepresented; the model degrades smoothly rather than failing loudly. In-plane (sin γ = 1) this is the familiar λ ≪ p; in a conical mount the grooves pose a transverse problem at λ/sin γ (the same decoupling the integral solver runs on), and judging validity on λ/p alone passes exactly the extreme off-plane mounts it should flag — see §7. |
| Kirchhoff (thin-element) boundary condition: the field just above the surface is the incident field times a local phase | No multiple scattering between groove facets, no field penetration into the material — that is what the [integral method](integral.md) exists for. Geometric shadowing *is* modelled (§9): self-shadowing by the local facet normal always, and cast shadows — one part of the groove occluding another along the ray — under `visibility="horizon"`. |
| Fraunhofer (far-field) observation | Valid for a telescope focal plane; not for near-field work. |
| Reflectivity is resolved **across the groove cycle** and carried inside the diffraction integral | Name a `coating` and efficiencies are **absolute**; leave it unset and they are **relative**, which is a correct answer to a different question rather than a deficiency. §9 has the three models and what each costs. |

> **This table used to say something else.** Until M16 the fifth row read
> "reflectivity applied separately, as one scale factor per wavelength", on the
> grounds that "scalar theory has no mechanism by which one order could reflect
> differently from its neighbour". There is such a mechanism, and the row was
> wrong. A groove whose reflectivity varies across the cycle is an **amplitude**
> grating as well as a phase grating; the two Fourier-transform together, and
> the result is order-dependent. On the reference geometry that is worth −51% to
> +12% order by order. See §9.

Because polarization is neglected, the solver advertises `rigorous=False` in
its [`Capabilities`](../../src/gratinglab/solvers/base.py) and records a
provenance warning whenever a result is labelled TE or TM.

---

## 2. The phase a groove imparts

This is the step where a factor is easiest to lose, so it is derived rather than
quoted.

Work in the frame of `conventions.md` §3: `d̂ = x̂` along the periodicity, `n̂ = ŷ`
the grating normal, grooves along `ĝ = −ẑ`. The wave vectors are

$$\mathbf{k}_i = k\left[-\sin\alpha\sin\gamma,\; -\cos\alpha\sin\gamma,\; \cos\gamma\right]$$

$$\mathbf{k}_m = k\left[\;\;\sin\beta_m\sin\gamma,\; \;\;\cos\beta_m\sin\gamma,\; \cos\gamma\right]$$

with $k = 2\pi/\lambda$. Note the normal components:

$$k_{y,i} = -k\cos\alpha\sin\gamma \qquad k_{y,m} = +k\cos\beta_m\sin\gamma$$

The incident one is negative — the wave travels *toward* the grating
(`conventions.md` §6).

α, β_m and γ are not independent — periodicity along $\hat{x}$ ties them together as
the **generalized (conical) grating equation**:

$$\boxed{\;\sin\alpha + \sin\beta_m = \frac{m\lambda}{p\,\sin\gamma}\;}
\qquad m = 0, \pm1, \pm2, \dots$$

with **α ∈ (−90°, 90°)** and **γ ∈ (0°, 90°]** — γ the polar angle from the groove axis
(it sweeps out the cone), α the azimuthal position around it. `sin γ = 1` (γ = 90°) is the
classical in-plane case; this is where the derivation below reduces to the familiar
in-plane phase. See `conventions.md` §3–§4 for the full derivation, the range
justification, and how these names translate to other papers' θ/φ.

Now consider a ray striking the surface at height $g(x)$ rather than at a flat
reference $y = 0$. Relative to that reference it accumulates $k_{y,i}\,g$ on the
way in and $k_{y,m}\,g$ on the way out, so the round trip contributes

$$\Phi_m(x) = \left(k_{y,m} - k_{y,i}\right) g(x)
            = k\,g(x)\,\sin\gamma\,\left[\cos\alpha + \cos\beta_m\right]$$

**This is where `sin γ` enters**, and it enters through the wave-vector
components rather than by assertion. The thesis (`Appendix-D.tex:651`) writes
this expression without `sin γ`, which is the in-plane special case; the ISSI
chapter's eq. (15) carries it. The general form above is what the code
implements — see errata item 1 in `conventions.md` §10.

**And this is where the order dependence enters.** $\Phi$ depends on $m$ through
$\cos\beta_m$, because different orders leave along different directions and
therefore accumulate different path lengths on the way out. Section 5 shows what
that costs.

---

## 3. From phase to efficiency

Write the grating's action on the field as a **transmittance function**

$$G(x) = \mathcal{A}(x)\, e^{i\Phi(x)}, \qquad G(x + p) = G(x)$$

where $\mathcal{A}$ is any amplitude modulation (unity for a pure phase grating
with perfect reflectivity). Under Fraunhofer observation the far field is the
Fourier transform of $G$, and periodicity turns that transform into a Fourier
*series*: power appears only at the discrete orders satisfying the grating
equation, $k_{x,m} - k_{x,i} = mK$ with $K = 2\pi/p$.

The coefficient of order $m$ is therefore

$$G_m = \frac{1}{p}\int_0^p G(x)\, e^{-imKx}\,dx$$

Substituting $t = x/p$ (the normalised coordinate the
[`Profile`](../../src/gratinglab/profiles.py) classes use) gives the form the
code evaluates:

$$\boxed{\;G_m = \int_0^1 e^{i\Phi_m(t)}\, e^{-2\pi i m t}\,dt\;}
\qquad
\boxed{\;\mathscr{E}_m = O_m\left|G_m\right|^2\;}
\qquad
O_m = \frac{4\cos\alpha\cos\beta_m}{(\cos\alpha + \cos\beta_m)^2}$$

Two implementation notes follow directly:

- Because $\Phi_m$ carries an $m$, this is **one integral per order**, not a
  single FFT. The cost is negligible and the alternative would be wrong.
- The integrand is periodic, so the rectangle rule converges *spectrally*. This
  is why `quadrature_points` reaches machine precision at modest values, and why
  it must still exceed $2|m|_{\max}$ to satisfy Nyquist.

### The obliquity factor: which one, and how we know

Thesis Appendix D carries an obliquity factor $\cos\beta_m/\cos\alpha$, from
projecting the Poynting flux onto a plane parallel to the surface, and then
renormalises so that $\sum_m \mathscr{E}_m = 1$.

**The renormalisation is simply an error** and is not offered as an option:
forcing the sum to unity erases the model's own signal about how far it has
strayed from its regime (§5).

**The obliquity factor was rejected too, and that went one step too far.** The
Appendix D form is *asymmetric* under $\alpha \leftrightarrow \beta_m$, so it
breaks Lorentz reciprocity — which is why this page rejected it, correctly. The
conclusion drawn from that, that no flux factor belongs at all, does not follow.

The factor that does belong is the **symmetric** one, $O_m$ above, and it is
settled by a theory derived from a different starting point. Expand the field in
Rayleigh plane waves, impose the boundary condition on the corrugated surface,
keep first order in groove height; for a Dirichlet (perfectly conducting, TE)
surface that gives

$$\eta_m = 4\,k_{z,0}\,k_{z,m}\left|\hat{g}_m\right|^2, \qquad
k_z = k\cos\theta\sin\gamma$$

Kirchhoff theory needs gentle slopes, perturbation theory needs shallow grooves,
and a grating that is both must be described correctly by both. In that overlap
the bare $|G_m|^2$ is wrong by exactly $1/O_m$ — verified to five significant
figures across in-plane and conical mounts in
[`test_perturbation.py`](../../tests/test_perturbation.py):

| Mount | $\lambda/p$ | $m$ | $\beta_m$ | $\vert G_m\vert^2$ / exact | $O_m$ |
|---|---|---|---|---|---|
| in-plane, $\alpha=10°$ | 0.80 | $+1$ | $38.8°$ | 1.0137 | 0.9865 |
| in-plane, $\alpha=10°$ | 0.80 | $-1$ | $-76.8°$ | **1.6375** | 0.6107 |
| off-plane, $\alpha=20°$, $\gamma=1.25°$ | 0.0127 | $+1$ | $13.9°$ | 1.0003 | 0.9997 |
| off-plane, $\alpha=20°$, $\gamma=1.25°$ | 0.0127 | $-1$ | $-67.5°$ | **1.2156** | 0.8227 |

Negligible near Littrow, 22% in this project's own grazing-incidence regime, 64%
at high dispersion. Three properties make $O_m$ admissible where the Appendix D
form was not:

- **Symmetric** in $\alpha \leftrightarrow \beta_m$, so reciprocity (§6) is
  untouched — still $10^{-16}$.
- **$\leq 1$**, by AM–GM, so it can only reduce an efficiency. It helps the
  energy bound rather than threatening it.
- **Exactly 1** when $\cos\alpha = \cos\beta_m$: Littrow, specular, and every
  $m=0$. The shallow-groove energy identity is therefore unchanged, and a
  perfect blaze *in Littrow* still reaches exactly unity.

Away from Littrow a perfect blaze now reaches $O_m$ rather than 1. That is a
correction, not a loss: unity-at-blaze was a property of the unfactored
$|G_m|^2$, and a rigorous calculation does not give exactly 1 there either.

> **Which polarization this matches.** The Neumann (TM) first-order result is
> $4k^2|\hat{g}_m|^2(1+\sin\alpha\sin\beta_m)^2/(\cos\alpha\cos\beta_m)$, which
> departs from the Dirichlet one away from Littrow — by a factor of 14 at the
> in-plane geometry above. A polarization-blind model cannot match both. TE is
> the conventional scalar correspondence, and the size of the TE/TM spread is a
> usable measure of where scalar theory stops meaning anything at all.

---

## 4. The closed forms are consequences, not code paths

Three profiles admit exact answers. **They are implemented nowhere** — the
solver evaluates the general integral above and the closed forms serve as tests
([`test_scalar.py`](../../tests/test_scalar.py)). Getting all three right from
one code path is the evidence the quadrature is sound.

Throughout, $\varphi_m \equiv k\,h\,\sin\gamma\left[\cos\alpha + \cos\beta_m\right]$
is the phase across the full groove depth $h$.

### Sawtooth (blazed)

$g(t) = h\,t$, so $\Phi_m(t) = \varphi_m t$ and

$$G_m = \int_0^1 e^{i(\varphi_m - 2\pi m)t}\,dt
      = \frac{e^{i(\varphi_m - 2\pi m)} - 1}{i(\varphi_m - 2\pi m)}$$

$$\mathscr{E}_m = \frac{2\left[1 - \cos(\varphi_m - 2\pi m)\right]}{(\varphi_m - 2\pi m)^2}
                = \operatorname{sinc}^2\!\left(\frac{\varphi_m}{2} - m\pi\right)$$

Peak efficiency in order $m$ occurs at $\varphi_m = 2\pi m$, which is the blaze
condition and reduces to $\beta_m = 2\delta - \alpha$ — see
`conventions.md` §9 and `geometry.blaze_wavelength`.

### Lamellar (binary)

$g(t) = h$ for $t < w$ and $0$ otherwise, with duty cycle $w$. For $m \neq 0$

$$G_m = \frac{\left(e^{i\varphi_m} - 1\right)\left(e^{-2\pi i m w} - 1\right)}{-2\pi i m}
\qquad
\mathscr{E}_m = 4\sin^2\!\left(\frac{\varphi_m}{2}\right) w^2 \operatorname{sinc}^2(mw)$$

The $\operatorname{sinc}^2(mw)$ factor is why a 50 % duty cycle extinguishes every
even order. The $m = 0$ case must be stated separately —
$G_0 = w e^{i\varphi_0} + (1 - w)$ — since the expression above is singular
there (errata item 4, `conventions.md` §10).

### Sinusoid (holographic)

$g(t) = \frac{h}{2}\left[1 - \cos 2\pi t\right]$, so
$\Phi_m(t) = \frac{\varphi_m}{2}\left[1 - \cos 2\pi t\right]$. Applying the
Jacobi–Anger expansion to $\exp\!\left[-i\frac{\varphi_m}{2}\cos 2\pi t\right]$,

$$G_m = e^{i\varphi_m/2}\,(-i)^m\,J_m\!\left(\frac{\varphi_m}{2}\right)
\qquad
\mathscr{E}_m = J_m^2\!\left(\frac{\varphi_m}{2}\right)$$

This one is worth having because it is **independent of Appendix D** — it comes
from a Bessel identity rather than from the same source as the other two, so it
checks the machinery from a different direction.

---

## 5. Energy is not conserved, and that is a deliberate trade

This is the most important caveat on the page, and the reasoning behind it is
not obvious.

For a **fixed** phase there is an exact identity

$$\sum_{m=-\infty}^{\infty} \operatorname{sinc}^2(x - m) = 1 \quad \text{for any } x$$

so the sawtooth result would conserve energy exactly. But $\Phi_m$ carries an
$m$ (§2), which makes $x = \varphi_m/2\pi$ depend on the summation index and
**breaks the identity**. Equivalently: the $G_m$ are not the Fourier
coefficients of any single function, so Parseval does not apply.

The obvious repair is to drop $\cos\beta_m$ from the phase, leaving
$\Phi = 2kg\sin\gamma\cos\alpha$ — the round-trip path of the incident wave,
with no diffracted direction in it. That does restore a genuine Fourier pair and
conserves energy to $10^{-10}$. **It is not what this solver does**, and the
reason is a three-way tension no scalar formulation escapes.

### The tension, measured

| Property | order-dependent $\cos\alpha + \cos\beta_m$ (this solver) | β-free $2\cos\alpha$ |
|---|---|---|
| Energy, $\sum_m \mathscr{E}_m \leq 1$ | ✗ violated, up to 1.61 | ✓ exact (Parseval) |
| Reciprocity, $\mathscr{E}_m(\alpha) = \mathscr{E}_m(\beta_m)$ | ✓ exact, $10^{-17}$ | ✗ violated, 0.44 |
| Blaze direction $\beta_b = 2\delta - \alpha$ | ✓ exact at every $\alpha$ | ✓ at Littrow only; 7.6° off at 10° from Littrow, and past ~15° the envelope peak goes evanescent |

The conflict is structural. Reciprocity **requires** the phase to be symmetric
under $\alpha \leftrightarrow \beta_m$, which forces both directions into it.
Parseval **requires** a single fixed function, which forbids $\beta_m$. They
cannot both hold. Physical optics is reciprocal but not energy-conserving; the
pure transmittance-function picture is the reverse.

**Reciprocity is the invariant this solver keeps**, so the energy deviation is
expected rather than repairable.

> A tempting third option is to reference the phase to the blaze direction
> $\beta_b = 2\delta - \alpha$, which is β-free *and* reproduces the blaze
> wavelength exactly. It is circular: $\beta_b$ is an **output** of the model —
> it is where the sinc² envelope peaks — so feeding it back in as an input
> assumes the answer. The apparent perfect agreement is a tautology, not a
> prediction.

### The implementation is not at fault

In the shallow-groove limit the transmittance tends to 1, all power lands in
order 0, and the sum must tend to 1 for any formulation. It does, exactly:

| depth / period | $\sum_m \mathscr{E}_m$ |
|---|---|
| 0.0001 | 1.00000 |
| 0.0010 | 0.99997 |
| 0.0200 | 0.98836 |
| 0.0500 | 0.94603 |
| 0.1000 | 0.92902 |

(Remeasured in M16-C. The obliquity factor of §3 is exactly 1 at $m=0$, where
all the power sits in the shallow limit, so the top row is still exactly 1 —
which is the point of the check. The worst-case excursion over a wide mount
sweep fell from 1.71 to 1.61, since $O_m \leq 1$ can only ever pull the sum
down.)

That is the check which would catch a genuine coding error, as distinct from the
formulation's known defect. The deviation scales with **phase excursion across
the groove** — depth relative to wavelength, hence working diffraction order —
*not* with λ/p and *not* with the number of propagating orders. This is why the
sum looks stubbornly constant across wildly different mounts at fixed profile,
and it is what the provenance warning reports.

The deviation is **reported and never rescaled**. Renormalising by the sum
(as thesis Appendix D does) would erase precisely the signal that says how far
the model has strayed from its regime.

---

## 6. Reciprocity — the check that constrains the model

Lorentz reciprocity requires that reversing the diffracted ray of order $m$
gives an incidence azimuth $\alpha' = \beta_m$ from which order $m$ exits along
$\alpha$:

$$\mathscr{E}_m(\alpha) = \mathscr{E}_m(\beta_m)$$

Scalar theory satisfies this **exactly**, and the reason is visible in §2: under
$\alpha \leftrightarrow \beta_m$ the phase $\cos\alpha + \cos\beta_m$ is
symmetric, and the Fourier kernel is unchanged, so $G_m$ is identical.

This makes it a far stronger test than the closed forms. Those compare the
solver against a formula derived the same way the solver computes it, so they
validate the *quadrature*. Reciprocity constrains the *structure* of $\Phi$ and
would fail immediately for $2\cos\alpha$, or for a subtraction. Measured with
[`checks.check_reciprocity`](../../src/gratinglab/checks.py):

| | max violation |
|---|---|
| as implemented | 5 × 10⁻¹⁸ |
| with the α↔β symmetry of Φ broken | 4 × 10⁻¹ |

Sixteen orders of magnitude of discrimination, needing no reference data.

---

## 7. Validity conditions

Reported on `Provenance.warnings` rather than enforced, because mapping where
the model breaks down is a deliverable in its own right.

| Condition | Meaning |
|---|---|
| $\lambda / (p \sin\gamma) \lesssim 0.1$ | Kirchhoff theory assumes structure ≫ wavelength — measured at the **reduced** wavelength $\lambda/\sin\gamma$, because that is the wavelength of the transverse problem the grooves actually pose (M&P eq. 4.65). In-plane this is the familiar $\lambda/p$. The distinction is not academic: the flagship off-plane geometry reads $\lambda/p \approx 0.007$ but a reduced ratio of $0.32$, and there the scalar-vs-integral disagreement is orderwise catastrophic — `tests/test_cross_method.py` holds $\lambda/p$ fixed and closes the cone to show the discrepancy growing along exactly this axis. Order positions and the blaze envelope stay useful past the threshold; per-order numbers need a rigorous cross-check. |
| $\zeta < \theta_c \approx \sqrt{2\,\text{decrement}}$ | Facet graze must stay below the critical angle for total external reflection, or reflectivity collapses. Evaluated when a `coating` is named; without optical constants there is nothing to compare against, so it is a check that does not apply rather than a warning. |
| $\lambda > 32 \sin(\zeta)\,\sigma$ | Fraunhofer smoothness: below this the facet is not optically smooth and scatter dominates. Checked for every profile — it was once gated on having a blaze angle, so a rough sinusoid was never checked at all. |
| No order fully shadowed | With reflectivity resolved across the groove (§9), an order into which no lit facet can radiate comes out at exactly zero. That is the model's answer, not the physical one, and it is reported by name so it is never mistaken for passing-off. |
| $r_p$ stays away from zero on the lit groove | At Brewster the p amplitude changes sign and its phase jumps by $\pi$; the geometric-mean symmetrisation of §9 does not carry that cleanly. Only reachable with steep grooves near normal incidence. |
| $R_s \approx R_p$ at the guarded graze (split $\lesssim 2\%$) | The 50/50 s/p average of §9 stands in for a groove-TE/TM → facet-s/p mapping nobody carries, and any basis rotation can move the result by at most half the split. At the grazes this project runs the split is under 1% (see §9); a steep groove near normal incidence pushes it to tens of percent, and the waiver is reported there instead of staying silently priced at zero. |
| $\sum_m \mathscr{E}_m \leq 1$ | Not a validity condition so much as a sanity bound; see §5. |

Both $\zeta$-based conditions are evaluated at the **worst local graze** across
the groove when the reflectivity model resolves it, and at the single facet
angle when it does not. A guard should describe the model that actually ran.

---

## 8. How many quadrature points

$G_m$ is evaluated as a rectangle rule on a uniform grid over one period, which
for a *periodic* integrand is the trapezoid rule and is spectrally accurate —
provided the integrand is smooth. Whether it is smooth depends entirely on the
profile:

`converged_at` below is what `check_convergence` reports at a $10^{-6}$
tolerance — the *coarsest* setting it could demonstrate adequate, which is what
a production run should use:

| Profile | Convergence in $n$ | `converged_at` |
|---|---|---:|
| Sinusoidal | spectral — machine precision from $n = 64$ | 256 |
| Lamellar | $O(n^{-2})$, clean 4× per doubling | 1 024 |
| Blazed | $O(n^{-2})$ asymptotically, **non-monotone below ~8 000** | 4 096 |

### Resolving reflectivity across the groove costs quadrature

The visibility mask of §9 puts a **jump** in the integrand at the shadow
boundary, where a pure phase grating had at worst a slope kink. A jump converges
at $O(n^{-1})$, not $O(n^{-2})$, and on a blazed profile the shadow boundary
sits at the apex — $t = 0.8331$, not a dyadic rational — so the error is
non-monotone on top of being slower. Measured on the reference geometry
(Blazed 29.5°/70.5°, Au, $\gamma = 1.25°$):

| `reflectivity_model` | `converged_at`, $10^{-4}$ | $10^{-5}$ | $10^{-6}$ |
|---|---:|---:|---:|
| `facet` | 256 | 512 | 4 096 |
| `average` | 2 048 | 8 192 | not reached by 131 072 |
| `local` | 2 048 | 32 768 | not reached by 131 072 |

**The default $n = 2048$ buys about $10^{-4}$ with a resolved model under
`visibility="facet-normal"`, not $10^{-6}$.** A smooth profile is unaffected —
a coated sinusoid still converges at 256, because it has no shadow boundary to
resolve. This is a real cost of the change and the honest way to spend it is on
the convergence harness rather than on a larger default: `check_convergence`
reports what a given case actually needs, and a case that never plateaus says
so.

**`visibility="horizon"` removes the penalty rather than paying it.** Its
shadow boundaries are geometric crossings, and
[`geometry.horizon_weights`](../../src/gratinglab/geometry.py) places each one
*inside* its cell — exactly, on a polygonal profile — with the lit-side sample
absorbing the sub-cell lit length and the incident and exit masks composed by
minimum (they share the apex boundary; a product counts it twice). Measured on
the coated blazed reference against a 16× finer run: $2\times10^{-7}$ at the
default $n = 2048$, worst rung of the whole ladder $1.3\times10^{-5}$ — against
$3.4\times10^{-4}$ for the same case with binary masks, and $5\times10^{-4}$
for `facet-normal`. `facet-normal` keeps its binary masks deliberately: its
boundary function jumps at profile corners, where a crossing estimate is biased
rather than sharp, and bit-for-bit reproducibility is that mode's purpose. So
the horizon mode is both the more physical treatment (§9) and the numerically
sharper one.

(256 for the sinusoid is the ladder's first rung, not a requirement — it was
already converged before the sweep began.)

The blazed case is the awkward one, and it is the primary application. The
groove has a slope discontinuity at the apex, and the apex sits at

$$t_{\text{apex}} = \frac{1}{1 + \tan\delta / \tan\delta'} = 0.8331$$

for the reference geometry — not a dyadic rational, so doubling $n$ changes how
near a node lands to the kink and the error does not fall monotonically. It can
*grow* by a factor of 2.5 on a refinement step. See
[`findings.md`](../findings.md#quadrature-error-is-not-monotone-in-the-number-of-points)
for the measured ladder.

Two consequences:

- **The default $n = 2048$ is not a converged setting for a blazed profile.**
  It is a reasonable starting point. `gratinglab.convergence.check_convergence`
  finds 4 096 adequate at a $10^{-6}$ tolerance, and proves it by sweeping to
  16 384.
- **A single agreement between two refinements is not evidence.** The harness
  requires three consecutive values to agree, which is what rejects the
  spurious plateau at $n = 512$.

There is also a hard floor, enforced rather than reported: $n$ must exceed
$2\max|m|$, or the highest order is aliased. That is Nyquist on the kernel
$e^{-2\pi i m t}$, and no accuracy argument can rescue a grid below it, so
`solve` raises instead of returning a plausible number.

---

## 9. Absolute efficiency: reflectivity across the groove cycle

A groove does not present one angle to the beam. Its local facet tilt
$\delta(t)$ varies across the period, so the local graze does too:

$$\sin\zeta(t) = \sin\gamma\,\cos(\delta(t) - \alpha), \qquad
\tan\delta(t) = +\frac{dy}{dt}$$

Three models are offered, selected with `reflectivity_model`.

### `"local"` — the default

$$\mathscr{E}_m = \frac{O_m}{2}\left(
\left|\int_0^1 r_s^{\text{gm}}(t)\,e^{i\Phi_m(t)}e^{-2\pi imt}\,dt\right|^2 +
\left|\int_0^1 r_p^{\text{gm}}(t)\,e^{i\Phi_m(t)}e^{-2\pi imt}\,dt\right|^2
\right)$$

with the complex Fresnel amplitude carried **inside** the integral. Because
$r(t)$ varies, the groove is an amplitude grating as well as a phase grating,
and the reflectivity becomes **order-dependent** — the thing §1 used to say was
impossible.

**Why a geometric mean.** The obvious weight $r(\zeta_i(t))$ depends on $\alpha$
alone, which destroys the $\alpha \leftrightarrow \beta_m$ symmetry and breaks
reciprocity — by 82%, measured. Weighting each facet by both the angle it
receives at and the angle it emits at,

$$r^{\text{gm}}(t) = \sqrt{r(\zeta_i(t))}\;\sqrt{r(\zeta_{d,m}(t))},
\qquad \sin\zeta_{d,m}(t) = \sin\gamma\,\cos(\delta(t) - \beta_m)$$

is symmetric by construction, holds reciprocity to $6\times10^{-15}$, and
collapses to plain $r(\zeta)$ wherever the facet is at specular — so a perfectly
blazed groove is left alone by the change.

**The two roots are taken separately**, not as $\sqrt{r_i r_d}$, so the only
discontinuities in the weight are the ones $r$ itself has. Measured, this is a
cheap defensive choice rather than a load-bearing one: where the product's
argument wraps *uniformly* across the groove — the grazing case — the two forms
differ by a global sign that cancels out of $|\int\cdots|^2$, and on Au at
1–6 nm they agree to $10^{-15}$. They diverge only where the wrap is
non-uniform, near normal incidence on a deep groove, where $\arg r_p$ sweeps
through Brewster; there the gap reaches 12%, on orders carrying $10^{-9}$, in a
regime already flagged below.

**Visibility.** Where the groove is shadowed in either direction it contributes
nothing; which shadows the masks can see is its own decision, `visibility`, with
its own subsection below.

**Per-order gain is expected.** Since $|r^{\text{gm}}| \leq 1$ pointwise but the
bare integral can cancel, an individual order can come out *above* its
perfectly-reflecting value — an amplitude grating diffracts into directions
where a pure phase grating cancels. The summed efficiency cannot, and that is
the conservation statement worth testing.

### `"average"` — the groove-cycle mean

$\langle R(\zeta(t))\rangle$ over the **whole** period, shadowed parts counted as
zero, applied as one factor per wavelength. It sees the shadowing and the
varying local angle but not the interference between them, so it stays
order-independent. The useful halfway diagnostic: the gap between `average` and
`local` *is* the amplitude-grating term.

**It breaks reciprocity too**, by up to $1.5\times10^{-2}$, and for exactly the
same reason `facet` does: $\zeta(t)$ is built from $\alpha$ alone. Resolving the
groove is not what repairs reciprocity — *symmetrising in the exit direction*
is, and only `local` does that. Use `average` to see how much of the change is
shadowing rather than interference, not as a physical model in its own right.

### `"facet"` — the M15 model, kept for reproducibility

One $R$ per wavelength at the active-facet angle
($\sin\zeta = \sin\gamma\cos(\delta-\alpha)$), or at the mean surface for a
profile with no single facet angle. **It breaks reciprocity**, by up to
$3\times10^{-2}$, because that $\zeta$ is built from $\alpha$ alone — a defect
that went unnoticed because `check_reciprocity` had only ever been pointed at
uncoated problems. Kept so an earlier run can be reproduced exactly and so the
size of the change is measurable rather than asserted.

### Visibility: which shadows the masks see

Both resolved models zero the reflection where the groove is shadowed, and
`visibility` selects what counts as shadowed.

**`"facet-normal"` (default)** is the local orientation test alone:
$\sin\zeta \leq 0$ means the point's own facet is turned away from the
direction in question. It cannot see a **cast** shadow — the groove apex
blocking surface beyond the trough that faces the ray perfectly well.

**`"horizon"`** adds them, via
[`geometry.horizon_visible`](../../src/gratinglab/geometry.py): a running-maximum
horizon scan in the transverse plane, run once for the incident direction and
once per diffracted order. Three pieces of geometry make it simpler than it
looks:

- **The half-cone angle drops out.** The surface is invariant along the
  grooves, so a 3D ray occludes exactly as its transverse projection does, and
  the scan runs at the azimuth alone — $\alpha$ in, $\beta_m$ out.
- **Reciprocity survives by construction.** Occlusion along a straight ray
  reads the same from either end, so the incident-at-$\theta$ and
  exit-at-$\theta$ masks are the same function and the pair is symmetric under
  $\alpha \leftrightarrow \beta_m$ — measured at $10^{-16}$, coated and
  uncoated.
- **Self-shadowing is a special case.** A back-facing point has a falling
  ray-adapted height and sits under the running horizon, so the horizon mask
  subsumes the orientation test.

On an ideal sawtooth the cast shadow past the trough has the closed form
$\Delta = (1-t_a)(s_a - s_r)/(s_b + s_r)$ with $s_r = \cot\theta$, zero until
the ray dips below the anti-blaze slope ($\theta > 19.5°$ for the reference
groove) — `tests/test_geometry.py` derives it independently and pins the scan
to it. What it is worth on the reference geometry splits sharply by side:

- **Incident** ($\alpha = 19.99°$): +0.38% of the period beyond the 16.7%
  anti-blaze facet the orientation test already sees. Nearly nothing.
- **Exit**: 9.5–52% of the period per order — $\beta$ is much steeper than
  $\alpha$ for the working orders. The blaze order at its blaze wavelength
  ($m=+3$, $\lambda = 2.226$ nm, $\beta = 39°$) has 14.7% of the period
  extra-masked and moves **0.514 → 0.348, −32%** — the same order of magnitude
  as the M16 groove-resolution change, invisible to the orientation test.

With **no coating**, `visibility="horizon"` applies the masks at unit
amplitude: geometric shadowing on a perfect reflector, which is the
configuration a perfect-conductor cross-check against the integral solver
wants on grooves deep enough to shadow themselves. The default uncoated run
stays the untouched pure phase integral. `"average"` sees the horizon on the
incident side only (its exit side is order-blind by construction), and
`"facet"` — which has no masks — refuses the combination rather than ignoring
it.

Two caveats, stated rather than implied. The mask is **ray optics**: near
passing-off it overestimates blocking, because diffraction bends around the
apex; the vanishing flux obliquity already suppresses those orders. And like
the geometric mean, it is **not yet validated against a rigorous
finite-conductivity method** — it moved the blaze order from 0.69 toward the
rigorous 0.21 of the (out-of-regime) perfect-conductor comparison, which is
suggestive rather than proof. It ships opt-in so existing results reproduce
bit-for-bit; the default flips once a finite-conductivity A/B lands in
[`findings.md`](../findings.md), the same road `facet` → `local` took.

### What it is worth

Reference geometry — Blazed 29.5°/70.5°, $\gamma=1.25°$, $\alpha=19.99°$, Au,
$\lambda = 3$ nm. The anti-blaze facet is **16.7% of the period and entirely
self-shadowed**, and `facet` applies the active facet's reflectivity across all
of it:

| $m$ | bare | `facet`/bare | `average`/bare | `local`/bare |
|---:|---:|---:|---:|---:|
| $-1$ | 0.00473 | 0.7146 | 0.5954 | 0.4790 |
| $0$ | 0.01454 | 0.7146 | 0.5954 | 0.5327 |
| $+1$ | 0.02750 | 0.7146 | 0.5954 | 0.5807 |
| $+2$ | 0.07385 | 0.7146 | 0.5954 | 0.6448 |
| $+3$ | 0.41571 | 0.7146 | 0.5954 | 0.7581 |

The `facet` column is flat by construction. The `local` column is not, and the
blaze order moves least — which is what "collapses to $r(\zeta)$ at specular"
looks like from the outside.

**This is not yet validated against a rigorous method.** The geometric mean is a
symmetrisation chosen because it preserves the invariant this solver is built
around, not because a vector calculation confirmed it. The corpus A/B in
[`findings.md`](../findings.md) is how that gets judged.

**Which polarization.** Unpolarized, always, whatever the illumination says.
The $s$ and $p$ of a Fresnel reflection are defined against the *facet's* plane
of incidence; `Illumination.polarization` is TE/TM defined against the
**grooves** (`conventions.md` §7). In a conical mount those frames differ, and
the rotation between them is deliberately unowned
([`fresnel.py`](../../src/gratinglab/materials/fresnel.py) refuses TE/TM labels
for exactly this reason). Resolving it here would be false precision on a model
that already reports TE and TM as identical — and the size of what is being
waived is measured rather than assumed. Fractional intensity split
$|R_s - R_p|/R_s$ on the vendored Au table:

| $\lambda$ (nm) | $\zeta = 0.5°$ | $1.25°$ | $2°$ | $5°$ |
|---:|---:|---:|---:|---:|
| 0.62 | $1.2\times10^{-4}$ | $3.9\times10^{-4}$ | $2.2\times10^{-3}$ | $2.8\times10^{-2}$ |
| 3.0 | $1.9\times10^{-3}$ | $4.9\times10^{-3}$ | $7.9\times10^{-3}$ | $2.4\times10^{-2}$ |
| 6.0 | $2.2\times10^{-3}$ | $5.5\times10^{-3}$ | $8.9\times10^{-3}$ | $2.5\times10^{-2}$ |

Since the unpolarized mean is already computed, any basis rotation moves the
result by at most **half the split** — under 0.5% at every graze the reference
groove actually presents (its worst local graze is 1.94°). The provenance
checks the split at the guarded graze and warns above 2% (§7), so the waiver is
priced exactly when it stops being free. Under `"local"` the two are combined
at *amplitude* level per order, which is the more defensible place.

**Roughness.** Névot–Croce by default, Debye–Waller and `none` on request, and
both are **intensity** factors — M16-B fixed a Névot–Croce that returned the
amplitude form and was being multiplied into a reflectivity. The limit that
pins it is $n \to 1$, where $k_{tz} \to k_{iz}$ and the two models must become
the same expression. Near and above $\theta_c$ Debye–Waller over-damps by
$\sim10^{-2}$; deep below it the ordering reverses and Debye–Waller retains
marginally more; far above $\theta_c$ they converge to $10^{-5}$.

**Energy.** Summed efficiency falls for three distinct reasons now: the
thin-element defect of §5, which is the model straying; ordinary absorption,
which is physics; and geometric shadowing — self-shadowing always, cast
shadowing under `"horizon"` — which is also physics. The provenance separates
the first from the other two.

---

## 10. Finite N and resolving power

Everything above treats the grating as infinite: each order is a delta
function in angle, which is the `N → ∞` limit of ISSI eq. (8),

```
[sin(Ns) / (N sin s)]²,    s ≡ [sin α + sin β]·sin γ·period·π/λ
```

The factor is deliberately **not** applied to order efficiencies — at an exact
order `s = mπ` it equals 1, so multiplying by it changes nothing. Its content
is the *angular line shape* around each order, and that is a different
question, answered by `gratinglab.resolution` rather than by this solver:

- `resolving_power` returns the Rayleigh closed form `R = |m|·N`, which follows
  from the first zeros at `s = mπ ± π/N` in two lines (`conventions.md` §9).
  It is mount-independent, because `sin γ` scales the line width and the
  dispersion alike.
- `line_profile` evaluates the interference function around `β_m`, and the
  closed form and the profile test each other — the same pattern as §4's
  closed-form efficiency anchors.

`N` is the **illuminated** groove count, `Problem.n_grooves`. A problem built
from a measured boundary carries the count the scan actually averaged
(`BoundaryProfile.to_problem`), which is how a measured surface flows through
to a spectrograph number in-process. A problem without `n_grooves` is refused,
not treated as infinite: "R is undefined" and "R is very large" are different
answers.

**Scope.** This is the ideal-grating result: `N` identical grooves at exact
spacing. Resolving power degraded by groove-placement and period errors —
computed from a measured groove ensemble rather than assumed away — is the
named successor in `roadmap.md`, and nothing here approximates it.

---

## References

- McCoy, *Scalar Treatment of Gratings*, PhD thesis Appendix D — the primary
  derivation. See `conventions.md` §10 for the errata this page corrects.
- Heilmann, Huenemoerder, McCoy & McEntaffer, *Diffraction Gratings for X-ray
  Spectroscopy*, Springer ISSI Scientific Reports (2024),
  [arXiv:2409.02297](https://arxiv.org/abs/2409.02297) — §2.1 conical framework,
  §4 blazed-RG scalar DE. **Normative where it disagrees with the thesis.**
- Harvey & Pfisterer, *Understanding diffraction grating behavior* I & II,
  Opt. Eng. **58**(8) 087105 (2019), **59**(1) 017103 (2020) — parametric
  treatment handling conical diffraction natively.
- Born & Wolf, *Principles of Optics* — Kirchhoff theory and the Fraunhofer limit.
