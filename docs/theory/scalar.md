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
| Structure is large compared to the wavelength, λ ≪ p | Sub-wavelength features are misrepresented; the model degrades smoothly rather than failing loudly. |
| Kirchhoff (thin-element) boundary condition: the field just above the surface is the incident field times a local phase | No multiple scattering between groove facets, no shadowing, no field penetration into the material. |
| Fraunhofer (far-field) observation | Valid for a telescope focal plane; not for near-field work. |
| Perfect reflectivity, or reflectivity applied separately as a scale factor | Efficiencies are **relative**, not absolute, until a materials layer supplies R_F. |

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
\boxed{\;\mathscr{E}_m = \left|G_m\right|^2\;}$$

Two implementation notes follow directly:

- Because $\Phi_m$ carries an $m$, this is **one integral per order**, not a
  single FFT. The cost is negligible and the alternative would be wrong.
- The integrand is periodic, so the rectangle rule converges *spectrally*. This
  is why `quadrature_points` reaches machine precision at modest values, and why
  it must still exceed $2|m|_{\max}$ to satisfy Nyquist.

### Two factors that are *not* here, and why

Thesis Appendix D carries an obliquity factor,

$$\mathscr{E}_m \propto \frac{\cos\beta_m}{\cos\alpha}\left|G_m\right|^2$$

from projecting the Poynting flux onto a plane parallel to the surface, and then
renormalises so that $\sum_m \mathscr{E}_m = 1$. **Both are removed.** For a
grating with $N \to \infty$ grooves the efficiency is the norm-squared Fourier
coefficient and nothing else — ISSI §2.1 is explicit about this, and its groove
function `f(x) = exp{ik[n(k)−1]y(x/p)}` carries no extra factor.

The renormalisation is the more damaging of the two: forcing the sum to unity
would erase the model's own signal about how far it has strayed from its regime
(§5). Neither is offered as an option, because neither is an alternative
convention — they are errors.

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
| Energy, $\sum_m \mathscr{E}_m \leq 1$ | ✗ violated, up to 1.71 | ✓ exact (Parseval) |
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
| 0.0010 | 0.99999 |
| 0.0200 | 0.99592 |
| 0.0500 | 0.97483 |
| 0.1000 | 0.90472 |

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
| $\lambda / p \lesssim 0.1$ | Kirchhoff theory assumes structure ≫ wavelength. Soft X-ray work runs ~0.005; visible gratings ~0.4. |
| $\zeta < \theta_c \approx \sqrt{2\delta_n}$ | Facet graze must stay below the critical angle for total external reflection, or the reflectivity assumption collapses. Needs a materials layer to evaluate. |
| $\lambda > 32 \sin(\zeta)\,\sigma$ | Fraunhofer smoothness: below this the facet is not optically smooth and scatter dominates. |
| $\sum_m \mathscr{E}_m \leq 1$ | Not a validity condition so much as a sanity bound; see §5. |

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
