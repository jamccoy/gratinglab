# Integral method (boundary integral equation)

Implemented in [`solvers/integral/`](../../src/gratinglab/solvers/integral/__init__.py).
Symbols and sign conventions follow [`conventions.md`](../conventions.md);
"M&P" below is Maystre & Popov, *Integral Method for Gratings*, ch. 4 of
Popov (ed.), *Gratings: Theory and Numeric Applications* (2012) — the
numerical recipe this implementation follows. The thesis derivation
(Chapter 2 §2.2) is the notation and sign reference, read together with the
corrections catalogued in [`conventions.md` §10](../conventions.md).

---

## 1. What the model assumes

This is a **rigorous** method: Maxwell's equations are solved exactly on the
true boundary curve, with no thin-element approximation, no staircase, no
smoothed facet. The milestone-1 scope makes exactly one physical idealisation:

| Assumption | Consequence |
|---|---|
| The grating is a **perfect conductor** | The field does not penetrate the surface; there is no absorption. Efficiencies are *relative to a perfect reflector* and the sum over propagating orders is **exactly 1** — a theorem (§5), not an aspiration. A named `coating` is ignored, with a provenance warning. |
| One boundary, one period, classical periodicity | Multilayer stacks and crossed gratings are out of scope. |

Everything else that limits other methods is *in* scope: conical (off-plane)
mounts are exact (§2), vertical facets are genuine geometry, and undercut
boundaries — representable by no height function — solve the same way,
because the method parametrises the curve itself
(`Capabilities.handles_undercut = True`, alone in the registry).

Finite conductivity — absorption, and with it absolute soft-X-ray
efficiencies — is the planned successor milestone, following Goray & Schmidt,
*Solving conical diffraction grating problems with integral equations*,
JOSA A 27 (2010). The `conductivity` option is the reserved seam; asking for
anything but `"perfect"` today raises `UnsupportedConfiguration` rather than
approximating. Folding a Fresnel reflectivity onto the perfect-conductivity
result is **not** a substitute — that is the transcription error catalogued
as `conventions.md` §10 item 5.

## 2. The conical reduction

For a perfect conductor, the off-plane problem at cone half-angle γ decouples
exactly into an in-plane problem at the **reduced wavelength**

    λ̄ = λ / sin γ

(M&P eq. 4.65; thesis Chapter-2.tex:931). Every field component separates as
`f(x, y) · exp(i k₀ cos γ · z)`, the conserved `k_z` drops out of the
transverse Helmholtz equation, and on a perfectly conducting boundary (whose
normal has no z-component) the two scalars never couple:

- **TE** (E ∥ grooves): `E_z` satisfies a **Dirichlet** condition — the total
  field vanishes on the boundary.
- **TM** (H ∥ grooves): `H_z` satisfies a **Neumann** condition — its normal
  derivative vanishes.
- Cross-polarization conversion is identically zero, and `"unpolarized"` is
  the plain average of two independent solves.

This decoupling is a perfect-conductivity privilege. With finite conductivity
the conical problem is a coupled 2×2 system — one reason that milestone is
separate.

The reduction is also what makes the X-ray case *cheap*. Discretisation cost
follows λ̄, not λ: at γ = 1.5°, a 1 nm photon becomes a 38 nm transverse
wavelength against a ~315 nm period, so a few hundred boundary points
resolve the field. The extreme off-plane mount — the reason this project
exists — is the integral method's easiest regime, not its hardest.

## 3. The boundary integral equations

Work in the transverse plane, period `d`, `K = 2π/d`, transverse wavenumber
`k = 2π/λ̄`, incident x-wavenumber `α₀ = −k sin α`, and
`α_m = α₀ + mK`, `γ_m = √(k² − α_m²)` with `Im γ_m ≥ 0` (so `α_m = k sin β_m`
reproduces `geometry.sin_beta` exactly). The quasi-periodic Green's function
is the Rayleigh sum

    G(X, Z) = (1 / 2id) Σ_m (1/γ_m) exp(i α_m X + i γ_m |Z|).

**TE** uses a single-layer ansatz: an unknown surface current σ on the
boundary with `∫ G σ ds' = −ψ_incident` on the curve (a first-kind equation,
M&P eq. 4.31). **TM** solves for the total boundary field ψ through the
double-layer jump relation, `ψ/2 + ∫ (∂G/∂n′) ψ ds' = ψ_incident` (the
second-kind equation; M&P eq. 4.39 states the equivalent in its own sign
convention). Rayleigh amplitudes then come from one quadrature per order
(M&P eqs. 4.33, 4.41) and efficiencies are `(γ_m/γ₀)|r_m|²`.

## 4. The numerics, and where the accuracy comes from

Everything below is [`_greens.py`](../../src/gratinglab/solvers/integral/_greens.py)
and [`_nystrom.py`](../../src/gratinglab/solvers/integral/_nystrom.py);
equation anchors are M&P §4.6.

- **Equal-arc-length Nyström mesh.** The rectangular rule on equally spaced
  nodes of a periodic integrand is spectrally accurate (M&P eqs. 4.70–4.74),
  which is why `Profile.boundary()` resamples every profile at equal arc
  length — and why `Blazed`/`Lamellar` supply exact vertex polygons, vertical
  segments included.
- **Kummer acceleration.** The raw Green's sum converges like 1/m. The
  large-|m| asymptote sums in closed form to two logarithms
  (eqs. 4.82–4.85); the paired remainder decays like m⁻³ on the diagonal
  (eq. 4.87) and at worst m⁻² on the boundary. `spectral_terms` defaults to
  `boundary_points / 2` so the single public knob converges both
  discretisations together (M&P's own Table 4.1 converges at P ≈ 2.2M).
- **Log singularity.** The kernel's integrable `ln` singularity is handled by
  the periodic comparison-function split (eqs. 4.97–4.100): the singular part
  integrates analytically to `−L/π`, and what the quadrature sees is
  continuous.
- **Neumann diagonal.** The double-layer kernel is continuous on a smooth
  curve; its diagonal (which carries the local curvature) is taken as the
  numerical principal-value limit — the average of the kernel at the two
  physically adjacent nodes — measured to converge as O(spacing²).

## 5. Energy balance is a theorem here

For a perfect conductor, `Σ_m ρ_m = 1` over propagating orders exactly (M&P
eqs. 4.34/4.42; thesis `eq:prop_order_unity`, by Green's theorem). This is
the sharpest self-check a solve can carry, and the solver **reports, never
rescales**: the deviation of the computed sum from 1 is the discretisation
error estimate, recorded in `provenance.notes["energy_balance_deviation"]`
and warned about above 0.5%.

Measured behaviour, honestly stated:

- **TE holds it to ~1e-5 or better on every profile tried**, corners
  included, at moderate `boundary_points`.
- **TM on smooth profiles is comparable.** On profiles with corners (ideal
  sawtooth, lamellar) the Meixner edge behaviour meets a plain equal-arc mesh
  and TM converges at **first order** — a percent-level deficit at P ≈ 100,
  still ~0.5% at P ≈ 1000. The convergence harness tells the truth about it;
  a graded corner mesh (M&P §4.6.6) is the known follow-up.

## 6. Validation evidence

- **Flat mirror**: unit specular efficiency, and the amplitude *signs* pin
  the boundary conditions (TE r₀ = −1, TM r₀ = +1).
- **M&P Table 4.1** (sinusoid, Littrow, their own P = 110, M = 50): TE
  0.4659/0.5341 and TM 0.9579/0.0421 reproduced to 1e-3
  (`tests/test_integral_core.py`).
- **Maréchal–Stroke theorem** (echelette, TM, exact blaze): approached
  first-order through the corner limitation, monotonically.
- **Lorentz reciprocity** (`checks.check_reciprocity`): TE to ~1e-8, TM to
  ~1e-5.
- **PCGrate, orderwise, on a measured groove**: the TASTE
  perfect-conductivity wavescan is reproduced to **7.8e-4 across every
  propagating order of the full 0.6–6.0 nm scan** — and the comparison
  *identified* the previously unconfirmed profile behind that run
  (`docs/findings.md`; `tests/test_integral_corpus.py`).

## 7. Choosing `boundary_points`

The solver refuses fewer than ~6 nodes per reduced wavelength along the
boundary arc — below that the field is unresolved and the answer would be
confidently wrong. Above the floor:

- Smooth profiles: convergence is fast; the harness typically certifies
  P in the low hundreds.
- Cornered profiles + TM: expect the energy-balance deviation to shrink only
  linearly; let `check_convergence` (which sweeps this knob) pick the plateau
  rather than trusting a single run.
- Cost is O(P²·M) to build the kernels and O(P³) to solve — roughly cubic in
  the knob. The default 400 at 200 wavelengths is minutes, not seconds; the
  GUI runs it on the worker thread with a live progress bar and a Cancel
  that actually cancels.

## References

- D. Maystre and E. Popov, "Integral Method for Gratings", ch. 4 of
  E. Popov (ed.), *Gratings: Theory and Numeric Applications*, AMU (2012).
- J. A. McCoy, PhD thesis, Chapter 2 §2.2 (read with `conventions.md` §10;
  §2.2.3 is superseded — see item 5 there).
- L. I. Goray and G. Schmidt, JOSA A 27, 585 (2010) — the finite-conductivity
  successor.
- H. A. Kalhor and A. R. Neureuther, JOSA 61, 43 (1971) — early
  perfectly-conducting integral treatment.
