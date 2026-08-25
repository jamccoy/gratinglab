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
smoothed facet. Two boundary conditions are offered through the
`conductivity` option:

| Mode | Physics | Normalisation |
|---|---|---|
| `"perfect"` (default) | The grating is a **perfect conductor**: the field does not penetrate, there is no absorption, and the sum over propagating orders is **exactly 1** — a theorem (§5). A named `coating` is ignored, with a provenance warning. | Relative to a perfect reflector — matching the PC PCGrate reference tables. |
| `"tabulated"` | One interface into a **semi-infinite material** whose complex index is read per wavelength from the named `coating` (falling back to `substrate`) via `materials.lookup` — the coupled conical system of Goray & Schmidt (2010), §8 below. Absorption is computed and recorded, and `R + A = 1` is the theorem (§5). | **Absolute.** |

Common to both: one boundary, one period, classical periodicity — multilayer
stacks and crossed gratings are out of scope. Asking for anything else (an
impedance/Leontovich rung, say) raises `UnsupportedConfiguration`: that
approximation requires |n| ≫ 1 and is invalid in the soft-X-ray regime this
solver targets. Folding a Fresnel reflectivity onto the perfect-conductivity
result is **not** equivalent to the tabulated mode — that is the
transcription error catalogued as `conventions.md` §10 item 5, and §8's
validation now measures the difference instead of assuming it.

Everything else that limits other methods is *in* scope: conical (off-plane)
mounts are exact (§2), vertical facets are genuine geometry, and undercut
boundaries — representable by no height function — solve the same way,
because the method parametrises the curve itself
(`Capabilities.handles_undercut = True`, alone in the registry).

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
the conical problem is a coupled 2×2 system — implemented as §8, where the
same reduced wavelength still governs the vacuum side.

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
eqs. 4.34/4.42; thesis `eq:prop_order_unity`, by Green's theorem). With
finite conductivity the theorem becomes `R + A = 1` (G&S §2.C) — and it
stays a *two-sided* check, because the absorption `A` is computed by an
independent boundary integral of the solved densities (§8), not as
`1 − R`; for a lossless substrate the same integral equals the transmitted
power and `R + T = 1` is the theorem. Either way this is the sharpest
self-check a solve can carry, and the solver **reports, never rescales**:
the deviation of the computed conserved sum from 1 is the discretisation
error estimate, recorded in `provenance.notes["energy_balance_deviation"]`
and warned about above 0.5%. `checks.check_energy_balance` consumes the
scan's recorded absorption automatically.

Measured behaviour, honestly stated:

- **TE holds it to ~1e-5 or better on every profile tried**, corners
  included, at moderate `boundary_points`.
- **TM on smooth profiles is comparable.** On profiles with corners (ideal
  sawtooth, lamellar) the Meixner edge behaviour meets a plain equal-arc mesh
  and TM converges at **first order** — a percent-level deficit at P ≈ 100,
  still ~0.5% at P ≈ 1000. The convergence harness tells the truth about it;
  a graded corner mesh (M&P §4.6.6) is the known follow-up.
- **Finite conductivity holds `R + A = 1` to ~1e-4 on smooth profiles**
  (measured: 1.2e-4 for Au at 2 nm, 1.5° graze, P = 256; 1e-5 on a flat
  interface). Corners now limit *both* polarizations at first order — the
  edge-singular densities meet the same equal-arc mesh — reported through
  the same deviation.

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

Finite conductivity adds its own ladder (`tests/test_integral_finite.py`):

- **Flat interface = Fresnel, amplitude and phase**: TE against `r_s`, TM
  against `r_p` from `materials/fresnel.py`, to 1e-4 in-plane; in the
  conical mount the full 2×2 (E_z, B_z) reflection matrix — polarization
  conversion included — matches the rotated-Fresnel closed form to 6e-5.
  This is the successor of the flat-mirror sign pins: every sign in the
  coupled system rests on it.
- **Energy theorems**: `R + T = 1` (glass sinusoid) and `R + A = 1` (gold
  sinusoid) to a few 1e-4; the boundary-integral `A` equals the transmitted
  Rayleigh sum for a lossless substrate — two unrelated post-processings of
  one solve.
- **Perfectly conducting limit**: at n = 100i the solve lands on the
  validated PC solver to 1.3e-3 worst order — the residual is the ~1 nm
  field penetration the mode is deliberately adding, not error — with
  absorption at numerical zero; along a lossy path (n = |n|·e^{iπ/4}) the
  absorbed fraction falls monotonically toward the conductor.
- **Goray & Schmidt Table 3** (dielectric sinusoid, conical, pure E_z
  incidence; itself a cross-method comparison against Li's CM): every
  tabulated order, reflected and transmitted, reproduced within 6e-5.
  Their Table 4 turned out to be a publisher's duplication of Table 3
  (`docs/findings.md`).
- **Lorentz reciprocity survives absorption**: measured ~5e-5 for Au at
  X-ray grazing incidence through `checks.check_reciprocity`, unchanged.

## 7. Choosing `boundary_points`

The solver refuses fewer than ~6 nodes per reduced wavelength along the
boundary arc — below that the field is unresolved and the answer would be
confidently wrong. The tabulated mode applies the same floor to the
**metal-side** transverse wavelength `λ / |√(n² − cos²γ)|`, which is the
shorter one whenever the mount sits below the critical angle (the field
then decays into the metal on the `λ/√(2δ)` scale — for Au at 2 nm and
1.5° graze that is ~21 nm against a 76 nm vacuum reduced wavelength, a
floor of ~210 points on a 600 nm sinusoid) and dramatically shorter for
visible-light metals (|n| of a few), which is why that regime costs nodes
here (M&P §4.6.4). Above the floors:

- Smooth profiles: convergence is fast; the harness typically certifies
  P in the low hundreds.
- Cornered profiles + TM: expect the energy-balance deviation to shrink only
  linearly; let `check_convergence` (which sweeps this knob) pick the plateau
  rather than trusting a single run.
- Cost: the kernels are assembled as separable BLAS products
  ([`_kernels.py`](../../src/gratinglab/solvers/integral/_kernels.py)) — a
  few `(P, 2M+1)·(2M+1, P)` GEMMs per wavelength plus a one-time O(P²·M)
  geometry pass (the cached Kummer tails) amortised over the scan. Measured
  at P = 400: ~100 ms per wavelength unpolarized, so the default 400 at 200
  wavelengths is ~20 s. (The kernels were previously evaluated pointwise
  over `(P, P, M)` broadcasts at ~14 s per wavelength — and contrary to
  what this document used to claim, the O(P³) dense LU was never the
  bottleneck: it is milliseconds at any realistic P, ~0.1% of the old
  runtime.) The GUI runs the solve on the worker thread with a live
  progress bar and a Cancel that actually cancels. The tabulated mode
  costs ~2× the perfect one per wavelength (more kernels and operator
  products, but the double/adjoint/tangential layers share one kernel
  evaluation per medium); both incident polarizations share one
  factorization, so `"unpolarized"` costs two triangular solves, not two
  factorizations.
- Formulation envelope: the Kummer acceleration cancels asymptote terms of
  magnitude `exp(|α₀|·|Z|)`, so for `|α₀|·h ≳ 30` (`h` the groove height —
  reachable only by driving a short *transverse* wavelength through a deep
  groove, e.g. an in-plane mount at hard-X-ray wavelengths on a large
  period) double precision has nothing left and the kernels are garbage in
  any formulation of this acceleration, old or new — measured to affect
  both equally. The grazing-conical mounts this solver targets keep
  `|α₀| = k sin γ · sin α` small precisely because `sin γ` is small; if the
  regime is ever hit, the energy-balance theorem check fails loudly (R+A
  deviations of order 1, reported in provenance) rather than silently.

## 8. Finite conductivity: the coupled conical system

Implemented in [`_finite.py`](../../src/gratinglab/solvers/integral/_finite.py),
following Goray & Schmidt (2010) — but **derived independently in this
project's conventions** rather than transcribed: G&S's normal points into
the metal and their potentials carry a factor 2, and sign-sensitive
equations do not survive PDF extraction (`references.md`). The
construction, in our operators (outward normal up, `_nystrom` scalings):

- The two coupled scalars are `E_z` and `B_z = (μ⁺/ε⁺)^{1/2} H_z`, each
  satisfying a transverse Helmholtz equation with its own side's
  wavenumber: `k_t⁺ = k sin γ` (exactly the reduced-wavelength problem of
  §2) and `k_t⁻ = k √(n² − cos²γ)` (the same Green's-function code with
  complex k — the branch pin in `_greens` was built for this).
- The substrate fields are single-layer potentials, `u⁻ = S⁻w`,
  `v⁻ = S⁻τ`: traces `V⁻w`, normal derivatives `(L⁻ − I/2)w` by the jump
  relation the `_nystrom` tests pin.
- On the vacuum side, Green's representation of the radiating scattered
  field plus the regularity identity of the incident field give the exact
  boundary relation `(I/2 + K⁺)ψ − V⁺ ∂_n ψ = ψ_i` for the **total**
  field — whose perfect-conductor limits are literally the two milestone-1
  equations (`∂_n ψ = 0` gives TM, `ψ = 0` gives TE), which is what pins
  its signs.
- The jump conditions (G&S eq. 6, mapped) eliminate the vacuum traces:
  `E_z`, `B_z` continuous, and

      ∂_n E_z⁺ = c_E ∂_n E_z⁻ + s ∂_t B_z,
      ∂_n B_z⁺ = c_B ∂_n B_z⁻ − s ∂_t E_z,

  with `c_E = n² k_t⁺²/k_t⁻²`, `c_B = k_t⁺²/k_t⁻²`, and
  `s = cos γ (1 − k_t⁺²/k_t⁻²)`. The tangential derivatives act on the
  *continuous* traces, so they are `D_t V⁻` on the densities — a
  geometry-only spectral differentiation matrix, no new singular
  quadrature (G&S §3 note the same option).

The result is a 2×2 block system in `(w, τ)`: diagonal blocks
`(I/2 + K⁺)V⁻ − c V⁺(L⁻ − I/2)`, cross blocks `∓ s V⁺ D_t V⁻`. In-plane
(`cos γ = 0`, exact — computed from `sin γ` so the cross blocks vanish
identically) the system decouples into the classical TE/TM transmission
problems. In the perfect-conductor limit `|n| → ∞` it degenerates to the
milestone-1 equations (`V⁻` and `c_B` die like 1/|n|).

Per reflected order the efficiency is `(γ_m/γ₀)(|E_m|² + |B_m|²)` for a
unit-power incident state — both components, because conical incidence
converts polarization (TE ≡ incident `(p_z, q_z) = (1, 0)`, TM ≡ `(0, 1)`
in the (E_z, B_z) basis; `conventions.md` §7). Absorption is the
boundary integral of G&S eq. 26 mapped to our operators — on a flat
interface it reduces *exactly* to `1 − |r|²` in both polarizations, which
is how its sign and normalisation were pinned analytically before any
test ran.

- D. Maystre and E. Popov, "Integral Method for Gratings", ch. 4 of
  E. Popov (ed.), *Gratings: Theory and Numeric Applications*, AMU (2012).
- J. A. McCoy, PhD thesis, Chapter 2 §2.2 (read with `conventions.md` §10;
  §2.2.3 is superseded — see item 5 there).
- L. I. Goray and G. Schmidt, JOSA A 27, 585 (2010) — the
  finite-conductivity formulation (§8), implemented at M19.
- L. Li, JOSA A 10, 2581 (1993) — the coordinate-transformation results
  behind G&S's comparison tables, the external anchor of §6.
- H. A. Kalhor and A. R. Neureuther, JOSA 61, 43 (1971) — early
  perfectly-conducting integral treatment.
