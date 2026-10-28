# Integral method (boundary integral equation)

Implemented in [`solvers/integral/`](../../src/gratinglab/solvers/integral/__init__.py).
Symbols and sign conventions follow [`conventions.md`](../conventions.md);
"M&P" below refers to Maystre & Popov, *Integral Method for Gratings*, ch. 4 of
Popov (ed.), *Gratings: Theory and Numeric Applications* (2012), whose numerical
recipe this implementation follows. The thesis derivation (Chapter 2 §2.2) is the
reference for notation and signs, read together with the corrections catalogued
in [`conventions.md` §10](../conventions.md).

---

## 1. What the model assumes

This is a rigorous method: Maxwell's equations are solved exactly on the true
boundary curve, with no thin-element approximation, no staircase and no smoothed
facet. Two boundary conditions are offered through the `conductivity` option:

| Mode | Physics | Normalisation |
|---|---|---|
| `"perfect"` (default) | The grating is a perfect conductor: the field does not penetrate, there is no absorption, and the sum over propagating orders is exactly 1, which is a theorem (§5). A named `coating` is ignored, with a provenance warning. | Relative to a perfect reflector, matching the PC PCGrate reference tables. |
| `"tabulated"` | One interface into a semi-infinite material whose complex index is read per wavelength from the named `coating` (falling back to `substrate`) via `materials.lookup`, using the coupled conical system of Goray & Schmidt (2010) described in §8 below. Absorption is computed and recorded, and `R + A = 1` is the theorem (§5). Into a lossless material the transmitted orders are recorded as the scan's `transmission` (since 2026-10-27), and the recorded absorption is the boundary integral minus their total, so the scan's conserved sum is `R + T + A`. | **Absolute.** |

Both modes assume one boundary, one period and classical periodicity; multilayer
stacks and crossed gratings are out of scope. A request for anything else (an
impedance/Leontovich rung, for example) raises `UnsupportedConfiguration`, since
that approximation requires |n| ≫ 1 and is invalid in the soft-X-ray regime this
solver targets. Folding a Fresnel reflectivity onto the perfect-conductivity
result is not equivalent to the tabulated mode. That shortcut is the
transcription error catalogued as `conventions.md` §10 item 5, and the validation
in §8 now measures the difference instead of assuming it.

Everything else that limits other methods is within scope. Conical (off-plane)
mounts are exact (§2), vertical facets are genuine geometry, and undercut
boundaries, which no height function can represent, are solved in the same way,
because the method parametrises the curve itself
(`Capabilities.handles_undercut = True`, the only such solver in the registry).

## 2. The conical reduction

For a perfect conductor, the off-plane problem at cone half-angle γ decouples
exactly into an in-plane problem at the *reduced wavelength*

    λ̄ = λ / sin γ

(M&P eq. 4.65; thesis Chapter-2.tex:931). Every field component separates as
`f(x, y) · exp(i k₀ cos γ · z)`, the conserved `k_z` drops out of the transverse
Helmholtz equation, and on a perfectly conducting boundary (whose normal has no
z-component) the two scalars never couple:

- **TE** (E ∥ grooves): `E_z` satisfies a Dirichlet condition, so the total
  field vanishes on the boundary.
- **TM** (H ∥ grooves): `H_z` satisfies a Neumann condition, so its normal
  derivative vanishes.
- Cross-polarization conversion is identically zero, and `"unpolarized"` is the
  plain average of two independent solves.

This decoupling holds only for perfect conductivity. With finite conductivity
the conical problem is a coupled 2×2 system, implemented as described in §8, in
which the same reduced wavelength still governs the vacuum side.

The reduction is also what makes the X-ray case inexpensive. The discretisation
cost follows λ̄ and not λ: at γ = 1.5°, a 1 nm photon becomes a 38 nm transverse
wavelength against a ~315 nm period, so a few hundred boundary points resolve the
field. The extreme off-plane mount, which motivates this project, is therefore
the easiest regime for the integral method.

## 3. The boundary integral equations

The equations are posed in the transverse plane, with period `d`, `K = 2π/d`,
transverse wavenumber `k = 2π/λ̄`, incident x-wavenumber `α₀ = −k sin α`, and
`α_m = α₀ + mK`, `γ_m = √(k² − α_m²)` with `Im γ_m ≥ 0` (so that
`α_m = k sin β_m` reproduces `geometry.sin_beta` exactly). The quasi-periodic
Green's function is the Rayleigh sum

    G(X, Z) = (1 / 2id) Σ_m (1/γ_m) exp(i α_m X + i γ_m |Z|).

TE uses a single-layer ansatz: an unknown surface current σ on the boundary with
`∫ G σ ds' = −ψ_incident` on the curve (a first-kind equation, M&P eq. 4.31). TM
solves for the total boundary field ψ through the double-layer jump relation,
`ψ/2 + ∫ (∂G/∂n′) ψ ds' = ψ_incident` (the second-kind equation; M&P eq. 4.39
states the equivalent in its own sign convention). The Rayleigh amplitudes then
follow from one quadrature per order (M&P eqs. 4.33, 4.41), and the efficiencies
are `(γ_m/γ₀)|r_m|²`.

## 4. Numerical method and sources of accuracy

Everything below is implemented in
[`_greens.py`](../../src/gratinglab/solvers/integral/_greens.py) and
[`_nystrom.py`](../../src/gratinglab/solvers/integral/_nystrom.py); the equation
numbers refer to M&P §4.6.

- **Equal-arc-length Nyström mesh.** The rectangular rule on equally spaced
  nodes of a periodic integrand is spectrally accurate (M&P eqs. 4.70–4.74). For
  this reason `Profile.boundary()` resamples every profile at equal arc length,
  and `Blazed`/`Lamellar` supply exact vertex polygons, vertical segments
  included.
- **Corner-graded mesh (`corner_grading=p`, opt-in).** At a corner the densities
  are singular (the Meixner edge condition), and on an equal-arc mesh that
  singularity caps convergence at first order. `Profile.boundary(n, grading=p)`
  instead splits a mesh parameter τ equally, gives each side between two
  corners a share of the nodes in proportion to its length, and maps each share
  through Kress's sigmoidal transformation of order p (M&P §4.6.6; Kress 1990),
  whose derivatives below order p vanish at both ends. The density times
  `ds/dτ` is then smooth in τ and the rectangular rule recovers high order. A
  corner is a polygon vertex turning by more than `CORNER_TURN_DEGREES` (10°),
  no node sits on one, and every node lies on a straight piece with its exact
  normal. A profile with no corner keeps the equal-arc mesh.

  Every quadrature is written in τ, with node weights `w_j = s'(τ_j)/P`, and
  reduces term by term to the equal-arc formulas when `s' = L`. Three places
  needed more than swapping the weights in, and each was wrong at first in a way
  only a solve exposed (`docs/findings.md`):

  - the single layer's logarithmic comparison moves into τ, which changes its
    diagonal to `w_j[R(0,0) + ln(2π s'_j/d)/2π − S̃] − s'_j/π`;
  - the tangential layer's Cauchy comparison acts on `σ s'`, never on σ, because
    a corner density is unbounded in τ;
  - the principal-value diagonals of K, L and D_t V are taken from the neighbour
    on the node's own side of a corner, since the neighbour across it sits where
    the kernel is of order 1/distance.

  The finite-conductivity system needs two more, because it is first kind in
  its densities and so does not absorb an error confined to the nodes next to a
  corner the way the perfectly conducting second-kind TM equation does (§5):

  - **K⁺ is singularity-subtracted at the corners**
    (`_nystrom.corner_subtraction`). For a node a distance a from a corner the
    double-layer kernel over the facet across it varies on the scale a, and the
    graded nodes there are spaced on the same scale, so the rule's error at that
    node is a fixed fraction of the integrand: refinement moves it to a closer
    node instead of shrinking it (up to 0.16 against the identity's 1/2 on a
    lamellar). Writing `K φ = K(φ − φ(s_i)) + φ(s_i) K1` and supplying the
    Laplace part of `K1` exactly, as minus the angle each polygon segment
    subtends over 2π, removes it. `L⁻` carries the same kind of error and needs
    nothing, because it appears only inside `V⁺`, whose weights damp it.
  - **D_t V⁻ is the composition `D_t (V⁻ w)`**, the spectral τ-derivative of the
    assembled trace, instead of the single operator of §8. The trace is smooth in
    τ on this mesh even where w is singular, while the single operator's Cauchy
    remainder is not continuous across a corner. It too appears only inside
    `V⁺`, and its corner error still stalled the conical solve. On the equal-arc
    mesh it is the other way round, and that mesh keeps the single operator.
- **Kummer acceleration.** The raw Green's sum converges like 1/m. The
  large-|m| asymptote sums in closed form to two logarithms (eqs. 4.82–4.85); the
  paired remainder decays like m⁻³ on the diagonal (eq. 4.87) and at worst like
  m⁻² on the boundary. `spectral_terms` defaults to `boundary_points / 2`, so
  that the single public knob converges both discretisations together (M&P's own
  Table 4.1 converges at P ≈ 2.2M).
- **Log singularity.** The integrable `ln` singularity of the kernel is handled
  by the periodic comparison-function split (eqs. 4.97–4.100): the singular part
  integrates analytically to `−L/π`, and the quadrature sees only a continuous
  remainder.
- **Neumann diagonal.** The double-layer kernel is continuous on a smooth curve.
  Its diagonal (which carries the local curvature) is taken as the numerical
  principal-value limit, the average of the kernel at the two physically adjacent
  nodes, which was measured to converge as O(spacing²).

## 5. Energy balance as a theorem

For a perfect conductor, `Σ_m ρ_m = 1` holds exactly over propagating orders
(M&P eqs. 4.34/4.42; thesis `eq:prop_order_unity`, by Green's theorem). With
finite conductivity the theorem becomes `R + A = 1` (G&S §2.C), and it remains a
two-sided check, because the absorption `A` is computed by an independent
boundary integral of the solved densities (§8) and is not taken as `1 − R`. For a
lossless substrate the same integral equals the transmitted power, and
`R + T = 1` is the theorem. In either case this is the sharpest self-check that a
solve can carry. The solver reports the deviation and never rescales: the
deviation of the computed conserved sum from 1 is the estimate of discretisation
error, recorded in `provenance.notes["energy_balance_deviation"]` and warned
about above 0.5%. `checks.check_energy_balance` uses the scan's recorded
absorption automatically.

The measured behaviour is as follows.

- **TE** holds the balance to ~1e-5 or better on every profile tried, corners
  included, at moderate `boundary_points`.
- **TM on smooth profiles** is comparable. On profiles with corners (ideal
  sawtooth, lamellar) the Meixner edge behaviour meets a plain equal-arc mesh and
  TM converges at first order, with a percent-level deficit at P ≈ 100 that is
  still ~0.5% at P ≈ 1000. The convergence harness reports this accurately.
- **TM on profiles with corners, graded mesh** (`corner_grading`), measured on
  2026-10-23 in-plane at period 1000 nm (lamellar depth 0.3 at α = 10°, 600 nm;
  ideal 30° sawtooth at α = 20°, 500 nm; 29.5°/70.5° blaze at α = 15°, 450 nm).
  The deviation of the sum from 1:

  | profile | mesh | P = 64 | 128 | 256 | 512 |
  |---|---|---|---|---|---|
  | lamellar | equal arc | 8.5e-3 | 4.3e-3 | 2.3e-3 | 1.3e-3 |
  | lamellar | graded, p = 4 | 5.9e-4 | 8.7e-5 | 1.3e-5 | 1.8e-6 |
  | sawtooth | equal arc | 7.5e-2 | 7.5e-2 | 2.2e-2 | 2.9e-2 |
  | sawtooth | graded, p = 4 | 3.2e-3 | 2.2e-4 | 4.9e-6 | 6.6e-6 |
  | 29.5°/70.5° | equal arc | 3.9e-2 | 2.8e-2 | 2.0e-2 | 1.4e-2 |
  | 29.5°/70.5° | graded, p = 4 | 6.8e-5 | 6.6e-5 | 2.0e-5 | 4.2e-6 |

  The efficiencies settle with it: the sawtooth's order −1 moves by 9e-6 from
  P = 256 to 512 on the graded mesh and by 3.7e-3 on the equal-arc one. Order 6
  is no better than 4 on these profiles, and orders above 6 put the two nodes
  nearest a corner within rounding of each other (NaN at order 8, P = 512), so
  the solver accepts 2 to 6 and refuses a mesh whose closest nodes it cannot
  separate. TE needed no help with the balance, but see below.
- **TE on profiles with corners passes the balance and is still wrong.** On the
  sawtooth the equal-arc TE balance reads 2.6e-10 at P = 512 while its order −1
  sits 6e-5 from the graded answer and moves by 2e-4 between rungs. The theorem
  holds for any density that solves the discrete system, so it cannot see an
  error that is common to every order. Read the convergence harness, not the
  balance alone, on a cornered TE solve; the graded mesh settles it to 5e-6.
- **Finite conductivity** holds `R + A = 1` to ~1e-4 on smooth profiles
  (measured: 1.2e-4 for Au at 2 nm, 1.5° graze, P = 256; 1e-5 on a flat
  interface). On the equal-arc mesh corners limit both polarizations at first
  order, since the edge-singular densities meet evenly spaced nodes, and this
  is reported through the same deviation.
- **Finite conductivity on a graded mesh**, measured on 2026-10-24 on the
  coupled system directly at fixed indices, with the two corrections of §4. The
  visible-light metal is n = 0.2 + 3.5i at 600 nm, period 1000 nm, α = 10°
  in-plane, on `Lamellar(0.2)` and the 29.5°/70.5° blaze. The X-ray cases are
  n = 0.998 + 0.002i at 2 nm, period 315.15 nm, off-plane at 1.5° graze and 25°
  azimuth, on `Lamellar(0.05)` and the same blaze. The deviation of `R + A`
  from 1 for the TM incident state:

  | case | mesh | P = 128 | 256 | 512 |
  |---|---|---|---|---|
  | visible, lamellar | equal arc | 1.3e-2 | 7.1e-3 | 9.3e-3 |
  | visible, lamellar | p = 4 | 1.3e-4 | 1.9e-5 | 3.4e-6 |
  | visible, lamellar | p = 6 | 1.6e-4 | 2.3e-5 | 4.0e-6 |
  | visible, blaze | equal arc | 1.0e-2 | 8.3e-3 | 4.9e-3 |
  | visible, blaze | p = 4 | 1.3e-3 | 1.7e-4 | 2.4e-5 |
  | visible, blaze | p = 6 | 1.7e-3 | 2.2e-4 | (3.1e-5) |
  | X-ray, lamellar | equal arc | 1.7e-1 | 1.3e-1 | 1.6e-1 |
  | X-ray, lamellar | p = 4 | (4.1e-4) | 1.6e-5 | 2.0e-6 |
  | X-ray, lamellar | p = 6 | (6.1e-4) | 2.7e-5 | 2.2e-6 |
  | X-ray, blaze | equal arc | 1.6e-1 | 1.5e-1 | 8.3e-2 |
  | X-ray, blaze | p = 4 | (2.3e-4) | 2.7e-5 | 2.9e-6 |
  | X-ray, blaze | p = 6 | (2.5e-4) | 3.2e-5 | (1.7e-6) |

  Bracketed rungs are ones the solver refuses: the graded X-ray meshes at
  P = 128 fall below the resolution floor (151 to 207 nodes needed), and p = 6
  at P = 512 on the 1000 nm blaze and on the X-ray blaze puts two nodes closer
  than the separation floor. That floor is earning its keep: past it the same
  solve degrades (the visible blaze reads 6.8e-2 at p = 6, P = 1024). The TE
  state reads the same as TM off-plane, where the two couple, and in-plane
  falls from 1e-4 and 3.9e-4 on equal arc to 2.7e-6 and 7.0e-6 at p = 4. Orders
  2, 3 and 5 converge with the rest (X-ray lamellar, p = 5: 2.0e-5, 2.1e-6,
  8.8e-7 at P = 256, 512, 768).

  Without the `K⁺` subtraction the visible-metal TM stalled between 3e-3 and
  5e-3, and the X-ray lamellar read 1.2 at P = 256, worse than the equal-arc
  mesh. Without the
  composition `D_t (V⁻ w)` the conical cases converged erratically, and at
  order 5 the X-ray lamellar stalled near 7e-3 from P = 256 to 768. The
  previous version of this section blamed oscillatory corner singularities from
  the negative permittivity. They are not the cause: the TM corner exponent at
  ε = n² ≈ −12.2 + 1.4i is ν = 0.600 − 0.009i at the lamellar's corners and
  0.564 − 0.010i at the blaze apex (against 2/3 and 0.643 for a perfect
  conductor), real to two decimals, and the grading resolves it like any other.

  **The balance cannot tell the two `D_t V⁻` paths apart, and the composition
  is the one that converges.** On the X-ray lamellar the single operator
  balances to 1.1e-6 at P = 1024 and the composition to 8e-7 at P = 768, yet at
  P = 1280 their efficiencies differ by 2.9e-5. Between rungs the single
  operator's efficiencies still move by 4.6e-5 (P = 768 to 1024) and 1.9e-5
  (1024 to 1280), while the composition's move by 3.2e-6 and 1.1e-6. This is the
  same blindness as the TE sawtooth above.

  The `K` subtraction is not applied to the perfectly conducting TM equation,
  which is second kind and was not limited by this error. Tried there, it helps
  the lamellar (1.8e-6 to 3.7e-7 at P = 512) and hurts the ideal sawtooth (6.6e-6
  to 1.6e-5, and 25 times worse at P = 256), so the perfectly conducting results
  stay exactly as measured above.

## 6. Validation evidence

- **Flat mirror**: unit specular efficiency, with the signs of the amplitudes
  pinning the boundary conditions (TE r₀ = −1, TM r₀ = +1).
- **M&P Table 4.1** (sinusoid, Littrow, their own P = 110, M = 50): TE
  0.4659/0.5341 and TM 0.9579/0.0421 reproduced to 1e-3
  (`tests/test_integral_core.py`).
- **Maréchal–Stroke theorem** (echelette, TM, exact blaze): approached
  monotonically, at first order because of the corner limitation.
- **Lorentz reciprocity** (`checks.check_reciprocity`): TE to ~1e-8, TM to
  ~1e-5 on smooth profiles. On cornered profiles at P = 256 the graded mesh
  takes TM from 2.3e-2 to 1.3e-4 (sawtooth) and from 6.1e-4 to 9.3e-6
  (lamellar), an independent check of the table in §5, since reciprocity is not
  the theorem the deviation measures.
- **An independent integral-method code, orderwise, on a measured groove**: the
  TASTE perfect-conductivity wavescan in the private reference corpus is
  reproduced across every propagating order of the full 0.6–6.0 nm scan, within
  the reference's own energy-balance error, and the comparison identified the
  previously unconfirmed profile behind that run (`docs/findings.md`;
  `tests/test_integral_corpus.py`). The figures remain with the private corpus.

Finite conductivity adds its own ladder of tests (`tests/test_integral_finite.py`):

- **Flat interface = Fresnel, amplitude and phase**: TE against `r_s` and TM
  against `r_p` from `materials/fresnel.py`, to 1e-4 in-plane. In the conical
  mount the full 2×2 (E_z, B_z) reflection matrix, polarization conversion
  included, matches the rotated-Fresnel closed form to 6e-5. This test succeeds
  the sign checks on the flat mirror, and every sign in the coupled system rests
  on it.
- **Energy theorems**: `R + T = 1` (glass sinusoid) and `R + A = 1` (gold
  sinusoid) to a few 1e-4. For a lossless substrate the boundary-integral `A`
  equals the transmitted Rayleigh sum, so two unrelated post-processings of one
  solve agree.
- **Perfectly conducting limit**: at n = 100i the solve lands on the validated
  PC solver to 1.3e-3 in the worst order, with absorption at numerical zero. The
  residual is the ~1 nm field penetration that the mode deliberately adds and is
  not an error. Along a lossy path (n = |n|·e^{iπ/4}) the absorbed fraction falls
  monotonically toward the conductor.
- **Goray & Schmidt Table 3** (dielectric sinusoid, conical, pure E_z incidence;
  itself a cross-method comparison against Li's CM): every tabulated order,
  reflected and transmitted, reproduced within 6e-5. Their Table 4 was found to
  be a publisher's duplication of Table 3 (`docs/findings.md`).
- **Lorentz reciprocity survives absorption**: measured at ~5e-5 for Au at X-ray
  grazing incidence through `checks.check_reciprocity`, unchanged.

## 7. Choosing `boundary_points`

The solver refuses fewer than ~6 nodes per reduced wavelength along the boundary
arc, since below that density the field is unresolved and the answer would be
wrong while appearing credible. The tabulated mode applies the same floor to the
metal-side transverse wavelength `λ / |√(n² − cos²γ)|`, which is the shorter of
the two whenever the mount sits below the critical angle. The field then decays
into the metal on the `λ/√(2δ)` scale; for Au at 2 nm and 1.5° graze that scale
is ~21 nm against a 76 nm vacuum reduced wavelength, giving a floor of ~210
points on a 600 nm sinusoid. The metal-side wavelength is much shorter still for
visible-light metals (|n| of a few), which is why that regime costs nodes here
(M&P §4.6.4). Above the floors, the following applies.

- Smooth profiles: convergence is fast, and the harness typically certifies P in
  the low hundreds.
- Cornered profiles + TM: on the default equal-arc mesh the energy-balance
  deviation shrinks only linearly. Pass `corner_grading=4`, which reaches 1e-5
  at P ≈ 256 on the perfectly conducting profiles in §5 and 2e-5 to 2e-4 on the
  finite-conductivity ones, against 7e-3 to 0.15 on equal arc. Either way,
  `check_convergence` (which sweeps this knob) should be left to pick the
  plateau, in preference to trusting a single run. On a graded mesh the
  resolution floor applies to the coarsest spacing, about twice the mean, so it
  asks for roughly twice the nodes.
- Cost: the kernels are assembled as separable BLAS products
  ([`_kernels.py`](../../src/gratinglab/solvers/integral/_kernels.py)), namely a
  few `(P, 2M+1)·(2M+1, P)` GEMMs per wavelength plus a one-time O(P²·M)
  geometry pass (the cached Kummer tails) amortised over the scan. Measured at
  P = 400, a solve takes ~100 ms per wavelength unpolarized, so the default 400
  at 200 wavelengths takes ~20 s. (The kernels were previously evaluated
  pointwise over `(P, P, M)` broadcasts at ~14 s per wavelength. Contrary to what
  this document used to claim, the O(P³) dense LU was never the bottleneck: it
  takes milliseconds at any realistic P, ~0.1% of the old runtime.) The GUI runs
  the solve on the worker thread with a live progress bar and a Cancel button
  that stops the solve itself. The tabulated mode costs ~2× the perfect one per
  wavelength: it needs more kernels and operator products, although the
  double/adjoint/tangential layers share one kernel evaluation per medium. Both
  incident polarizations share one factorization, so `"unpolarized"` costs two
  triangular solves and not two factorizations.
- Formulation envelope: the Kummer acceleration cancels asymptote terms of
  magnitude `exp(|α₀|·|Z|)`, so for `|α₀|·h ≳ 30` (`h` being the groove height,
  which is reachable only by driving a short transverse wavelength through a deep
  groove, e.g. an in-plane mount at hard-X-ray wavelengths on a large period)
  double precision has no significant digits left and the kernels are
  meaningless in any formulation of this acceleration, old or new; this was
  measured to affect both equally. The grazing-conical mounts this solver targets
  keep `|α₀| = k sin γ · sin α` small precisely because `sin γ` is small. If the
  regime is ever reached, the energy-balance theorem check fails loudly, with R+A
  deviations of order 1 reported in provenance, instead of failing silently.

## 8. Finite conductivity: the coupled conical system

Implemented in [`_finite.py`](../../src/gratinglab/solvers/integral/_finite.py),
following Goray & Schmidt (2010). The formulation was derived independently in
this project's conventions and was not transcribed, because G&S's normal points
into the metal and their potentials carry a factor 2, and because sign-sensitive
equations do not survive PDF extraction (`references.md`). The construction, in
this project's operators (outward normal up, `_nystrom` scalings), is as follows.

- The two coupled scalars are `E_z` and `B_z = (μ⁺/ε⁺)^{1/2} H_z`, each
  satisfying a transverse Helmholtz equation with the wavenumber of its own side:
  `k_t⁺ = k sin γ` (exactly the reduced-wavelength problem of §2) and
  `k_t⁻ = k √(n² − cos²γ)` (the same Green's-function code with complex k; the
  branch pin in `_greens` was built for this case).
- The substrate fields are single-layer potentials, `u⁻ = S⁻w` and
  `v⁻ = S⁻τ`, with traces `V⁻w` and normal derivatives `(L⁻ − I/2)w` given by
  the jump relation that the `_nystrom` tests pin.
- On the vacuum side, Green's representation of the radiating scattered field,
  together with the regularity identity of the incident field, gives the exact
  boundary relation `(I/2 + K⁺)ψ − V⁺ ∂_n ψ = ψ_i` for the total field. Its
  perfect-conductor limits are exactly the two milestone-1 equations
  (`∂_n ψ = 0` gives TM, `ψ = 0` gives TE), and this correspondence is what pins
  its signs.
- The jump conditions (G&S eq. 6, mapped) eliminate the vacuum traces: `E_z`
  and `B_z` are continuous, and

      ∂_n E_z⁺ = c_E ∂_n E_z⁻ + s ∂_t B_z,
      ∂_n B_z⁺ = c_B ∂_n B_z⁻ − s ∂_t E_z,

  with `c_E = n² k_t⁺²/k_t⁻²`, `c_B = k_t⁺²/k_t⁻²`, and
  `s = cos γ (1 − k_t⁺²/k_t⁻²)`. The tangential derivatives act on the
  continuous traces, so they are `D_t V⁻` on the densities. On the equal-arc
  mesh that product is one operator, the derivative taken analytically on the
  kernel (`tangential_layer_matrix`), because a spectral derivative of an
  assembled `V⁻` leaves an O(1) error there. On a corner-graded mesh it is that
  composition after all, the spectral τ-derivative of the assembled trace,
  which is smooth in τ there (§4, §5).

The result is a 2×2 block system in `(w, τ)`, with diagonal blocks
`(I/2 + K⁺)V⁻ − c V⁺(L⁻ − I/2)` and cross blocks `∓ s V⁺ D_t V⁻`. In-plane
(`cos γ = 0` exactly, since it is computed from `sin γ` so that the cross blocks
vanish identically) the system decouples into the classical TE/TM transmission
problems. In the perfect-conductor limit `|n| → ∞` it degenerates to the
milestone-1 equations (`V⁻` and `c_B` vanish like 1/|n|).

For each reflected order the efficiency is `(γ_m/γ₀)(|E_m|² + |B_m|²)` for a
unit-power incident state. Both components are included because conical
incidence converts polarization (TE ≡ incident `(p_z, q_z) = (1, 0)`, TM ≡
`(0, 1)` in the (E_z, B_z) basis; `conventions.md` §7). Absorption is the
boundary integral of G&S eq. 26 mapped to this project's operators. On a flat
interface it reduces exactly to `1 − |r|²` in both polarizations, and that
reduction is how its sign and normalisation were pinned analytically before any
test ran.

- D. Maystre and E. Popov, "Integral Method for Gratings", ch. 4 of
  E. Popov (ed.), *Gratings: Theory and Numeric Applications*, AMU (2012).
- J. A. McCoy, PhD thesis, Chapter 2 §2.2 (read with `conventions.md` §10;
  §2.2.3 is superseded, see item 5 there).
- L. I. Goray and G. Schmidt, JOSA A 27, 585 (2010): the finite-conductivity
  formulation (§8), implemented at M19.
- L. Li, JOSA A 10, 2581 (1993): the coordinate-transformation results behind
  G&S's comparison tables, the external anchor of §6.
- H. A. Kalhor and A. R. Neureuther, JOSA 61, 43 (1971): an early
  perfectly-conducting integral treatment.
