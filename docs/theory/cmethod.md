# C-method (Chandezon coordinate transformation)

Implemented in [`solvers/cmethod/`](../../src/gratinglab/solvers/cmethod/__init__.py).
Symbols and sign conventions follow [`conventions.md`](../conventions.md). The
method is that of Chandezon, Maystre & Raoult (1980) and Chandezon, Dupuis,
Cornet & Maystre (1982), in the formulation Li, Chandezon, Granet & Plumey
(1999) recast in terms of Li's Fourier factorization rules. The derivation in
§3 was written out for this implementation and is checked numerically (§6)
rather than transcribed; where it agrees with the literature that is stated,
and where it was only checked by agreement with other methods that is stated
too.

---

## 1. What the model assumes

The groove is a single-valued height function `y = a(x)` with a **finite
slope everywhere**. The method is rigorous for any such profile: Maxwell's
equations are solved exactly in a coordinate system in which the boundary is
flat, and the only approximation is the truncation of the Fourier series.

| Mode | Physics | Normalisation |
|---|---|---|
| `"perfect"` (default) | Perfect conductor, any mount. TE is `E_z = 0` on the boundary, TM is `∂H_z/∂n = 0`. | Relative to a perfect reflector; the sum over orders is 1 in the continuum. |
| `"tabulated"` | One interface into a semi-infinite material read per wavelength from `coating` (else `substrate`). Absorbing, or lossless with transmitted orders. | **Absolute**. The flux into the material is computed from the boundary fields. Into a lossless material the transmitted orders are recorded as the scan's `transmission`, with their own order axis, and `absorption` keeps the remainder, so `R + T + A = 1` is checkable. |

The two options are named exactly as the integral solver names them, so one
`Problem` and one option dictionary reach either method and their results are
comparable without an adapter.

What it refuses, with `UnsupportedConfiguration`, rather than approximating:

| Refused | Why |
|---|---|
| an undercut profile | the coordinate `u = y − a(x)` needs `a` to be a function |
| a vertical facet or sidewall: `Lamellar`, `Blazed(antiblaze_angle=90)`, a `FromProfileData` with two heights at one `t`, or ends that differ in height | `a'` is infinite there; any finite stand-in would be a smoothed facet, the failure `solvers/base.py` names |
| a wall too steep for the scan: one whose floor (§6, "Steep walls") exceeds the largest truncation accepted, 256 | the answer is wrong at every truncation the solver takes. Issue #7: an 88.8° wall at 1.55 nm off-plane returned efficiencies of 11 at N = 128 |
| a truncation below a steep wall's floor (`ValueError`, naming the floor) | measured wrong by 3.5e-3 to 6e2 there |
| any `conductivity` other than the two above | as for the integral solver |

A named `coating` is not a layer: as in the integral solver, it is the
semi-infinite material under the boundary. Coatings and layer stacks are
deferred (§5).

## 2. The coordinate transformation

Grooves run along `z` and every field carries `exp(i k_z z)` with
`k_z = k cos γ` conserved (`conventions.md` §3). In a homogeneous medium of
permittivity `ε`, both `E_z` and `H_z` (H scaled by the vacuum impedance)
satisfy the transverse Helmholtz equation

    ∇_t² F + κ² F = 0,    κ² = k² ε − k_z²,

and every other component follows from these two. The vacuum side has
`κ = k sin γ`, which is why the conical problem looks like an in-plane one at
the reduced wavelength `λ / sin γ`.

The C-method's change of variables is the translation

    x' = x,    u = y − a(x),

which maps the boundary onto `u = 0`. The Jacobian is 1, and with the metric
`g^xx = 1`, `g^xu = −a'`, `g^uu = 1 + a'²` the Helmholtz equation becomes

    ∂_x(∂_x F − a' ∂_u F) + ∂_u(−a' ∂_x F + (1 + a'²) ∂_u F) + κ² F = 0.

The second bracket is the quantity matched at the boundary,

    Φ = −a' ∂_x F|_y + ∂_y F = √(1 + a'²) ∂F/∂n,

the normal derivative weighted by the line element. It is what the TM
condition sets to zero, what flux continuity is written in, and what enters
the Poynting flux through the boundary.

The profile parameter runs against the dispersion direction, `t̂ = −d̂`, so
the solver places the curve at `x = (1 − t) P`, as the integral method does
(`solvers/cmethod/_series.py` is the one place this reflection happens).

## 3. The eigenproblem, and the factorization

Floquet periodicity gives `F = exp(iρu) Σ_n F_n exp(iα_n x)` with
`α_n = −k sin α sin γ + 2πn/P`. The derivation keeps three unknowns that
are **continuous along a line `u = const` even across a facet corner**: `F`,
`P = ∂_y F` and `Q = ∂_x F|_y`. A field in a homogeneous medium is smooth in
`(x, y)`; only the direction of the line `u = const` jumps at a corner. In
these unknowns

    ∂_x' F = Q + a' P,    ∂_u F = P,    ∂_x' Q − a' ∂_u Q + ∂_u P + κ² F = 0,

and every product is of the discontinuous `a'` with a continuous factor. By
Li's rules (Li 1996) Laurent's rule, the plain Toeplitz product, is then
correct in every term, and no inverse rule is needed anywhere. With `A` the
Toeplitz matrix of `a'` and `α = diag(α_n)`, eliminating `P` and `Q` gives

    ρ² (I + A A) F − ρ (α A + A α) F − (κ² − α²) F = 0.

`A A` is the **matrix product**, not the Toeplitz matrix of `a'²`. The two
agree for a smooth profile in the limit and differ at every finite truncation
on a profile with corners; this formulation produces the former, and §6 shows
the result converging to the integral method's answer on sawtooth profiles.
That is a numerical argument for the factorization, not a citation of one.

It is linearised as `ρ [F; G] = [[0, I], [(I+AA)⁻¹(κ²−α²), (I+AA)⁻¹(αA+Aα)]] [F; G]`,
a dense `2(2N+1)` eigenproblem per medium per wavelength. The flux of a mode
in the same factorization is `Φ = −i A α F + i ρ (I + A A) F`.

**The propagating orders use exact plane waves.** A plane wave
`exp(i(α_m x + β_m y))` is an exact solution, and in the new coordinates it is
`exp(iβ_m u) Σ_n L_{n−m}(β_m) exp(iα_n x)` with `L_j(β)` the Fourier
coefficients of `exp(iβ a(x))`. From `a' e^{iβa} = (iβ)⁻¹ d/dx e^{iβa}` its
flux is exactly `Φ_n = i L_{n−m} (κ² − α_m α_n)/β_m`. So the real eigenvalues
of the truncated problem, which approximate `±β_m`, are discarded and replaced
by these. Each reflected amplitude is then the coefficient of a true outgoing
plane wave, and the efficiency is read directly:
`e_m = (β_m/β_0)(|E_m|² + |H_m|²)` per unit incident `|E|² + |H|²`
(`conventions.md` §6, §7).

The split relies on counting. The vacuum problem must have exactly two real
eigenvalues per propagating order; if it has more, which a truncated problem
can produce on a deep or steep profile, the solve is refused with a
`ValueError` naming the count, rather than guessing which is spurious. The
evanescent modes kept are those with `Im ρ > 0` above the boundary and
`Im ρ < 0` below it.

The Fourier coefficients are computed in closed form, segment by segment, for
a polygon (`Blazed`, `FromProfileData`), so a corner carries no sampling
error; and by FFT for a smooth profile (`Sinusoidal`), where the trapezoidal
rule converges exponentially.

## 4. The boundary conditions

At `u = 0` the tangential fields are continuous. With the tangent
`t̂ ∝ (1, a')` and the transverse fields expressed through `E_z` and `H_z`,

    √(1+a'²) E_t = (i/κ²)(k Φ^H + k_z ∂_x' E_z)
    √(1+a'²) H_t = (i/κ²)(k_z ∂_x' H_z − k ε Φ^E)

where `∂_x'` is the derivative along the boundary, `iα_n` per harmonic. So
four series are matched, harmonic by harmonic: `E_z`, `H_z` and the two
brackets divided by `κ²`. In a conical mount `k_z ≠ 0` couples `E_z` and `H_z`
and polarization converts; in-plane the couplings vanish identically and TE
and TM separate into the familiar `E` and `Φ^E` continuity for TE and `H` and
`Φ^H/ε` for TM.

For a perfect conductor `E_z = E_t = 0` gives `E_z = 0` and `Φ^H = 0`, decoupled
in any mount. That is the conical reduction theorem the integral method also
uses, falling out of the same equations.

**Absorption** is the Poynting flux through the boundary into the material,
`−Re Σ_n (E_z,n (√ H_t)_n* − (√ E_t)_n H_z,n*)`, normalised by the incident
flux `k β_0 (|E|²+|H|²)/κ_1²`. It is computed from the boundary values, not as
`1 − R`, so `R + A − 1` is a genuine check. On a smooth profile it is a good
error estimate; on a cornered one it is not (§6).

## 5. What is deferred

- **Coatings and layer stacks.** The C-method extends to multilayer gratings
  with parallel or differing interfaces (Li, Chandezon, Granet & Plumey 1999),
  propagated by an S-matrix; none of that is here. A `Problem` names one
  material and this solver treats it as a semi-infinite substrate.
- **Adaptive spatial resolution** and other improvements for steep or deep
  profiles. On a sawtooth convergence is roughly first order (§6), and the
  integral method with `corner_grading` is the better tool there.
- **A GUI tab.** The solver is registered and appears in the Help menu; no
  tab is written for it yet (`_TAB_FACTORIES` skips a solver without one).

## 6. Validation evidence

Measured on 2026-10-27; `tests/test_cmethod.py` pins each.

| Check | Case | Result |
|---|---|---|
| Fresnel, amplitude and phase | flat boundary, gold (0.2+3i) and X-ray index, TE and TM, two angles | exact to rounding |
| Goray & Schmidt (2010) Table 3, itself Li's C-method | dielectric sine, n = 2, conical (θ = 60°, φ = 15°), all 4 reflected and 8 transmitted orders | within 1.4e-5 of every tabulated order; R + T = 1 to 1e-13 |
| Maystre & Popov Table 4.1 | perfectly conducting sine at Littrow | TM 0.95789 against 0.9579. TE 0.466441 against a printed 0.4659: the integral solver gives 0.466441 too, at every mesh from P = 110 to 512, so the gap is the table's |
| Maréchal–Stroke | 30°/60° echelette, TM at Littrow | unit efficiency at any truncation. The true field is two plane waves, so this pins the plane-wave machinery and the blaze side, not convergence |
| Integral solver, perfect conductor | sinusoid, TE and TM | to 1e-6 (TE) and the integral's own 1e-5 (TM) |
| Integral solver, finite conductivity | gold (0.2+3i) sinusoid in-plane and conical at γ = 60°; CXRO gold off-plane at γ = 1.5°, 2 and 5 nm | 2e-6 to 1.7e-5 per order for visible gold. Off-plane 8.8e-5 at 512 integral points and 2.6e-5 at 768: the C-method is converged there (N = 24 and 48 agree to 1e-14), so the gap is the integral's discretization, closing as its mesh refines |
| Energy theorem | every smooth case above | `Σ e` or `R + A` = 1 to 1e-10 |

**On a cornered profile the energy balance is not an error estimate.** On
Blazed 20°/70°, perfect conductor, TE, at N = 160 the balance reads 1e-8
while the efficiencies still differ from the integral method's graded-mesh
answer by 1e-4 per order, approaching it at roughly first order in `N`. The
solver says so in a provenance warning whenever the profile has a corner. The
convergence harness (three agreeing rungs) and a cross-method comparison
are the evidence to use there, not the balance.

**On a deep profile a high truncation makes the answer worse.** The vacuum
side's evanescent modes, written on the boundary, span about
`exp(|Im ρ| · depth)`, and once that passes double precision the solve loses
digits. Measured on sinusoids, perfect conductor, λ = 0.9 P, against the
integral solver at 768 points:

| depth / P | last good N | first bad N, error |
|---|---|---|
| 0.3 | 128 and beyond | none found |
| 0.5 | 32 | 48, 1e-5 to 2e-4 |
| 0.8 | 24 | 32, 3e-5 to 1e-3 |
| 1.0 | 24 | 32, 2e-3 to 3e-2; at 48 TM is wrong by 2.4 |

The estimated condition number of the boundary system tracks this (below
1e12 every case was right; from 1e14 the deep ones were wrong) but over-warns
on shallow profiles, where 0.3 P stays right to six digits at 1e17. So the
solver warns above 1e12 and does not refuse, and on a smooth profile the energy
balance, which tracks this error closely, is the better alarm. A naive remedy,
dropping the fastest-decaying modes and solving the boundary conditions by
least squares, was tried and made every case worse (errors of 0.1 to 0.9), so
it is not here.

**A steep wall needs harmonics in proportion to its slope and to
`P sin γ / λ`, and below a floor the answer is wrong.** Issue #7: on a groove
with an 88.8° wall, at 1.55 nm and γ = 1.73°, the C-method returned
efficiencies of 1.9e4 at N = 16 and 11 at N = 128 where the integral solver
gives 0.26. Before the fix only an exactly vertical wall was refused. The
mechanism is resolution, not the deep-profile precision loss above: a plane
wave written on the boundary, `exp(iκ a(x))`, oscillates
`(P κ / 2π) |a'|` times per period along a wall of slope `a'`, the modes on
either side must be resolved as finely, and every answer did converge to the
integral solver's once the truncation passed that. In the off-plane X-ray
mount `κ = k sin γ`, so 1.55 nm at γ = 1.73° on a 400 nm period acts like
an in-plane wavelength of 0.13 P, and a wall that is easy in visible light
is not easy there.

`wall_floor` reads two numbers off the profile's own Fourier series, the
same one the solve uses, so it applies to every profile type and to the
fillets of a rounded `Trapezoidal` (which are tangent to the wall and so
carry its slope): `B`, the harmonics that hold all but 1e-3 of the power of
`exp(iκ a)` at the largest transverse wavenumber in the scan, and the
steepest slope `s`. The floor is `0.75 max(B, 3 s)`. The `3 s` term is for
long reduced wavelengths, where `B` shrinks while a very steep wall still
needs its own harmonics. A lossless material's transverse wavenumber counts
where it is larger than the vacuum's. An absorbing one's does not: gold's
X-ray answers needed no more harmonics than a perfect conductor's, and
counting it took the floor of the default blazed groove in gold at 4 to
4.5 nm from 11 to 21, where N = 16 was already within 2e-3 of the integral
solver and N = 8 within 1.2e-2.

Measured against the integral solver (`corner_grading = 4`, 1200 points,
which moved by under 3e-5 from 800) on 24 configurations and 140
truncations: a symmetric trapezoid (depth 0.3 P, crest 0.3 P) and the
issue's asymmetric shape (a 27° facet, a 77 nm land, the steep wall) on a
400 nm period, walls from 45° to 89.5°, perfect conductor in-plane at 100,
500 and 700 nm and CXRO gold at γ = 1.73° and 1.55 to 6 nm. Largest
difference over the orders:

| truncation | rungs | difference from the integral solver |
|---|---|---|
| below the floor | 62 | 3.5e-3 to 5.8e2; under 1e-2 only four times, each one harmonic below the floor |
| floor to twice the floor | 25 | 1.2e-4 to 5.4e-2 |
| beyond twice the floor | 53 | 4.0e-3 or less |

The asymmetric groove in the issue's mount, gold, 1.55 to 6 nm:

| wall | floor | N = 32 | 64 | 128 | 192 | 256 |
|---|---|---|---|---|---|---|
| 80° | 44 | 4.8e-2 | 1.7e-3 | 1.3e-4 | 3.3e-5 | 1.6e-5 |
| 85° | 84 | 6.7 | 5.3e-2 | 2.1e-3 | 4.0e-4 | 1.6e-4 |
| 87° | 135 | 1.1e2 | 2.2 | 2.5e-2 | 4.0e-3 | 1.0e-3 |
| 88° | 195 | 5.8e2 | 28 | 0.35 | 2.3e-2 | 1.1e-2 |
| 88.8° | 311, refused | | | | | |

Swept without the floor, the C-method needed the same truncation on this
groove as a perfect conductor as it did in gold (the step at which its
energy balance fell below 1e-2 was the same at every wall from 65° to
88.8°), which is the conical reduction to `λ / sin γ` at work. In visible light the same walls
are cheap: the symmetric trapezoid at 500 and 700 nm has floors of 13, 26,
65 and 108 at 80°, 85°, 88° and 88.8°, and 89.5° (floor 258) is refused at
any wavelength, its error still 1e-2 at N = 192.

So the solver refuses a floor beyond 256 as `UnsupportedConfiguration`
(no truncation would answer it), refuses a truncation below the floor with
`ValueError` naming the floor (as it does for too few orders), and warns
between the floor and twice it. **There the energy balance is not an error
estimate**: at 87° and N = 128 it read 2e-3 against a true 2.5e-2. The GUI's
default scan sits just inside: its 70.5° facet needs 29 at 1 nm, and N = 24,
29 and 32 were 2.5e-2, 7.6e-3 and 3.0e-3 from the integral solver.

The floor depends on the wavelength, so it cannot sit in
`Capabilities.check`, which sees no scan. It runs at the start of `solve`,
beside the vertical-wall refusals, before any eigenproblem; the GUI's form
puts a floor the field can meet on the truncation field, and a floor beyond
256 reaches the C-method and Compare tabs as the solver's refusal.

Not yet validated: a measured (AFM) profile with many short segments.

## 7. Choosing `truncation`

`truncation = N` keeps harmonics `c − N .. c + N`, centred on the middle `c`
of the propagating orders at the shortest wavelength, so the off-plane X-ray
case, whose orders are far from symmetric about 0, does not waste half of
them. It must leave at least one evanescent harmonic beyond each end of the
propagating orders, counting the transmitted ones into a lossless material,
which for n > 1 reach further than the reflected ones; a truncation that
would drop a transmitted order is refused rather than leaving it out of the
record. Cost is dominated by two dense eigenproblems of size
`2(2N+1)` per wavelength, so it grows as `N³`; the accepted maximum is 256.

On smooth profiles convergence is exponential: a sinusoid of depth 0.3 P
converges to 1e-8 by N ≈ 10 to 20. On a sawtooth it is slow, and the default
32 is a starting point for a sweep, not an answer. A steep wall sets a floor
of its own, `wall_floor` (§6): a truncation below it is refused, and from it
to twice it the solver warns that the energy balance can understate the
error tenfold.

Two things a sweep can walk into, both deliberate. A rung above 256 raises
`ValueError`, and the convergence harness re-raises a refusal that follows a
successful rung, so a sawtooth swept on the default ladder (4, 8, 16, ... from
the default 32) stops with that error rather than reporting a plateau it never
reached. And on a deep profile the rungs past the instability (§6) disagree
with one another, so the harness finds no plateau there, which is the correct
report: the answer is the plateau below.

## References

- J. Chandezon, D. Maystre & G. Raoult, "A new theoretical method for
  diffraction gratings and its numerical application", J. Opt. 11, 235 (1980).
- J. Chandezon, M. T. Dupuis, G. Cornet & D. Maystre, "Multicoated gratings: a
  differential formalism applicable in the entire optical region", J. Opt.
  Soc. Am. 72, 839 (1982).
- L. Li, "Use of Fourier series in the analysis of discontinuous periodic
  structures", J. Opt. Soc. Am. A 13, 1870 (1996).
- L. Li, J. Chandezon, G. Granet & J.-P. Plumey, "Rigorous and efficient
  grating-analysis method made easy for optical engineers", Appl. Opt. 38,
  304 (1999).
- L. I. Goray & G. Schmidt, "Solving conical diffraction grating problems with
  integral equations", J. Opt. Soc. Am. A 27, 585 (2010), Table 3.
