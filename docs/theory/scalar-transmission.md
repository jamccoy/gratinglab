# Scalar transmission: a design draft (#83)

**Status: proposal, 2026-11-20. Nothing here is implemented.** It answers the
three open questions in issue #83 and records the measurements they rest on.
The app half is magmacrunch-media/gratinglab#44, which draws transmitted
orders today and has no efficiency to give them.

## 1. Two axes, not one

A transmission grating is not one problem for a scalar model. Two questions
decide what scalar theory can say about it:

- **Continuous relief or discrete bars.** A relief cut into a substrate
  (a blazed or holographic grating in glass) has one interface, and the light
  ends up in the substrate. Bars standing free in vacuum, perhaps on a thin
  support membrane (an X-ray phase grating, HETG-style bars), have two
  interfaces per bar, and the light ends up back in vacuum. The phase profile
  is similar. The flux factor, the material each one needs, and what
  `Problem` has to express are all different.
- **Thin or thick: can propagation inside the groove region be neglected?**
  A thin element changes the phase at a single plane and nothing else. A thick
  one makes the ray cross the structure on a slant, refract, guide or reflect
  off sidewalls inside it, and the orders couple to each other on the way
  through (the Bragg regime). Scalar theory is a thin-element theory, so thick
  is not a model choice. It is the boundary of the method.

| | thin | thick |
|---|---|---|
| **continuous relief** | slice 1: scalar, checkable against RCWA today | RCWA / C-method only |
| **discrete bars** | slice 2: scalar, needs `Problem` work first (§6) | RCWA only (e.g. CAT gratings) |

### Thin versus thick, measured: the Klein–Cook parameter

The usual criterion is Klein and Cook's

```
Q = 2 π λ h / (n Λ²)
```

where h is the groove depth, Λ the period and n the index of the medium the
light crosses the grooves in. Q ≪ 1 is thin (Raman–Nath), Q ≫ 1 is thick
(Bragg). It reads correctly on both of our home regimes:

| grating | λ | h | Λ | Q | regime |
|---|---|---|---|---|---|
| HETG-like bars | 1 nm | 500 nm | 200 nm | 0.08 | thin, scalar applies |
| CAT grating | 2 nm | 4 µm | 200 nm | 1.3 | thick, which is how it blazes: sidewall reflection |
| relief in glass, below | 500 nm | 1 µm | 2.5 µm | 0.34 | already failing |
| relief in glass, below | 500 nm | 1 µm | 10 µm | 0.021 | thin |

In a conical mount the reduced wavelength λ / sin γ goes in place of λ, the
same reduction `_REDUCED_RATIO_WARN` already makes.

Note that **thin is not the same as "neglect refraction"**. A single relief
interface still refracts. The thin approximation neglects what happens inside
the groove region, not the bend at the boundary. §2 is about that bend.

## 2. Which phase: does the order's direction enter?

Two thin-element phases are available for a relief between vacuum and a
substrate of index n, for order m leaving at β_m inside the substrate
(n sin β_m = sin α + mλ/Λ):

- **order-free**: Φ(t) = k g(t) (n cos β₀ − cos α), with β₀ the Snell angle of
  the zero order. One transmittance function, so a single FFT, and the sum of
  efficiencies conserves energy (Parseval).
- **order-dependent**: Φ_m(t) = k g(t) (n cos β_m − cos α). This is the
  transmissive twin of `ScalarSolver`'s reflective k g (cos α + cos β_m). It is
  symmetric under swapping the incident and diffracted waves, so it is
  reciprocal, and it puts the blaze where a refracting prism facet sends the
  light at every α, not only at one angle.

This is the same three-way conflict as in [`scalar.md`](scalar.md) §5, carried
over unchanged: energy conservation, reciprocity and blaze direction cannot
all hold at once. **The proposal is the order-dependent phase**, so that the
scalar solver keeps one invariant, reciprocity, in both directions. The
energy deviation would then be reported, never renormalised away, exactly as
for reflection. The alternative is to offer both. That is a judgement call
and is left open in §7.

At normal incidence the two phases differ only through cos β_m against
cos β₀, and the measurements below cannot tell them apart. At α = ±30° they
can: they differ by up to 0.13 on the blaze peak. That comparison is not
conclusive yet (§5).

## 3. The flux factor

Proposed: the flat-interface TE Fresnel transmittance,

```
O_m = 4 cos α · n cos β_m / (cos α + n cos β_m)²
```

It has the same form as the reflective solver's symmetric obliquity and is
reciprocal for the same reason. In the shallow-groove limit it gives the
zero order the plain Fresnel transmittance. Measured: 0.8683 against RCWA's
0.8678 for a 100 nm lamellar groove at Λ = 10 µm.

Scalar theory ignores polarization, and TE and TM Fresnel differ at
oblique incidence. The candidates are: TE for every polarization, which is
what was measured; the polarization asked for; or the mean for unpolarized.
Using the polarization asked for is the least surprising choice. The
solver's degeneracy note on the provenance would then say that only the
interface factor depends on polarization.

For discrete bars in vacuum at X-ray wavelengths, n ≈ 1 − δ + iβ, so the
interface factor is 1 to within δ. The bars' absorption is what matters, and
it enters through the imaginary part of the phase. A bar is then an amplitude
grating as well as a phase grating, the transmissive counterpart of the
`"local"` reflectivity model.

## 4. Material index: does a table come back?

- **Continuous relief** needs the substrate's real index n. `Transmission`
  already defines a transmitted order only for a lossless substrate, so the
  relief case is real-n by construction. The package takes it through the
  materials layer as it does now. The app can take one typed number, n at the
  wavelength, and ship no table.
- **Discrete bars** need a complex n, so in practice δ and β, which for X-rays
  means CXRO/Henke-type data. The package already has the answer: a user
  table, `typical:` or `factors:`. The app can again take typed values and
  embed nothing. So **no table has to come back into either half**. The
  public-domain-data rule for the app holds without a decision.

## 5. Measurements

Run on 2026-11-20 against `rcwa` (truncation 40 at Λ = 2.5 µm, 90 at 10 µm):
λ = 500 nm, TE, a relief cut into a lossless n = 1.5 table, transmitted into
the glass. "Thin" is the order-free phase and "refr" the order-dependent one,
both with the §3 factor. `python benchmarks/scalar_transmission.py`
reproduces both tables in under a minute.

Normal incidence, efficiency of the orders named:

| grating | Q | order | RCWA | thin | refr |
|---|---|---|---|---|---|
| lamellar, h 100, Λ 2.5 µm | 0.034 | 0 / ±1 | 0.8645 / 0.0377 | 0.8683 / 0.0372 | 0.8683 / 0.0353 |
| lamellar, h 500 (π step), Λ 2.5 µm | 0.17 | ±1 / ±2 | 0.3824 / 0.0142 | 0.3898 / 0 | 0.3891 / 0 |
| lamellar, h 1000 (2π step), Λ 2.5 µm | 0.34 | 0 | **0.6897** | 0.9600 | 0.9600 |
| blazed into 1st, h 1000, Λ 2.5 µm | 0.34 | blaze | **0.7827** | 0.9617 | 0.9594 |
| lamellar, h 100, Λ 10 µm | 0.002 | 0 / ±1 | 0.8678 / 0.0372 | 0.8683 / 0.0372 | 0.8683 / 0.0370 |
| lamellar, h 500, Λ 10 µm | 0.010 | ±1 | 0.3889 | 0.3891 | 0.3891 |
| lamellar, h 1000, Λ 10 µm | 0.021 | 0 | 0.8799 | 0.9600 | 0.9600 |
| blazed into 1st, h 1000, Λ 10 µm | 0.021 | blaze | 0.9181 | 0.9601 | 0.9601 |

Blaze peak at α = ±30°, same blazed gratings:

| Λ | α | RCWA | thin | refr |
|---|---|---|---|---|
| 2.5 µm | −30° | 0.6014 | 0.9270 | 0.9365 |
| 2.5 µm | +30° | 0.6161 | 0.9050 | 0.8085 |
| 10 µm | −30° | 0.8231 | 0.9167 | 0.9344 |
| 10 µm | +30° | 0.8669 | 0.9112 | 0.8886 |

What these say:

- **Shallow and thin, both phases agree with RCWA to about 0.5%.** That is
  where slice 1 can promise agreement.
- **Q is the right axis, and the depth of phase matters too.** At Q = 0.34
  both phases miss the blaze peak by 0.18. At Q = 0.021 they are still 0.04 to
  0.08 over on the 2π-deep gratings, where scalar theory predicts a groove
  that is invisible or perfectly blazed and the real edges scatter. A
  warning threshold therefore wants Q and the phase excursion together, the
  way the reflective solver's energy deviation scales with phase excursion.
- **The oblique rows are not yet evidence for either phase.** RCWA's ±30°
  asymmetry is small (0.60 against 0.62). The two scalar phases split by
  more and in the opposite direction. But the benchmark's sign conventions
  for α and for the phase were set by which order peaks, to agree with
  RCWA's, and were not checked against `docs/conventions.md`. Settle that
  before reading anything into those rows.
- The benchmark puts a negative sign on the phase, which puts the scalar blaze
  peak on the RCWA's order (−1). With a positive sign it lands on +1. The
  sign is a convention to take from the reflective solver when implementing,
  not a physics result.

## 6. What the API needs

- **Slice 1 (continuous relief).** `ScalarSolver.solve` fills
  `scan.transmission` whenever the problem's substrate is lossless at the
  wavelength, exactly as RCWA does. The meaning of `Problem` is unchanged:
  ridges of `coating` (or `substrate`) on that same material. The reflected
  side for such a problem then wants the substrate's Fresnel reflectivity in
  place of the coating's, which the existing reflectivity models already
  supply.
- **Slice 2 (discrete bars).** `Problem` cannot say "free-standing" today.
  RCWA's own docstring defines a problem as ridges on a semi-infinite
  substrate of the same material. A vacuum substrate, or a thin support
  membrane, needs a `Problem` field. RCWA's `Stack` already has homogeneous
  layers above and below the grating region, so the rigorous side looks like
  wiring, not new numerics. That has not been verified.
- **Thick.** No scalar code. Warn above a Q threshold, the way
  `_REDUCED_RATIO_WARN` warns, and say in the warning which solver to use.

## 7. Open, for a person to decide

1. The order-dependent phase alone, or both phases behind an option.
2. Interface factor: TE throughout, the polarization asked for, or the mean.
3. The warning threshold in Q, and whether the phase excursion joins it.
   In §5, every case with Q at or below 0.17 and a phase excursion of π or
   less agrees with RCWA to within 0.01 on its strong orders. Every 2π case
   misses by 0.04 or more, even at Q = 0.021.
4. Whether slice 2 is wanted at all, or whether X-ray transmission gratings are
   always RCWA territory in practice. HETG-like bars at Q ≈ 0.08 suggest
   scalar would be useful there.

## 8. How it would be checked

- Shallow limit: zero order to Fresnel T, and R₀ + T₀ → 1 with the substrate's
  Fresnel R on the reflected side.
- Closed forms as tests, as in `scalar.md` §4: a lamellar thin grating gives
  T±1 = (4/π²) sin²(φ/2) times the factor, with φ = k(n − 1)h.
- Reciprocity to machine precision for the order-dependent phase, by the same
  swap test the reflective solver uses.
- Against RCWA over a grid in Q and depth. The table in §5 becomes a test,
  and its agreement gate is the threshold chosen in §7.
- `tools/fixture.json` in the app gains transmitted rows from the same run,
  so the port is held to the package as it is for reflection.
