# The typical-material model, and how far to trust it

`materials.typical("Au")`, or `coating="typical:Au"` on a `Problem`, is an X-ray
index built from physical facts rather than read from a table, because the
repository commits no tabulated optical constants
(`src/gratinglab/materials/data/SOURCES.md`). The derivation is in the docstring
of `src/gratinglab/materials/typical.py`; this page records what it is good for.

It covers 30 eV to 30 keV (0.0413 to 41.3 nm) and refuses anything outside that,
the way a table refuses wavelengths beyond its range. It knows Au, Pt, Ni, C, Si
and O, so SiO2 and other oxides of those, and takes a density for anything without
a default.

## The model in one paragraph

Each atom scatters with *f*₁ + i*f*₂ electrons, and a material is the sum of its
atoms. Every core shell, and the valence electrons as one more shell at the first
ionization energy, absorbs by Kramers' law (*f*₂ ∝ *E*⁻² above its edge, nothing
below), with a strength fixed by that shell's sum rule. The Kramers–Kronig
relation then gives *f*₁ in closed form, *Z** + Σ *N*ₛ *x*⁻² ln|*x*² − 1| with
*x* = *E*/*E*ₛ, which tends to *Z** far above every edge and dips at each one.
That dip is what the free-electron formula, *f*₁ = *Z*, misses: at 1 nm gold's
K, L and M shells are still bound, and the free-electron decrement is about 1.5
times too large.

## What is checked in CI

`tests/test_typical.py`, against the model's own physics, since there is no
public-domain table to compare with:

- *f*₁ recomputed by numerical principal-value integration of *f*₂ matches the
  closed form to 1e-6, for C, Si and Au at five energies;
- (2/π)∫*f*₂ d*E*/*E* = *Z* for every element, and *f*₁ → *Z** above all edges
  and *Z** − *Z* far below;
- linearity in density, *λ*² scaling far above the edges, additivity over the
  atoms of SiO2, and absorption above zero everywhere in range;
- the range, the names, the refusals, and that a caveat reaches a solver's
  provenance.

Breaking the Kramers–Kronig sign, the sum-rule factor or the relativistic term
fails between 11 and 17 of those tests (checked 2026-11-11, before the two
delayed-onset tests were added).

## Measured against CXRO, 2026-11-11

CXRO tables for each default material were downloaded to a scratch directory,
500 points from 30 eV to 30 keV at the model's default densities, compared, and
**not committed**. Figures are |model/CXRO − 1|. "Clear" is everywhere a scan
would get no caveat: more than 10% (or 30 eV) from every edge, outside the
delayed-onset band of a 4d or 4f level, and above 100 eV.

| | clear points | decrement, median | decrement, 90th pct | absorption, median | absorption, 90th pct |
|---|---|---|---|---|---|
| C | 397 | 0.2% | 9% | 25% | 83% |
| SiO2 | 339 | 0.5% | 7% | 9% | 48% |
| Si | 354 | 1.5% | 17% | 21% | 69% |
| Ni | 346 | 2.4% | 24% | 36% | 85% |
| Pt | 151 | 4.3% | 9% | 28% | 68% |
| Au | 145 | 3.7% | 12% | 27% | 66% |

So where it gives no caveat the decrement, which sets the critical angle, is
good to a few percent, and the absorption is a factor-of-two estimate. Gold at
8 keV reads about 6% low in decrement and 20% low in absorption.

The caveats cover most of gold's and platinum's soft X-ray range (354 and 348
of their 499 points). Kramers' law puts a shell's absorption at its edge; for
the 4d and 4f shells of a heavy atom the centrifugal barrier delays it to
several times the binding energy, so gold's 4f (84–88 eV) and 4d (335–353 eV)
shells absorb far too much just above their edges and too little further up,
and between about 150 and 500 eV the decrement is out by up to a factor of
three. The caveat reaches four times each 4d and 4f binding energy, which is
what covers platinum's 4f, and so includes 1 nm for both metals: there gold's
decrement is 6% low and its absorption about half the table's. Gold's deeper 3d
levels (M4,5, 2.2 keV) show no such error and are not flagged beyond their
edges.

Two of the caveat's settings were chosen by reading this comparison, so it is
not an independent test of them: the delayed band's reach of four times the
binding energy, and leaving the 3d levels out of it. The model's physics and
its facts were not fitted to anything.

These are measurements, not the accuracy the package claims; that wording is
still to be decided.
