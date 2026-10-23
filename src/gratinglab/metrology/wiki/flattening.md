# Flattening

AFM data arrives with a background, because the sample sits at an angle under
the scanner and the scan lines drift in Z relative to one another. Removing this
background is called *flattening*. This software flattens at three separate
points, which are easily confused.

| | Acts on | When | Changes the blaze angle? |
|---|---|---|---|
| **Image flattening** | the 2-D scan | before rows are averaged | **No** (for the methods offered) |
| **Profile flattening** | one averaged 1-D trace | per row group | **Yes, about 0.49°** |
| **Endpoint flattening** | the normalised groove | only on PCGrate export | n/a — a different job |

The first two are set in the Import tab, where the profile preview draws the
background that is about to be removed. The third is not a setting: it is part
of writing a `.ggp` file, and it is described at the end of this page.

## Why a background has to be removed at all

A blaze angle is a facet slope *measured against the substrate plane*. The
substrate plane is not in the data — the scanner records height above whatever
plane the sample happened to be sitting in — so flattening is what supplies the
reference the angle is quoted against. An error of δ in that reference is an
error of δ in every facet angle on the scan, in the same direction. This is why
the choice of profile-flattening method matters as much as the fit that follows
it, and why it is recorded in every summary file.

## Image flattening (2-D)

Image flattening acts on the image as it came off the instrument and corrects
the scan lines relative to each other.

| Method | What it does |
|---|---|
| `none` | Leave the image alone. The behaviour before August 2026. |
| `plane` | Subtract one tilted plane fitted to the whole image. Corrects the sample sitting at an angle. |
| `align_rows` | Subtract each row's median, putting every scan line on a common level. **The default.** |

`align_rows` uses the median because a grating profile is not symmetric about
its own mean: a row that catches more groove bottom than land would shift the
mean noticeably, while the median barely moves.

The Import tab reports the *row-offset spread*, defined as the standard
deviation of the per-row medians, before and after flattening, so that the
choice of method can be based on evidence. Across this project's data the spread
ranges over 0.44–2.74 nm, against groove depths near 100 nm.

### Why it doesn't change the answer

This result is stated explicitly because it could otherwise be mistaken for a
bug.

On `20250820_280C_00004.txt`, the scan with the largest row-offset spread in the
dataset, `none`, `plane` and `align_rows` all give 30.8524°, agreeing to within
5×10⁻¹⁵. Across all eight compare-mode samples, enabling `align_rows` changed
the mean by 0.0000°.

The reason is structural, and it is worth following through, because it is what
makes the default free.

1. Rows are averaged into a profile. `align_rows` removes a constant from each
   row, so after averaging a band it has removed **one constant**. `plane`
   removes `a + bx + cy`, so after averaging over rows it has removed **a
   constant and a linear ramp in x**. Both are polynomials in x of degree at
   most 1.
2. Profile flattening then fits a polynomial to that averaged profile and
   subtracts it. A least-squares polynomial fit is linear in the data, and it
   reproduces exactly any polynomial up to its own order. Every profile-
   flattening method fits at least a first-order polynomial.
3. So adding `a + bx` to the profile adds exactly `a + bx` to the fitted
   background, and the difference — the flattened profile — is unchanged.

Whatever the image stage removes, the profile stage would have removed anyway.
For this reason `align_rows` could be made the default without any stored number
changing, and this is guarded by `test_affine_methods_give_identical_blaze_angles`
in `tests/metrology/test_image_flatten.py`. That test asserts a tolerance rather
than exact equality on purpose: subtracting medians perturbs the last bits, and
the averaging that follows is not associative, so the results differ by around
5×10⁻¹⁵ degrees against a sigma near 2.

### Where that argument stops

Two limits are worth knowing, because the guarantee is narrower than "flattening
is harmless".

**It covers affine methods only.** A per-row *quadratic* fits coefficients to
each row's own data and changes the shape within the row, not merely its offset
and tilt. That survives the averaging: it was measured at 0.0003°, which is
small but not zero. The registry in `core/image_flatten.py` is where such a
method would be added, and it would need its own measurement rather than
inheriting this one.

**Step 3 assumes the same features are found.** `level_grooves` and
`groove_peaks` fit through detected features, and detection uses topographic
prominence. Prominence is unchanged by adding a constant, so `align_rows` cannot
disturb it. It is *not* unchanged by adding a tilt, so with `plane` a large
enough sample tilt could in principle shift which peaks clear the prominence
threshold, and the two methods would then diverge for that reason rather than
through the arithmetic above. On this project's data it does not happen — hence
the 5×10⁻¹⁵ agreement — but on a markedly more tilted scan it is the thing to
check first.

### So what is it for

- **Seeing the data.** A raw image with sample tilt renders as a gradient in
  which the grating is barely visible. After levelling, it is recognisable as a
  grating.
- **Methods that aren't affine.** As above: the registry is the place for one,
  and the invariance result does not extend to it.
- **Being explicit.** The step is now performed in a named and recorded place.
  Previously it was done in Gwyddion before export, and nothing recorded that it
  had been done.

## Profile flattening (1-D)

Profile flattening acts on each row group's averaged height trace, immediately
before grooves are detected and facets are fitted. This is the flattening choice
that affects the result.

| Method | What it fits the background to |
|---|---|
| `linear` | A straight line through the whole profile |
| `polynomial` | A polynomial of `FLATTEN_POLY_ORDER`, optionally ignoring the edges |
| `groove_peaks` | A polynomial through the detected peaks only |
| `level_grooves` | A polynomial through chosen features. **The default.** |

The settings that go with it, with the values they ship as:

| Setting | Default | Meaning |
|---|---|---|
| `FLATTEN_METHOD` | `level_grooves` | which of the four methods above |
| `FLATTEN_FEATURE` | `peaks` | for `level_grooves`: `peaks` (the lands), `troughs` (groove bottoms), or `both` |
| `FLATTEN_POLY_ORDER` | `2` | order of the fitted background: 1 linear, 2 quadratic, 3 cubic |
| `FLATTEN_EXCLUDE_EDGES` | `0.05` | fraction of the trace ignored at each end when fitting |

Note that the default *feature* is `peaks`, not `both`: the background is fitted
through the lands alone unless asked otherwise.

On the same scan:

| Method | Mean angle | N |
|---|---|---|
| `level_grooves` | 30.8524° | 102 |
| `linear` | 30.9619° | 102 |
| `polynomial` | 31.1098° | 105 |
| `groove_peaks` | 31.3466° | 106 |

The methods span 0.49°, which is around 250× the image-flattening effect and
comparable to the differences between samples that this software is intended to
detect. The choice also changes how many grooves are found, because groove
detection runs on the flattened profile.

### Why the lands are the reference

`level_grooves` is the default because levelling on the lands places the
reference where the surface is genuinely flat. The lands are the original wafer
surface, the part the etch did not remove, so they are parallel to the substrate
by construction and are the one part of the profile whose height is not a
measurement of the thing being measured.

The alternative is to fit through a profile whose own asymmetry and curvature
are the quantity being measured. A sawtooth is not symmetric about its own mean
line, so a fit through the whole trace is pulled by the facet, and the
background absorbs part of the slope the facet fit is about to report. That is
the mechanism behind the 0.49° in the table: `linear` and `polynomial` are not
noisier than `level_grooves`, they are *biased*, and in a direction set by the
groove shape rather than by anything about the instrument.

`groove_peaks` fits through peaks too, but it finds them with a bare separation
criterion (one twentieth of the trace) and no prominence requirement, so it
accepts noise maxima on a facet as readily as lands. It is kept because it was
the historical method and its number is needed to interpret older output.

### How the features are found

For `level_grooves`, each requested feature set is detected with `find_peaks`
on the profile — troughs by running the same detection on the negated trace —
under two criteria:

- **Separation** of at least half a period, `0.5 · period_nm / dx`, so that one
  groove cannot contribute two lands.
- **Prominence** of at least 20% of the profile's full height range, which is
  what excludes noise maxima part-way up a facet.

The surviving points, from one or both sets, are fitted with a polynomial of
`FLATTEN_POLY_ORDER` and that polynomial is the background.

### The period comes first

`level_grooves` requires an estimated period: it raises rather than guessing if
one is not supplied, because the separation criterion above is expressed in
periods. So the order is period estimate → flattening → groove detection →
measured period, and the period used to level the profile is not independent of
the period the run reports.

This is not vicious — the separation criterion is a loose gate at half a period,
and a period estimate wrong by a few percent does not change which maxima clear
it — but it is a real dependency, and it is the reason a badly wrong period
estimate degrades the flattening rather than only the detection.

### When detection fails

If fewer than two features survive, `level_grooves` falls back to `linear` and
says so; `groove_peaks` falls back to `polynomial` below three peaks. Either
fallback means the background for that row group was fitted by a method other
than the one requested, and from the table above that is worth up to a few
tenths of a degree.

The fallback is announced on standard output rather than recorded in the
summary, so in a batch run over many files it is easy to miss. If a row group's
angle is an outlier, the console is the place to check whether its background
was fitted the way the others' were.

### Judge it by eye

The Import tab draws three curves: the averaged profile, the background that is
about to be removed, and the result. A background that follows the grating
instead of the underlying tilt is removing signal, and this is immediately
visible in the plot in a way that a number in a configuration file is not.

Because this parameter affects the answer, the value used is written into every
summary file alongside the result. What matters for a comparison is that the
same choice is applied to every sample being compared.

## Endpoint flattening (export only)

The boundary profile written for PCGrate is levelled a third way, and this one
is not a choice. The straight line through the first and last sample is
subtracted, so that the trace starts and ends at exactly zero.

The reason is that a `.ggp` profile is one period, repeated by the solver. A
background fitted the way the blaze-angle path fits one leaves the two ends at
slightly different heights, which on repetition is a step at every period
boundary. That step makes the computed efficiency curve wrong while the file
still looks entirely reasonable, which is the worst combination. Substituting
one kind of flattening for the other here is a known way to produce a confidently
wrong answer, and the boundary-export path keeps them apart deliberately.
