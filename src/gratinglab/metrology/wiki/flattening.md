# Flattening

AFM data arrives with a background, because the sample sits at an angle under
the scanner and the scan lines drift in Z relative to one another. Removing this
background is called *flattening*. This software flattens at two separate
stages, which are easily confused.

| | Acts on | When | Changes the blaze angle? |
|---|---|---|---|
| **Image flattening** | the 2-D scan | before rows are averaged | **No** (for the methods offered) |
| **Profile flattening** | one averaged 1-D trace | per row group | **Yes, about 0.49°** |

Both are set in the Import tab, where the profile preview draws the background
that is about to be removed.

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

The reason is structural. Rows are averaged into a profile, and profile
flattening then fits and removes a background from that profile. Subtracting a
per-row constant shifts a band's average by a constant, and subtracting a plane
adds a constant and a linear ramp. Every profile-flattening method fits at least
a first-order polynomial, so it removes exactly those terms again, and whatever
is removed at the image stage is removed at the profile stage in any case.

For this reason `align_rows` could be made the default without any stored number
changing, and this is guarded by a test.

### So what is it for

- **Seeing the data.** A raw image with sample tilt renders as a gradient in
  which the grating is barely visible. After levelling, it is recognisable as a
  grating.
- **Methods that aren't affine.** A per-row quadratic changes the shape within
  each row as well as its offset or tilt, and its effect does survive the
  averaging: it was measured at 0.0003°, which is small but not zero. The
  registry in `core/image_flatten.py` is where such a method would be added.
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
| `level_grooves` | A polynomial through chosen features: `peaks` (the lands), `troughs` (groove bottoms), or `both`. **The default.** |

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

`level_grooves` is the default because levelling on the lands places the
reference where the surface is genuinely flat. The alternative is to fit through
a profile whose own curvature is the quantity being measured.

### Judge it by eye

The Import tab draws three curves: the averaged profile, the background that is
about to be removed, and the result. A background that follows the grating
instead of the underlying tilt is removing signal, and this is immediately
visible in the plot in a way that a number in a configuration file is not.

Because this parameter affects the answer, the value used is written into every
summary file alongside the result. What matters for a comparison is that the
same choice is applied to every sample being compared.
