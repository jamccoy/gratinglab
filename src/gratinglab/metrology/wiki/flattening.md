# Flattening

AFM data arrives with a background: the sample sits at an angle under the
scanner, and scan lines drift in Z relative to one another. Removing that is
"flattening", and this software does it at **two stages** that are easy to
confuse.

| | Acts on | When | Changes the blaze angle? |
|---|---|---|---|
| **Image flattening** | the 2-D scan | before rows are averaged | **No** (for the methods offered) |
| **Profile flattening** | one averaged 1-D trace | per row group | **Yes — about 0.49°** |

Both are set in the **Import** tab, where the profile preview draws the
background about to be removed.

## Image flattening (2-D)

Acts on the image as it came off the instrument, correcting scan lines relative
to each other.

| Method | What it does |
|---|---|
| `none` | Leave the image alone. The behaviour before August 2026. |
| `plane` | Subtract one tilted plane fitted to the whole image. Corrects the sample sitting at an angle. |
| `align_rows` | Subtract each row's median, putting every scan line on a common level. **The default.** |

`align_rows` uses the median rather than the mean because a grating profile is
not symmetric about its own mean — a row catching more groove bottom than land
would be dragged by a mean where the median barely moves.

The Import tab reports the **row-offset spread** (the standard deviation of the
per-row medians) before and after, so the choice can be made on evidence. Across
this project's data that spread runs 0.44–2.74 nm, against groove depths near
100 nm.

### Why it doesn't change the answer

This is worth stating plainly, because it looks like a bug otherwise.

On `20250820_280C_00004.txt` — the scan with the worst row-offset spread in the
dataset — `none`, `plane` and `align_rows` all give **30.8524°**, agreeing to
within 5×10⁻¹⁵. Across all eight compare-mode samples, switching `align_rows` on
moved the mean by **0.0000°**.

The reason is structural. Rows are averaged into a profile, and profile
flattening then fits and removes a background from *that*. Subtracting a per-row
constant shifts a band's average by a constant; subtracting a plane adds a
constant and a linear ramp. Every profile-flattening method fits at least a
first-order polynomial, so it removes exactly those terms again. What you take
out here comes out there regardless.

That is why `align_rows` could become the default without a single stored number
moving, and it is guarded by a test.

### So what is it for

- **Seeing the data.** A raw image with sample tilt renders as a gradient with
  the grating barely visible. Levelled, it looks like a grating.
- **Methods that aren't affine.** A per-row *quadratic* changes the shape within
  each row rather than just its offset or tilt, and does survive the averaging —
  measured at 0.0003°. Small, but not zero. The registry in
  `core/image_flatten.py` is where such a method would go.
- **Being explicit.** The step now happens somewhere named and recorded, rather
  than in Gwyddion before export where nothing captured that it had happened.

## Profile flattening (1-D)

Acts on each row group's averaged height trace, immediately before grooves are
detected and facets fitted. **This is the choice that matters.**

| Method | What it fits the background to |
|---|---|
| `linear` | A straight line through the whole profile |
| `polynomial` | A polynomial of `FLATTEN_POLY_ORDER`, optionally ignoring the edges |
| `groove_peaks` | A polynomial through the detected peaks only |
| `level_grooves` | A polynomial through chosen features — `peaks` (the lands), `troughs` (groove bottoms), or `both`. **The default.** |

On the same scan:

| Method | Mean angle | N |
|---|---|---|
| `level_grooves` | 30.8524° | 102 |
| `linear` | 30.9619° | 102 |
| `polynomial` | 31.1098° | 105 |
| `groove_peaks` | 31.3466° | 106 |

A **0.49° spread** — around 250× the image-flattening effect, and comparable to
the differences between samples this software exists to detect. It also changes
how many grooves are found, because the flattened profile is what groove
detection runs on.

`level_grooves` is the default because levelling on the lands puts the reference
where the surface is genuinely flat, rather than fitting through a profile whose
own curvature is the thing being measured.

### Judge it by eye

The Import tab draws three curves: the averaged profile, the background about to
be removed, and the result. A background that follows the grating rather than the
underlying tilt is removing signal, and that is visible immediately in a way a
number in a config file never is.

Since this parameter moves the answer, the value used is written into every
summary file alongside the result. What matters is that the same choice is
applied across every sample being compared.
