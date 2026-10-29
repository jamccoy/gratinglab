# Groove metrology: from an AFM scan to a boundary profile

Implemented in [`metrology/`](../../src/gratinglab/metrology/). Units and sign
conventions follow [`conventions.md`](../conventions.md), §1 in particular. The
metrology package spans four unit regimes, and the table in that section is
normative.

An efficiency computed from a measured profile inherits every assumption made
while producing that profile, and those assumptions are not visible in the
resulting curve. This page records them so that anyone reading such a curve can
find them.

---

## 1. What the pipeline assumes

| Assumption | Consequence when it fails |
|---|---|
| The scan is of a periodic grating, and the period estimate is close enough for peak detection to find groove minima | Grooves are missed or spurious; the measured period is wrong, and with it the groove depth (§3) |
| Every groove in the scan is a sample of *the same* groove | Real variation across the scan is averaged into a shape that no individual groove has |
| Row-averaging the 2-D image to 1-D loses nothing that matters | True only if the grooves run parallel to the image rows. A rotated scan smears the facet |
| The averaged groove tiles seamlessly into a periodic boundary | Enforced rather than assumed (see §2), but enforcement is itself an operation on the data |
| The AFM tip is small compared to the features | The steep facet is the tip shape, not the grating, and the apex and trough are rounded off; this removes depth while leaving intact the mid-facet slope on which a line fit lands |
| The groove is faceted, so a fitted facet angle and the groove depth describe the same shape | False on a rounded groove, and then the two disagree. The solver uses the **depth** |

The last two assumptions are not mitigated anywhere in this package and have no
diagnostic. On the one measured grating in this repository both are false: see
`findings.md`, "The measured groove is rounded rather than faceted", where the fitted
facet angle (27.91° ± 2.13°) and the angle implied by the depth (20.33°) differ
by 7.6°.

The least expensive diagnostic available has not yet been implemented. It would
compute the blaze angle implied by the measured depth and compare it with the
fitted angle. The two agree on a sharp groove and diverge on a rounded one, and
the divergence is the same error that reaches the efficiency. Both numbers are
already available within the pipeline.

Such a diagnostic should be thresholded in degrees. The scatter of the fitted
angle is small on clean data, so a gap of a few tenths of a degree already
amounts to several sigma on a groove that is perfectly sharp: the synthetic
fixture shows 0.33°, which is 3.8σ and carries no meaning. The gap on the real
scan is 7.6°. The absolute size of the gap discriminates between the two cases,
and its size in sigma does not.

The diagnostic would also have limits. A rounded groove and a blunt tip produce
the same profile, so a divergence indicates that the groove is not a sharp
sawtooth without indicating that the tip is at fault. A *flat land*, a genuinely
unfaceted flat within the period that is common on real blazed gratings, also
produces a depth deficit and is the first alternative to exclude. It appears as
a spike at zero in the local-slope distribution, which rounding does not
produce.

---

## 2. What `normalize_profile` does, in order

The function performs four operations, in
[`boundary/average.py`](../../src/gratinglab/metrology/boundary/average.py).
They are kept in the order used by the original script, because changing the
order changes the answer.

1. **Endpoint flattening.** The straight line through the first and last sample
   is subtracted, so that the trace starts and ends at zero. This differs from
   the background flattening used for blaze-angle work, which removes a fitted
   background so that a facet can be measured. Substituting one for the other
   leaves a step at the period boundary, a defect that makes a PCGrate
   efficiency curve wrong while the file still appears correct.
2. **Shift and normalise.** The profile is referenced to the mean of the two
   half-minima, and both axes are then divided by the period.
3. **Roll the trough to the boundary.** The deepest point is placed at `x = 0`
   and `x = 1`, and the apex lies somewhere between them. Only `y` is rolled;
   `x_norm` remains a uniform ramp.
4. **Optional smoothing**, by default `uniform_filter1d(size=5, mode='wrap')`,
   after which the endpoints are forced back to exactly zero.

Handedness is preserved and never mirrored. The facet orientation matches
increasing scan-x in the source data. Whether that matches the sense the solver
expects is a question about the *mount* and lies outside this package; see
`conventions.md` §3, "Direction of the profile parameter".

### The stretch

Step 2 divides by the *averaged window width*, and the output of step 3 spans
exactly `[0, 1]` whatever that width was. If a groove near the scan edge narrows
the window that every groove shares, the profile is consequently stretched to
fill a period it never occupied. This effect is real and was found empirically:
`test_edge_rule_reduces_x_stretch` was written because of it, and
`ggp_min_half_width` mitigates it by discarding clipped grooves.

A stretched profile does not look obviously wrong in a plot. Its error appears in
the efficiency, because the facet angle it implies is wrong.

---

## 3. Why the period is load-bearing

The exported profile is dimensionless. Physical depth is `y_norm × period`, and
that depth enters the scalar phase term directly:

    phase = (2π / λ) · height · sin γ

so a fractional error in the period is a fractional error in the phase, at every
wavelength. The period therefore enters the result directly and cannot be
treated as metadata.

The period is measured as the mean spacing of detected groove centres, stored as
`period_nm` on `BoundaryProfile`. The initial `period_est` from settings is used
only to seed peak detection, and it is replaced as soon as there are two grooves
between which to measure.

A `.ggp` file cannot carry this number.
[`BoundaryProfile.to_problem`](../../src/gratinglab/metrology/boundary/pipeline.py)
exists so that a profile reaching a solver in-process retains it. The file-based
route requires the period to be supplied again from another source, and
`benchmarks/corpus.toml` serves that purpose.

---

## 4. What is measured, and what is not

| Quantity | Measured? | Where |
|---|---|---|
| Period | **yes**, from groove spacing | `BoundaryProfile.period_nm` |
| Groove depth, peak-to-valley | **yes**, as a fraction of the period | `metrics` |
| Blaze / anti-blaze facet angle | **yes**, by linear fit with a propagated uncertainty | `core.analysis.extract_blaze_angle`, the *other* branch of the pipeline |
| Groove-to-groove form spread | **yes** | `y_std_nm`, the ±1σ band |
| Max sidewall angle | derived, not fitted | `metrics['max_angle_deg']`, the steepest local slope of the normalised profile |
| Apex angle | **no** | not computed anywhere |
| Groove residual about the average | **yes**, after sub-pixel registration | `metrics['groove_residual_rms_nm']`, per groove in `metrics['groove_residual_per_groove_nm']` |
| **RMS surface roughness** (the Névot–Croce σ) | **no** | (none; the residual above is related, and is not it) |

### What an electron micrograph measures instead

The table above is about an AFM scan. A cleaved cross-section, read by
`metrology/io/sem.py`, is subject to none of the same limits and to several of
its own, and the two are complementary rather than ranked:

| Quantity | AFM scan | SEM cross-section |
|---|---|---|
| Groove depth | yes, absolutely calibrated | yes, but only of the grooves the cleave happens to cross |
| Facet or sidewall angle | yes, with a propagated uncertainty, subject to tip dilation | yes, free of any tip, subject to the plane of the cleave |
| Width of a crest or a land | dilated by the tip, which rounds exactly this feature | yes, directly |
| Undercut | never; a tip cannot reach beneath an overhang | yes |
| Averaging over many grooves | yes, hundreds per scan | no; one cleave, a few dozen grooves, destructively |
| Absolute height scale | from the Z sensor | from the pixel size, and foreshortened if the stage was tilted |

The last row is the SEM's counterpart of the tip problem, and it is handled the
same way: `SemImage.foreshortening` reports `cos` of the stage tilt when the
instrument's own tilt correction was off, and the reader declines to apply it,
because which image direction lies across the tilt axis depends on how the
specimen was loaded and the file does not record that.

The destructive row is why a cross-section is an arbiter and not a routine
input. It settles a question about a sacrificial piece of a wafer; it cannot be
run on the grating that is about to fly.

### Roughness is not computed, and related quantities are not substitutes

`Problem.roughness` denotes high-spatial-frequency surface roughness in nm, the
quantity with which Névot–Croce damps a Fresnel coefficient. Nothing in this
package measures that σ, and `to_problem` passes 0.0, meaning "not supplied",
unless the caller passes a number.

Two quantities appear related but do not measure it:

- `y_std_nm` is the spread *between grooves*: form variation at the spatial
  scale of a whole groove. Névot–Croce concerns structure far below that scale.
- `metrics['rms_slope']` is a slope statistic of the normalised profile, in
  dimensionless units. It is not a length.

Passing either into `Problem.roughness` would produce a confidently wrong
efficiency where an obviously missing one is expected, so `to_problem` refuses
to do so, and a test pins the refusal.

### The groove residual: what it is

`average.groove_residuals` takes the grooves that went into the average and
measures how far each departs from it: the part of the profile that does not
repeat with the period. It is in the metrics, the metrics sidecar, the
boundary `.json` and the Boundary tab, as `groove_residual_rms_nm`.

A plain difference against the average does not measure that. Each groove is
centred on the sample nearest its minimum, so the set is misregistered by up
to half a pixel, and on a 70° facet half a pixel of shift is several
nanometres of height. On the synthetic sawtooth, with no noise in it at all,
the plain difference reads 2.1 nm. So each groove is first registered to the
average by a sub-pixel shift and a height offset, the average is recomputed
from the registered grooves, and the residual is taken at each groove's own
samples, over the average's aligned window less one pixel at each end. It is
scaled by √(N/(N−1)) because the average contains the groove it is compared
with. The exported average is unchanged: it is still the unregistered mean.

What the number contains, measured on synthetic scans of the fixture's sawtooth
(`tools/metrology/make_synthetic_scan.py`), the floor removed in quadrature:

| Term | How it enters | Measured |
|---|---|---|
| Non-periodic structure coherent along the groove | in full | 1-D noise of σ added identically to every scan line reads 0.93 to 0.94 σ, over eight seeds each at 0.3, 0.5 and 0.8 nm |
| AFM noise, and anything that varies along the groove | divided by √(scan lines), because the profile is the mean of the lines | per-pixel noise of 4 nm on 64 lines reads about 0.46 nm |
| Corners sharper than a pixel | a sampling floor | 0.21 nm on the sharp synthetic at 5.2 nm per pixel, 0.03 nm at 2.6 |

The 6% shortfall on noise is most likely the shift and offset fitting a little
of it away, which was not separately measured. The floor is the spline failing to place a corner it never sampled. A real,
tip-rounded groove has no such corner, so its floor is lower, though the tip
also hides whatever is finer than its radius.

### Why it is not passed as the Névot–Croce σ by default

It is passed only when the caller asks:
`profile.to_problem(roughness=profile.groove_residual_rms_nm)`. Four reasons,
each of which makes it a different number from the σ the model wants:

1. **It is measured on the row-averaged profile.** Roughness that varies
   along a groove averages away over the scan lines, and for isotropic
   roughness with a correlation length well under the scan height that is
   nearly all of it. What survives is the part coherent for the height
   of the scan, so the number is a lower bound on the surface's roughness in
   its band, not an estimate of it.
2. **Its band is not the σ's band.** An AFM scan resolves nothing finer than
   the pixel pitch or the tip, a few nanometres. The roughness that damps a
   soft X-ray reflection extends down to wavelengths near the X-ray
   wavelength itself. At the other end the residual includes groove-scale form
   error, whose correlation length is about a period, while Névot–Croce
   assumes a short correlation length.
3. **Form error acts on the orders, not the facet.** Groove-to-groove
   variation removes light from each order by scattering it into the space
   between orders, a Debye–Waller factor on that order's own momentum
   transfer. `Problem.roughness` damps the Fresnel coefficient of the facet.
   Folding one into the other applies the wrong factor to a real effect.
4. **It carries instrument terms**: the noise and sampling floor in the table
   above.

So an absolute efficiency from a measured profile still has no roughness
damping by default, and is therefore an overestimate. Passing the residual
gives a damping that is too small for reason 1 and 2 and of the wrong form for
reason 3. Which error is worse depends on the grating and the energy, which is
why the choice is left to the caller and written in their code rather than
made silently here.

A per-scan-line residual, taken before the lines are averaged, would recover
reason 1 at the cost of carrying the AFM's full noise, which a single scan
cannot separate from roughness. It is not computed.

---

## 4b. Tip correction: what erosion can and cannot give back

An AFM image is the surface *dilated* by the tip. At every pixel the tip
descends until first contact, and the recorded height is the position at which
the apex stopped, which lies above the true surface wherever contact occurred on
the flank. Villarrubia (1997) gives the exact morphological inverse, which this
pipeline implements in `core/tip.py`. It is enabled through
`AnalysisSettings.tip_correction = 'erosion'` with a parametric tip: a spherical
apex of `tip_radius_nm` on two straight flanks, `tip_half_angle_deg` on the -x
side and `tip_half_angle_right_deg` on the +x side (both measured from the axis;
the second defaults to the first, giving a symmetric cone). The two flanks are
the cross-section of a pyramid in the scan plane, so they depend on how the
cantilever lay against the fast axis. A ScanAsyst-Air, for example, shows its
17.5° sides across the cantilever and its 15° front and 25° back along it. The
correction runs on the 2-D array immediately after image flattening, in both
branches, so that the blaze fit and the boundary profile see the same surface.

With unequal flanks the erosion uses the reflected tip,
`rec(y) = min_u [img(y - u) + tip(u)]`. A symmetric tip conceals the
reflection, which is why it is stated here.

Three properties of the correction determine how its output may be used:

- **Erosion gives a least upper bound.** Wherever the apex made contact, the
  bound coincides with the surface and the recovery is exact: to machine
  precision on clean data, and to the noise floor on the fixture
  (`tests/metrology/test_tip.py`). Elsewhere the tip physically could not reach
  the surface, and no algorithm can recover what lies there. The trough wedge
  between facets is the canonical case: a sphere of radius R stands off a corner
  by a distance that scales with R, and that depth is absent from the data.
- **The certainty map is a primary output.** A pixel is certain when it is the
  unique contact point for some image position. The certain fraction is stored
  in `metrics['tip_certain_fraction']` and in the metrics sidecar, and the
  uncertain remainder must be read as an upper bound on the height, not as the
  height itself. A corrected profile reported without that number would create
  the false confidence that the correction is intended to remove.
- **A facet steeper than the flank facing it cannot be recovered at all.** A
  flank rises at `90 − half_angle` degrees, and a wall rising towards +x is
  followed by the +x flank. An 18° tip clears a 70.5° anti-blaze facet by 1.5°.
  The 25° back face of a ScanAsyst-Air, a 65° flank, does not clear it at all,
  so along the cantilever the orientation decides whether that facet exists in
  the data. At these geometries it is this near-parallelism, and not the
  nanometre apex, that limits blazed-grating metrology.

The correction is off by default because it changes measured numbers, depth in
particular, and a corrected depth and an uncorrected depth are different
measurements. When it is enabled, the settings and the certain fraction are
recorded with the output (in the metrics, the sidecar, and the log line of the
CLI).

The correction is validated against the dilated fixture: the same synthetic
sawtooth as the sharp fixture, with the identical noise field, imaged through a
known R = 20 nm tip by `dilate`. The truth is committed next to the image, so
recovery is asserted pointwise. The result of applying it to the real TASTE scan
is recorded in `findings.md` ("A nominal tip does not explain the rounded
groove"): erosion with the nominal probe recovers nothing, and the rounding
would require an apex of ~80 nm to be a tip artefact at all.

---

## 5. Two branches, one front-end

Loading, image flattening, row-averaging and groove detection are shared. After
these steps the pipeline divides into two branches:

- **Blaze angle.** `core.analysis.extract_blaze_angle` fits a line to each facet
  independently, per groove, with row-group statistics and an ICC correction for
  the correlation between rows of the same scan.
- **Boundary profile.** `boundary.build_boundary_profile` averages the grooves
  into one and normalises it.

The two branches were once separate scripts with their own copies of the
front-end, and these copies drifted apart until a bug fixed in one persisted in
the other. The front-end is shared for that reason.

The two branches do not exchange results. `build_boundary_profile` does not run
the facet fit, so `to_problem(blaze_angle=...)` takes the angle as an argument
instead of obtaining it. Connecting them would require substantial work and has
deliberately not been done: the fit has settings of its own (facet trimming in
particular), and adopting their defaults silently would conceal a choice that
changes the answer.

---

## 6. From an SEM cross-section to a groove

Implemented in
[`boundary/cross_section.py`](../../src/gratinglab/metrology/boundary/cross_section.py)
as `trace_cross_section`, which takes a `SemImage` (or an array and a pixel
size) and returns the traced boundary, the averaged groove as a
`FromProfileData`, and a fitted `Trapezoidal` with uncertainties. It is the
third row of the table in §4, "What an electron micrograph measures instead",
turned into a number.

### What it does, in order

1. **Levels.** The median of the region's top rows is the vacuum level and of
   its bottom rows the material level. The iso-level sits halfway between
   (`level_fraction`), or is given outright (`level`). A region whose two
   levels differ by less than four times the vacuum noise is refused rather
   than traced. Material darker than vacuum is handled by negating the image.
2. **Edge, per column.** After a light Gaussian smoothing (`smoothing_px`,
   0.5 by default), each column is scanned downward from the vacuum and the
   *first* pixel at or above the level, followed by at least two more, is the
   surface. Its row is interpolated linearly between that pixel and the one
   above it. The first crossing is the point of the method: an
   in-lens or secondary-electron image puts a bright rim just inside every
   surface, and the bulk beneath it may be darker than the rim and even darker
   than the level. The outer flank of the rim is the surface; the inner flank
   is not, and a contour tracer would follow both. The run of three rejects an
   isolated charging speck in the vacuum.
3. **Edge, across walls.** Where two neighbouring columns cross more than a
   pixel apart, each row between them is crossed horizontally too, so a steep
   wall is located to a fraction of a pixel across as well as down. Below its
   crossing a column is treated as solid, but only a pixel reading *below* the
   level is raised to the material level. Raising every pixel there also
   lifted the partly covered pixel at the wall and moved each wall 0.17 nm
   into the vacuum on a clean render; that was found by the synthetic tests and
   is why the rule is as narrow as it is.
4. **Grooves.** A first period comes from the strongest spatial frequency of
   the traced height. Each ridge's centre is the midpoint of its two
   half-height crossings, found with hysteresis at a quarter and three quarters
   of the depth so noise on a wall cannot count as a second wall. A straight
   line through the centres against their index gives the measured period, its
   standard error and the placement scatter. A tilted cleave edge is levelled
   by a line through the *per-groove mean heights*: a least-squares line
   through the whole trace is biased by the groove shape (the covariance of a
   periodic function with `x` over whole periods is not zero), and the mean
   over one period is not.
5. **Fit.** Every groove with a whole period about its own centre is folded
   into one period, ridge on `u = 0.5`, heights in periods, and a
   `Trapezoidal` is fitted to all the folded points at once by **orthogonal**
   distance (`scipy.optimize.least_squares`), with a horizontal and a vertical
   offset. Orthogonal, because a vertical residual on a 72° wall is three times
   the real one and would let the walls outvote the flats. The parameters are
   chosen so every point inside the bounds is a valid groove: crest width, the
   base width as a share of the room left, so the floor never closes, and the
   radius as a share of the largest that fits.
6. **Uncertainties.** The same fit is run on each groove alone, and every
   quoted error is the standard deviation across grooves over the square root
   of their number. That carries the real variation of the etch, not only the
   pixel noise. Depth and radius in nm also carry the period's error. With a
   single complete groove the least-squares covariance is used instead and
   `TrapezoidFit.uncertainty` says `"covariance"`.
7. **The averaged groove.** Starting from the fitted shape, the orthogonal
   residuals are binned along it by the arc position of their foot (one bin
   per pixel of arc), the curve is moved along its normal by each bin's mean,
   and the projection is repeated on the result until a pass moves nothing by
   a hundredth of a pixel. One pass is enough where the trapezoid fits. It is
   not where it does not: an asymmetric groove came back from one pass with its
   steep wall twice as wide as rendered, and the repeated passes restore it.
   The averaged groove is therefore model-free, and its distance from the
   fitted trapezoid (`mean_profile_rms_nm`, `mean_profile_max_nm`) is the
   honest statement of how well a trapezoid describes the grating.

### Inputs that are inputs, not guesses

| Input | Why the file cannot decide it |
|---|---|
| `foreshortened` | `SemImage.foreshortening` is `cos` of the stage tilt when the tilt correction was off, and which image axis it compresses depends on how the piece was mounted. Any image with a factor below one is refused until this is `"vertical"`, `"horizontal"` or `"neither"`, and the factor is then undone on that axis alone. Naming the wrong axis is not silent: it shows as a period wrong by `1 / cos`, which the tests pin |
| `material_below` | which side of the cleave edge is solid |
| `image_x_along_dispersion` | which way the columns run against `d̂`. `to_profile` uses `FromProfileData.along_dispersion` when it is true (the default, matching the AFM pipeline's `x_norm`), so `t = 1 - u` as `conventions.md` requires. A symmetric groove does not care; an asymmetric one does, and the tests use one |
| `period_nm` | optional. A period known from diffraction is better than a few dozen grooves; given, it normalises the profile, and the lattice period is still reported beside it |
| `region` | the levels are read from the region's top and bottom rows, so it must hold vacuum above and solid below and nothing else |

### Validation

No real micrograph is in the repository, and none may be. The tests render
cross-sections of a known `Trapezoidal` (`tests/metrology/xsec_synth.py`):
pixel coverage from an 8×8 supersample, vacuum 40 and material 110 on the
8-bit scale, an optional bright rim decaying over 4 nm inside the solid, a
Gaussian beam, Gaussian noise and 8-bit quantisation, at the scale of
`ogre_xsec/xsec_050` in `findings.md` (2.233 nm pixels, a 163.3 nm period,
8 whole grooves in 640 columns). The truth is that table's groove: 68 nm deep,
a crest 0.30 of the period, a 72° wall and a 6 nm radius. Five noise
realisations per row, measured 2026-10-29 with the defaults:

| case | depth error | wall position, per side | angle error | radius, nm | reported σ depth / angle | scatter between realisations |
|---|---|---|---|---|---|---|
| clean, no beam | +0.10 nm | −0.01 nm | −0.01° | 6.7 | noise-free | noise-free |
| 0.7 px beam | +0.10 | −0.02 | −0.03 | 7.7 | noise-free | noise-free |
| + noise σ = 6 | +0.11 | −0.02 | −0.04 | 7.6 | 0.028 nm / 0.048° | 0.028 nm / 0.051° |
| + rim of 40 | −0.15 | +0.56 | −0.05 | 7.2 | 0.024 / 0.038 | 0.025 / 0.035 |
| + rim of 80 | −0.27 | +0.93 | −0.05 | 7.1 | 0.021 / 0.032 | 0.016 / 0.035 |
| + rim 80, edge tilted 0.03 | −0.04 | +0.87 | −0.02 | 7.7 | 0.059 / 0.033 | 0.031 / 0.039 |

The reported errors match the scatter between independent images, which is
what they claim to be. They are **statistical only**: the two systematic rows
below are larger than any of them, and nothing in the result adds those in.

A vertical wall (the resist line, 223.5 nm deep, 90°) is found at 89.9° and
223.5 nm; a 45° wall at 45 ± 0.5°.

### What moves the answer, and is not corrected

**The beam rounds every corner.** The level set of a blurred edge moves with
its curvature, so a corner traced through a beam is rounder than the corner
was, and the fitted radius reads high. Measured with no noise, true radius
against beam σ:

| true radius | no beam | 0.7 px (1.6 nm) | 2 px (4.5 nm) |
|---|---|---|---|
| 0 | 2.9 | 4.8 | 11.4 |
| 6 | 6.7 | 7.7 | 12.7 |
| 12 | 12.5 | 12.9 | 16.1 |
| 20 | 20.3 | 20.6 | 22.4 |

Close to `R_fit² = R² + R₀²`, with `R₀` about two and a half times the beam σ.
So a fitted radius is an upper bound, and one below about `2.5 σ` of the beam is
not resolved at all. The depth, the angle and the widths are not affected:
the walls are straight and a straight edge does not move.

**A bright rim moves each wall into the vacuum** by 0.5 to 0.9 nm for a rim
0.6 to 1.1 times the material contrast, the ridge reads that much wider on each
side and the depth a few tenths of a nanometre shallower. The first crossing
stops the rim's inner flank being traced; it cannot stop the beam spreading
the rim's light outward. A threshold method has no way to tell that light from
the specimen without a model of the detector, and none is attempted.
Raising `level_fraction` pulls the edge back in, but on an image without a rim
it moves the edge by the same amount the other way.

### What it cannot do yet

- **Undercut walls.** The trace is a height function (a column has one first
  crossing), so an overhang is traced as a vertical wall. The averaged groove
  drops any point that steps back in `u` for that reason, so `to_profile()`
  is always single-valued. Tracing an undercut wants a true contour follower,
  which is the one thing a column scan cannot become.
- **A protective cap.** A cross-section prepared with a deposited cap over the
  grooves shows grating against cap, not grating against vacuum, and the
  outermost-surface rule will trace the cap. `level` and a `region` cropped
  below the cap may get through it; nothing here has been tested on one.
- **A tilted view of the top surface.** An image holding both the cleave face
  and the grating's top surface beyond it must be cropped to the face.
- **Rotation beyond a shear.** A tilted cleave edge is levelled, which is a
  shear; a rotation of θ also shortens the period by `cos θ`, which is second
  order and not corrected.
- **The real micrographs.** `findings.md`'s three OGRE cross-sections are on a
  machine this was written without access to, so the method has not yet been
  run on them; its numbers there were read from brightness thresholds by hand.
