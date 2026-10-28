# Changelog

Notable changes, newest first. The project is pre-1.0, so a minor version may
change an interface; any such change is stated here when it is made. The
versioning policy, including what 0.2 and 1.0 will mean, is in
[`docs/roadmap.md`](docs/roadmap.md#versioning).

The format loosely follows [Keep a Changelog](https://keepachangelog.com). The
version is read from `gratinglab.__version__`, and a release commit moves the
entries under "Unreleased" into a section named for the new version;
`tests/test_changelog.py` fails if it does not.

## Unreleased

### Added

- The rigorous solvers, `integral` and `rcwa`, cap BLAS at one thread for
  the length of a solve and restore the process's setting on the way out,
  however the solve ends. OpenBLAS's thread-per-core default is slower at
  these matrix sizes, and on a busy machine much slower: RCWA paid 20 to 135x
  and the integral solver 1.1 to 3.7x. `blas_threads=None` opts out. The cap
  is process-wide, as a BLAS pool is, and reference-counted so overlapping
  solves on different threads share it. New dependency: `threadpoolctl`
  (pure Python, BSD-3), the maintained way to reach NumPy's and SciPy's
  separate OpenBLAS copies.
- `workers=N` on both rigorous solvers solves wavelengths on `N` threads.
  The default, 1, is the serial loop. The scan is identical to the serial
  one, and the progress contract holds: every call on the caller's thread,
  `(0, n)` first, `done` stepping by one, and `SolveCancelled` re-raised
  unchanged once the rows already running finish, with nothing left running.
  The integral solver gives each thread its own kernel memo over the scan's
  shared geometry (`OperatorAssembler.fork`).
- RCWA solves an in-plane mount as two decoupled `N`-square systems, TE and
  TM, instead of one coupled `2N` system, and skips a polarization nothing
  excites. The backend also forms the flux at the top plane only, which is
  all `absorption` needs; `solve_stack(planes="all")`, still the core's
  default, gives every plane. Both agree with the paths they skip to 1e-12
  (`tests/test_solver_execution.py`).
- `benchmarks/solver_speed.py` times each case in a fresh interpreter under
  inherited and single-thread BLAS. Before and after, on MC1 (16 cores) with
  another session's solver running alongside, best of 3, one wavelength
  unless noted:

  | case | before, default BLAS | before, 1 thread | after |
  |---|---|---|---|
  | integral, 512 points | 2.94 s | 1.43 s | 1.07 s |
  | integral, 1024 points | 10.93 s | 7.80 s | 9.63 s |
  | integral, tabulated conical, 512 | 3.12 s | 2.66 s | 2.65 s |
  | rcwa lamellar in-plane, M=30 | 4.93 s | 0.057 s | 0.024 s |
  | rcwa lamellar in-plane, M=80 | 15.01 s | 0.733 s | 0.232 s |
  | rcwa lamellar conical, M=30 | 3.16 s | 0.046 s | 0.059 s |
  | rcwa lamellar conical, M=80 | 14.18 s | 0.761 s | 0.947 s |
  | rcwa sinusoid in-plane, M=30, 16 slabs | 64.00 s | 0.470 s | 0.322 s |
  | integral, 8 wavelengths, 512 | 13.10 s | 3.55 s | 5.49 s; 1.87 s at `workers=8` |
  | rcwa conical sinusoid, 8 wavelengths | 324.3 s | 4.99 s | 5.18 s; 0.66 s at `workers=8` |

  "After" is with the inherited environment, so it is what a caller now
  gets without setting anything. The conical and integral rows run the same
  arithmetic as the single-thread column before and differ from it by the
  machine's load, which moved between runs (12% to 25% at the start of each,
  and higher while they ran); read differences under about 30% there as
  noise. The in-plane RCWA gain over the single-thread column, 2 to 3x, is
  the decoupling.

- An RCWA backend, `rcwa`, written in-house: the first working slice of
  roadmap item 2. It is the Fourier modal method with Li's factorization
  rules (the inverse rule for the field component normal to the ridge walls),
  Rumpf's symmetric S-matrix rather than a transmittance matrix, and the full
  vectorial conical formulation, so off-plane incidence and TE-TM conversion
  are exact. `truncation` is the accuracy knob. A lamellar groove is one exact
  slab; any other single-valued profile is a staircase whose `slices` follow
  the truncation by default, and an undercut is refused. A `Problem` means
  what it means to the integral solver's tabulated mode, so efficiencies are
  absolute and `absorption` is the power that enters the structure. The core
  solves a layer stack with films above and below the grating region;
  `Problem` cannot describe one yet, so that is tested at the core only.
  Validated without the corpus: Fresnel on a flat interface to 1e-13, the
  Airy formula for a film to 1e-12, reciprocity falling with truncation, Li's
  rule converging on metallic TM where Laurent's is off by a fifth, and the
  integral solver agreeing order by order (to 5e-7 on a dielectric and 8e-4 in
  metallic TM, the gap shrinking as both refine). It registers, so it appears
  in the conformance suite and the Help menu; it has no GUI tab yet. The
  design note and what is not yet validated are in `docs/theory/rcwa.md`.

- Groove orientation in the AFM pipeline: a sample mounted a degree or two
  off square is measured and, when it matters, corrected, with no rotating by
  eye in Gwyddion. `metrology.core.orientation` finds the grooves band by band,
  links them into tracks and fits one slope to all of them, giving the angle
  and its standard error; the correction shifts each row so the grooves line
  up and takes the pixel along the groove normal, which is exact for a
  grating. On a synthetic 30° grating rotated 3° the blaze angle comes back
  as 29.998° against 29.966° uncorrected, the period 315.00 nm against 315.43,
  and the whole-image groove depth to within 0.2% where uncorrected it reads
  32% shallow. The `groove_orientation` setting is `auto` by default: correct
  at 0.2° or more when that is three standard errors from zero; `always`,
  `detect` and `off` are the others. The angle is reported whether or not it
  was applied, in the analysis result, the results panel, the Import and
  Boundary tabs, a boundary profile's metrics (so the `.json` export) and the
  metrics sidecar. A new wiki page, *Groove orientation*, explains it.

- `load_scan_y_size` reads the slow-axis scan size from `.spm`, `.ibw` and
  Gwyddion text files, scaled to the lines actually recorded, and returns
  `None` rather than assuming a square scan when the file does not record it.
  A scan without one is reported as not measured.
- AFM scan-line scar repair, the step Gwyddion calls "Remove scars": a run of
  one to four lines standing off the lines either side of it is found after
  image flattening and either interpolated down its columns or dropped from row
  grouping (`scar_repair = 'none' | 'interpolate' | 'drop'`, default `none`,
  with `scar_threshold`, `scar_min_length`, `scar_max_width`). It is its own
  stage in `metrology.core.scars`, not an image-flattening method, because it
  is used together with one. The number of lines repaired is shown in the
  Import tab and written to the results panel, the summary file and the
  boundary metrics file. On the synthetic scan with seven lines scarred over a
  third of their length the mean blaze angle reads 29.731° against 29.896°
  clean and the averaged groove is 1.01 nm off; interpolated, 29.895° and
  0.006 nm. Lines are judged levelled by their median step rather than by
  `align_rows`, and a scar's shadow on a rotated grating is set aside.
- A region mask: drag on the image in the AFM window's Import tab to exclude a
  rectangle, whole scan lines or whole columns. It is stored as settings
  (`mask_rect`, `mask_rows`, `mask_cols`, in pixels; `MASK_*` in `config.py`) so
  a batch or CLI run reproduces it, and it is recorded wherever scar repair is.
  Masked lines join no row group and masked columns yield no grooves. Masked,
  a 150 nm particle that moved the mean blaze angle by 0.49 degrees has no
  influence left: the result equals the clean scan's under the same mask.
- Wiki page "Scars and masking".

- Each window carries a **Go** menu: it switches to the other mode and returns
  to the chooser, with `Ctrl+Shift+M` for the latter. The chooser was a
  one-way door — picking a mode closed it with no route onward — which is half
  a fix for the problem it was added to solve.
- Switching reuses a window that is already open instead of constructing a
  second one. The metrology window holds a loaded scan and its settings, so a
  switch that produced an empty copy would read as the work having been lost.
  A window that was closed comes back with its state, since closing does not
  destroy it while the launcher holds a reference.

### Changed

- The default analysis now corrects a groove rotation of 0.2° or more. A
  synthetic scan with none, including the committed fixture, is unaffected
  (it measures -0.002 ± 0.002° and is left alone). **The real scans have not
  been re-measured with it**, since none are on the machine this was written
  on, so the real-scan numbers quoted in the wiki (the 0.49° flattening span,
  the ICC range) predate it. Set `GROOVE_ORIENTATION = 'detect'` to see what it
  would do to a scan without changing anything.
- The mask is applied before image flattening, so `align_rows` and `plane`
  now skip NaN pixels; a scan with nothing masked levels exactly as before.
- The Import tab's y axis counts scan lines from the top, where line 0 was
  always drawn, so it reads the same line numbers the mask stores. It used to
  count from the bottom.
- With both of the above, a scan is prepared in this order: region mask,
  image flattening, scar repair, groove orientation, tip correction. Scars
  are repaired before the orientation correction shifts rows, and the mask
  is drawn and stored in pixels of the scan as loaded, which is why the
  Import tab draws the image before that correction and the profile after.

### Fixed

- The integral solver read the profile parameter backwards and solved every
  asymmetric groove as its mirror image. It placed the curve at `x = t P`
  where `docs/conventions.md` puts `t` along `-d̂`, so a `Blazed` grating
  sent its blaze to the wrong side of the normal while the scalar solver, on
  the same problem, sent it where `geometry.blaze_direction` says. At
  `Blazed(9)`, α = 40°, 550 nm, p = 1667 nm the integral put 0.56 in m = −1
  and now puts 0.80 in m = +1, beside scalar's 0.89. Symmetric profiles,
  including every sinusoid, are unchanged, which is why the cross-method
  tests never saw it; an asymmetric case now guards it.
- `.ggp` files and the metrology pipeline's `x_norm` advance along `+d̂`,
  which the PCGrate TASTE match fixes, so reading or writing one now reflects
  the points into `t` (`FromProfileData.along_dispersion` and
  `dispersion_axis`). A measured groove therefore reaches both solvers facing
  the way PCGrate sees it. The integral solver's results on a `.ggp` are
  unchanged, to 2e-14, because its two reflections cancel; the scalar
  solver's change, because it had been solving every measured asymmetric
  groove mirrored. `read_profile(path).t` is now `1 − x` reversed rather than
  the file's `x`, and the boundary `.json` keeps storing `x_norm`.

- **The slow-axis scan size ignored `\Aspect Ratio`, so the measured groove
  rotation was too small by that ratio on most scans.** Nanoscope writes the
  same number twice in `Scan Size` even for a rectangular scan and records the
  shape in `\Aspect Ratio`; `io/spm.py::_slow_size_um` took the second value
  at face value. 1013 of the 1245 scans in this project's archive are not
  1:1 -- 911 at 4:1, 68 at 8:1, 26 at 2:1 -- and since a rotation is an
  arctangent of a ratio of the two axes, each of those reported a rotation too
  small by its aspect ratio. Some then declined an orientation correction they
  needed, because the understated angle fell below the threshold.

  **This changes measured numbers.** Any slow-axis length, and anything
  derived from one, differs for a non-1:1 scan: a 1.5 um scan at 4:1 now
  reports 0.375 um of slow axis instead of 1.5. Blaze angles taken from the
  fast axis are unaffected, as are all 1:1 scans -- which is why nothing
  caught this, the TASTE fixture and the synthetic scan both being 1:1.

  Checked against Gwyddion, which reads the field and states the result as
  `# Height:` in its export: of the 70 archive scans with an export beside
  them, 3 agreed before the fix and 70 agree after. The 3 were the 1:1 ones.

  One trap is pinned by its own test. Two different fields are called
  `\Aspect Ratio`: the scan list states the shape as a ratio with a colon,
  while every plane section carries a bare number of the same name that is not
  a shape at all (0.352542 on the TASTE file). Reading the section's value in
  preference -- which is how every other field in this reader works -- turns
  that file's 1.57 um slow axis into 4.44 um.

## 0.1.12 (2026-10-24)

### Added

- Asylum Research `.ibw` files open directly, in the AFM window and through
  `load_afm_data`, with no Gwyddion export step: `metrology.io.ibw` reads Igor
  binary waves (version 5, either byte order) without a new dependency, takes
  the ZSensor channel (Height when there is none) in the chosen scan direction,
  and refuses a channel that is not in metres. **Not yet checked against a file
  from an instrument**: the tests build files to the specification, and a
  synthetic 25° blaze at 400 nm comes back as 24.81° and 400.29 nm through the
  full analysis.

- `gratinglab.measured.MeasuredEfficiency`, measured efficiency with an
  uncertainty on every point: sparse per order, as parallel `order`,
  `wavelength`, `value` and `sigma` float arrays, with a
  `MeasurementProvenance` naming the facility, campaign, polarization and
  mount. It is a separate type rather than an `EfficiencyScan` with error
  fields, because measured data has none of a model's dense-grid shape, and
  the `uncertainties` package does not become a dependency. A polarization
  that was not recorded stays `None` and is never defaulted, and a slightly
  negative value from a weak order is kept as measured.
  `MeasuredEfficiency.from_energies` converts photon energies with the
  constant the materials code uses.
- `compare.align` and `compare.records` take measurements beside scans, so
  measured, scalar and integral results share one `Comparison`.
  `Comparison.measured_order(m)` gives the points with their error bars for
  plotting, and `Comparison.residuals(measurement, method)` the model at each
  measured point with the residual, the pull in units of sigma and a
  `summary()` that reports chi-squared with its point count. The model must
  have been swept on the measured wavelengths (`wavelength_grid()`, or its
  union with a denser grid): nothing interpolates, and an off-grid point is
  refused by name. `sweep` refuses a measurement and says to pass it to
  `align`, since resampling sparse data onto a grid would invent points.
  Existing comparisons are unchanged.
- The PCGrate XML results dialect, ported from `als_eff_meas`:
  `io.efficiency_table.read_xml_table`, `xml_to_scan` and `read_xml_scan`.
  TE, TM and unpolarized all read; TM never could in the reader it came
  from. The file is read by position, as the original did, and its scan
  axis is not read, so the wavelengths come from the caller and must pair one
  to one with the scan points. An order a point does not list is recorded as
  non-propagating there.
- `tests/test_measured.py` holds a worked comparison against the ported ALS
  campaign (OGRE laminar grating I, October 2021). It needs `als_eff_meas`
  and the private beamline tree, and skips without them; neither MC1 nor
  MAC1 holds that tree, so it has not yet run.

## 0.1.11 (2026-10-24)

### Added

- `corner_grading=p` on the integral solver: a corner-graded boundary mesh
  (Kress's sigmoidal transformation, M&P §4.6.6) for perfectly conducting
  gratings with corners. Cornered TM converged at first order on the
  equal-arc mesh; with `corner_grading=4` its energy-balance deviation is about
  1e-5 at 256 points on a lamellar, an ideal sawtooth and a 29.5°/70.5° blaze,
  against 2e-3 to 3e-2 before, and reciprocity improves by the same factor. It
  is opt-in, the default results are unchanged bit for bit, and a profile with
  no corner keeps the equal-arc mesh whatever is passed. It serves
  `conductivity="tabulated"` too, where `R + A` reaches 2e-6 to 3e-5 at 512
  points on a visible-light metal and on conical X-ray lamellar and blazed
  gratings, against 5e-3 to 0.16 on the equal-arc mesh. That needed two
  corrections the perfectly conducting solve does not use, a singularity
  subtraction of the double layer at corners and the tangential derivative
  taken on the assembled trace. The measurements are in
  `docs/theory/integral.md` §5 and `docs/findings.md`.
- `BoundaryCurve.polygon`: a graded mesh carries the vertices it lies on.
- `Profile.boundary(n, grading=p)`, which every built-in profile accepts, and
  `BoundaryCurve.graded`. A profile written against the older `boundary(n)`
  signature still works with the default mesh.
- The Setup tab acts on the material library instead of only describing it.
  **Add material…** picks a CXRO, refractiveindex.info or PCGrate table, asks
  its coating name, reads it to prove it is usable and copies it into your
  library; **Rescan library** re-reads the directories after a change made
  outside the window; **Open library folder** opens your library, creating it
  first. Both changes reach the geometry dock's Coating dropdown at once, which
  before needed a restart, and a selected coating that disappears falls back to
  "(none)" and says so. The page's links now open in the system browser, and
  its provenance section cites the kinds of table actually installed rather
  than Henke for every one.
- `materials.install(path, name=None, overwrite=False)`, the same install from
  Python. An unreadable table is refused before it reaches the library.
- Themes: **View → Theme** offers System (follows the operating system's light
  or dark setting, live), Light, Dark, Solarized Light, Solarized Dark and High
  contrast, remembered between runs and shared by both windows. Light is the
  look the windows had before, colour for colour. Every theme's text meets
  WCAG AA (4.5:1) and every plotted line 3:1, which `tests/test_gui_theme.py`
  asserts. The palette, every plot, the provenance panel, the Setup and theory
  pages and the geometry views all restyle in place.

- A **Corner mesh** selector on the Integral tab: off (the equal-arc mesh,
  still the default) or a grading order from 2 to 6, with a note on what it
  does. The form checks the graded mesh's own resolution floor, whose
  coarsest spacing is about twice the mean, and refuses an order that would
  put nodes closer than double precision separates, each on the field that
  caused it, using the solver's own check.

### Changed

- Diffraction orders are coloured by sign, cool for negative, warm for
  positive and neutral for m = 0, instead of from matplotlib's ten-colour
  cycle, which drew a sixteen-order scan with repeats (m = -3 and m = +7 were
  the same orange). An order keeps its colour when others are hidden, and the
  order list shows each line's colour beside its entry. The legend moved
  outside the axes, where it no longer covers the curves it names.
- The geometry dock's four groups share one label column, so every input
  starts at the same offset (they sat up to 64 px apart), and switching mount
  no longer moves them.
- The Grating Geometry tab: the 3D box fills more of its panel (a uniform
  zoom, so angles stay exact) with five ticks a side instead of nine; the
  cross-section's order labels are pushed outward along their own rays until
  none overlaps another or the title, where they used to stack at the top of
  the arc; the 3D cone-axis label moved off the tip where the orders
  cluster; and the Wavelength box is two compact rows with the scan's ends
  shown, giving the captions beside it the width.

## 0.1.10 (2026-10-23)

### Added

- `gratinglab-app`, a chooser that asks which window you want, and the single
  application bundle that launches it. Both halves have always been this one
  package, but the desktop layer shipped two console scripts and so two icons,
  the AFM one still called "AFM Blaze Meas" with the bundle id
  `org.jamccoy.afmblazemeas`, all of it inherited from the repository the
  metrology was merged from. Spotlight indexed it as `AFMBlazeMeas`, which
  sorts nowhere near GratingLab and names a project that no longer exists, so
  the half of the package in regular use was the half that could not be found
  — including by its author, looking for the AFM analysis inside the
  efficiency window. `tools/apps.py bundle --app gratinglab` now builds one
  `GratingLab.app` (`org.gratinglab.app`) that opens the chooser, and the
  chooser opens either window in the same process. The two direct scripts are
  unchanged.
- Bundles that skip the chooser are still available as `--app efficiency` and
  `--app metrology`, renamed "GratingLab Efficiency" and "GratingLab AFM" with
  `org.gratinglab.*` identities. Neither is built by default, so the Dock gets
  one icon.

- `flatten_profile` takes `return_report=True` and returns a
  `FlatteningReport` beside the flattened profile: the method and polynomial
  order actually used, how many features the background was fitted through,
  whether edge exclusion was skipped, and the RMS scatter of those features
  about the background they produced. The keyword follows
  `find_groove_positions`' `return_n_edge_rejected` and is off by default, so
  two-value unpacking keeps working and the GUI needed no change.
- The blaze-angle analysis asks for a report per row group and the per-sample
  summary records the aggregate. `analysis_data_*.csv` gains `Flatten_method`,
  `Flatten_fallback_groups`, `Flatten_residual_nm` and
  `Flatten_residual_groups`, appended so earlier columns keep their positions,
  and the text summary says whether every row group got the method asked for.
  A row group levelled by a fallback can now be found in a file afterwards
  instead of only in the console of the run that produced it.
- `Flatten_residual_nm` is the new diagnostic: how far the lands lie from the
  background fitted through them, in nanometres, and so how much the plane
  every blaze angle is quoted against deserves to be trusted. It is the worst
  row group's rather than the mean, since the weakest reference is the one that
  matters. It is absent for `linear` and `polynomial`, whose residual would
  describe the groove depth rather than the reference, and groups whose fit had
  no spare features are excluded and counted separately — three lands and a
  quadratic give zero by arithmetic, not by flatness.

### Changed

- `flatten_profile`'s fallbacks raise a `UserWarning` instead of printing to
  standard output. A print scrolls past in a batch run with nothing able to
  catch it; a warning can be turned into an error, counted, or routed to a log
  by whatever drives the pipeline. The suite runs with `filterwarnings =
  error`, so no test can now pass on a background fitted by a method it did not
  ask for — which is also how this change was checked: the full suite still
  passes, so nothing in it was flattening by fallback.

## 0.1.9 (2026-10-23)

### Changed

- The PyPI keywords no longer lead with x-ray work or name PCGrate. Nothing in
  the solvers or the metrology assumes a wavelength band — the grating equation
  carries none, and the materials layer reads CXRO x-ray tables and
  refractiveindex.info ultraviolet-to-infrared tables alike — so the old list
  told a reader the in-plane half of the package was not meant for them.
  `x-ray` is still there, among the bands rather than ahead of them, and
  `optics`, `diffraction grating`, `diffraction efficiency`, `integral method`,
  `groove profile`, `ultraviolet` and `infrared` are new. `pcgrate` is gone:
  it is a program this project validates against, not a feature of this one.
  `CITATION.cff` gains the same balance.
- Three classifiers added: `Operating System :: OS Independent`, which CI
  demonstrates on three platforms, and `Topic :: Scientific/Engineering ::
  Information Analysis` and `:: Visualization` beside Physics. There is no
  trove classifier for optics, and no Qt environment classifier is claimed,
  because a classifier cannot be conditional on an extra and a default install
  pulls no toolkit.
- The README names [endstation](https://pypi.org/project/endstation/) as what
  gratinglab is built on and links its PyPI page rather than only its
  repository, in the dependency paragraph and in the licence table. A reader on
  gratinglab's project page who wants the shared scaffolding wants the package,
  not the source tree.
- The README's licence section, which is the PyPI project page, was the single
  line "BSD-3-Clause." It now states the copyright, what the licence permits
  (commercial and closed-source use, no copyleft reach), what it requires
  (notice retention, and the no-endorsement clause), and that there is no
  warranty — which matters more than usual while the solvers are experimental.
  It adds how to cite a version, a table of the licences each extra brings in,
  the LGPL position of the `gui` extra and why PySide6 rather than PyQt6, and
  the two things the licence does not cover: optical constant tables, whose
  terms belong to whoever published them, and the measurement data, which is
  not ours to publish.

- The metrology window's Flattening page explains the levelling it documents
  rather than only tabulating it: why a background has to be removed at all
  (the substrate plane is the reference a blaze angle is quoted against, and
  an error in it is an error in every angle), why fitting through the lands is
  not the same as fitting through the whole trace, and where the proof that
  image flattening is free stops holding. It also documents the third
  levelling, the endpoint flattening applied on PCGrate export, which was
  described only in `docs/theory/metrology.md`; the page previously said there
  were two stages.
- The same page now states the `level_grooves` feature-detection criteria
  (half a period of separation, 20% of the height range in prominence), the
  shipped defaults including `FLATTEN_FEATURE = 'peaks'`, the fallbacks taken
  when too few features are found, and the fact that `level_grooves` needs a
  period estimate, so the period used to level a profile is not independent of
  the period the run reports. `tests/metrology/test_wiki.py` pins each of
  these against the code.

### Fixed

- `flatten_profile`'s `level_grooves` method no longer fits its background to
  whatever survived edge exclusion while reporting that it had used every
  feature. Edge exclusion runs after feature detection, and on a trace whose
  features sit near the ends it could leave one point, or none; the arrays had
  already been overwritten, so the fit used the remnant and the message saying
  otherwise was wrong. The exclusion is now skipped, with an accurate warning,
  whenever it would leave fewer features than the requested order needs.
- The same method no longer asks numpy for a polynomial the feature set cannot
  determine. The guard admitted two features while `FLATTEN_POLY_ORDER`
  defaults to 2, and a degree-k fit needs k+1 points, so the background came
  back rank-deficient behind a `RankWarning` that nothing was reading. The
  order is reduced to what the features support and the reduction is reported.
  Both paths produced a background that looked like a fit and was not one,
  which matters because flattening sets the plane every blaze angle on the
  scan is quoted against. `tests/metrology/test_flatten_profile.py` covers
  them; the four tests that pin the fixes fail against the previous code.

- The Flattening page's method table marked `level_grooves` as the default in
  a row that also listed its three feature choices, which read as though
  `both` were the default feature. It is `peaks`.

- The README and the docstrings name `tools/apps.py`, the script that builds
  the icons and the macOS bundles. They still named `tools/make_app.py` and
  `tools/make_icon.py`, which were merged into it, so every command the
  documentation gave for building a bundle or regenerating an icon failed. The
  README also now says that a bundle execs the console script in the
  checkout's virtualenv rather than freezing one, so it does not survive being
  copied to another machine, and names `tools/apps.py check`, which reports
  that breakage.

## 0.1.8 (2026-10-22)

### Added

- This changelog, linked from the README and the PyPI project page, and a
  versioning policy in `docs/roadmap.md` saying what 0.2 and 1.0 will mean.
  `tests/test_changelog.py` fails if the current version has no section.

### Changed

- The wiki's generated Home page matches the README's opening and no longer
  says its pages come from `docs/` alone.

### Fixed

- On Python 3.11 the `gui` extra no longer installs PySide6 6.12.0, which
  aborts the interpreter at exit (`Fatal Python error: bool_dealloc`) after
  heavy use of the windows. Versions 6.9 to 6.11 are unaffected, and on Python
  3.12 and later the bug is harmless, so the exclusion applies to that one
  release on 3.11 only.

## 0.1.7 (2026-10-22)

### Changed

- The README, every page under `docs/`, the metrology window's help pages and
  the docstrings are rewritten in plain prose. No code changed: each module's
  syntax tree with its docstrings blanked is identical to 0.1.6's.
- The README describes the AFM metrology as the part in regular use and marks
  the scalar and integral solvers as experimental, since they pass their test
  suite but have not yet been used for real analysis.
- Figures derived from PCGrate's output are no longer quoted in the public
  documentation. The comparisons against it are still described, as against an
  independent integral-method code on the private reference corpus.

### Fixed

- `docs/metrology/README.md` said the SEM of the row-group analysis scales as
  1/√N and that the row groups are treated as independent. Comparisons use the
  ICC-corrected SEM and the effective sample size, and the page now says so.
- The disagreement between the fitted facet angle and the depth-implied blaze
  angle of the TASTE scan was quoted as 3.6σ in three places, although
  `findings.md` explains that sigma is the wrong measure for it. It is now
  given as the 7.6° gap.
- Four docstrings were wrong: an unclosed code span in `gui/integral_options`,
  a claim in `gui/qt/scalar_tab` that no second solver tab exists, "both
  actions" for four buttons in the Boundary tab, and a cross-reference to a
  renamed section.

## 0.1.6 (2026-10-22)

### Added

- Asymmetric tip shapes. Each flank of the tip has its own half angle
  (`tip_half_angle_right_deg`; `None` keeps the symmetric cone), so a pyramidal
  probe's cross-section in the scan plane can be modelled.
  `gratinglab.metrology.core.tip.TIP_PRESETS` holds Bruker's ScanAsyst-Air
  scanned across the cantilever and along it in both orientations, and the
  metrology window has a Tip shape picker.
- The GitHub wiki also carries the metrology window's help pages, published as
  `AFM-<slug>`.

### Fixed

- Tip erosion now uses the reflected tip. The difference is invisible for a
  symmetric tip; with an asymmetric one the previous form gave a looser bound
  whose re-imaging did not reproduce the image.

## 0.1.5 (2026-10-22)

### Added

- Swath statistics (`gratinglab.metrology.boundary.swath_statistics`): the
  averaged groove is built separately for each horizontal band of scan rows
  and compared, reporting the shape spread between swaths, the groove-to-groove
  spread within one, per-swath metrics, the lateral shift down the scan, and how
  far the whole-scan groove lies from the swaths' mean. The Boundary tab has a
  "Swath statistics…" dialog with a per-swath CSV export.

## 0.1.4 (2026-10-22)

### Added

- A versioned boundary format, `gratinglab.io.boundary`: `.json` carrying the
  period, groove count, metrics and provenance, plus a `.csv` table. The
  Boundary tab exports `.ggp`, `.json` or `.csv`, and "Send to efficiency
  window" opens the groove there with its period.
- Tip correction controls in the metrology window.

### Fixed

- The Qt windows run Python's cycle collector on the GUI thread only. Garbage
  collection on a worker thread could free Qt objects there, which segfaulted on
  macOS and aborted the interpreter at exit on Python 3.11.

## 0.1.3 (2026-10-22)

### Fixed

- Groove averaging skips grooves that the scan edge cuts into
  (`average_grooves(..., require_full_period=True)`). A partial groove at the
  edge of a scan had been included in the average and distorted the exported
  shape.

## 0.1.2 (2026-10-22)

### Changed

- The version is written once, in `gratinglab.__version__`; `pyproject.toml`
  reads it from there.

## 0.1.1 (2026-10-21)

### Changed

- Requires endstation 0.2.2, which fixed its own encoding warnings.
- The README shows how to install with uv and links to PyPI, the
  documentation, the releases and the issue tracker.

### Fixed

- The critical angle is reported only where the refractive index is near 1,
  and is computed exactly there.

## 0.1.0 (2026-10-21)

First release on PyPI.

- AFM groove metrology (`gratinglab.metrology`): blaze angles per groove from
  Nanoscope `.spm` or Gwyddion text scans, row-group statistics with the
  correlation (ICC) correction, an averaged groove exported as a PCGrate `.ggp`
  boundary, and tip erosion with a certainty map.
- Scalar and integral-method efficiency solvers, the latter with perfect and
  finite conductivity, a comparison harness that also accepts imported
  efficiency tables, and a convergence harness.
- A materials layer reading CXRO exports, refractiveindex.info files and
  PCGrate `.ari` tables from user libraries; no optical constant tables are
  distributed.
- Resolving power from the measured groove count.
- Qt windows for the efficiency solvers and for the metrology.
