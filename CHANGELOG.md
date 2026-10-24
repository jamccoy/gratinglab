# Changelog

Notable changes, newest first. The project is pre-1.0, so a minor version may
change an interface; any such change is stated here when it is made. The
versioning policy, including what 0.2 and 1.0 will mean, is in
[`docs/roadmap.md`](docs/roadmap.md#versioning).

The format loosely follows [Keep a Changelog](https://keepachangelog.com). The
version is read from `gratinglab.__version__`, and a release commit moves the
entries under "Unreleased" into a section named for the new version;
`tests/test_changelog.py` fails if it does not.

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
