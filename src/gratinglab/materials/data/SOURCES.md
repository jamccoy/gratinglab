# Vendored optical constants — provenance

Every file here is a **verbatim CXRO/Henke export**, not a re-derivation. They
are committed in their original format so that the reader shipped in
`optical.py` is the one exercised in anger, and so the header line carries the
density the tables were computed at.

| File | Material | Range | Points | Origin |
|---|---|---|---|---|
| `Au.txt` | Gold, density 19.32 g/cm³ | 200 eV – 2 keV (0.620 – 6.171 nm) | 498 | CXRO, <https://henke.lbl.gov/optical_constants/> |

## Citation

The underlying data is the Henke–Gullikson–Davis compilation. Work using these
constants should cite:

> B.L. Henke, E.M. Gullikson, and J.C. Davis, *X-ray interactions:
> photoabsorption, scattering, transmission, and reflection at E = 50–30000 eV,
> Z = 1–92*, Atomic Data and Nuclear Data Tables **54** (2), 181–342 (1993).

## Redistribution — read before the first public release

CXRO makes these tables freely available for scientific use and they are widely
redistributed, but **the terms have not been formally confirmed for
redistribution inside a BSD-3 package**, and this file is not a licence
determination. Before `gratinglab` is published:

1. Confirm CXRO's current terms for bundling the tables.
2. If bundling is not clearly permitted, delete this directory and fall back to
   `materials.from_cxro_file(path)` — the reader works on a user-supplied
   download, the tests that need a table skip the way the reference-corpus tests
   already do, and nothing else in the package changes.

The layer is deliberately built so that step 2 is a deletion rather than a
rewrite. Nothing in `gratinglab` requires a vendored file to exist.

## Adding a material

Download an `(Energy, Delta, Beta)` table from the CXRO page above, drop it in
here as `<Element>.txt`, and add a row to the table. `available()` reads the
directory, so no code changes.

Do **not** hand-edit these files. If a table needs trimming or resampling, do it
at load time — an edited copy stops being the thing this file claims it is.
