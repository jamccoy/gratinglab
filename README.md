# gratinglab

An open comparison platform for diffraction grating efficiency analysis.

Rigorous grating efficiency analysis today means either paying for PCGrate (integral
method) or assembling a patchwork of open RCWA codes that each carry their own
conventions, normalizations, and blind spots. Nothing open lets you run the *same*
physical problem through several methods and compare honestly.

`gratinglab` is built around one idea: **the physical problem and the numerical method
are separate, and the problem spec is serializable.** One `Problem` goes to every solver,
and the disagreement between methods is made visible rather than hidden.

```
Problem       (period + profile + materials — no solver fields, ever)
Illumination  (wavelength + direction + polarization)
                      |
                      v
Solver (plugin)  ->  Result (+ provenance)
```

## Status

Pre-alpha. Phase 0: problem spec, conventions, PCGrate interop, materials layer.

## Scope

Both mounts are first-class:

- **In-plane / classical** — UV–Vis–NIR optical engineering.
- **Extreme off-plane (conical), grazing incidence** — soft X-ray spectroscopy, where
  PCGrate is the community standard and general-purpose open RCWA codes are weakest.

Planned method backends, in order:

| Method | Status |
|---|---|
| Scalar (Harvey parametric, Kirchhoff, X-ray blaze form) | planned |
| RCWA (in-plane, conical, off-plane; Li factorization, S-matrix) | planned |
| C-method (Chandezon) | planned |
| Integral method | planned |

An already-computed PCGrate table is a first-class "method" in the comparison harness,
so scalar-vs-RCWA-vs-integral plots are possible without a live PCGrate license.

## Conventions

`docs/conventions.md` is **normative**. Read it before writing a solver or filing a bug
about a sign. The short version: time convention `exp(-iωt)`, lossy *n* = n′ + ik,
grating equation `sin α + sin β_m = mλ/(p sin γ)`, efficiencies **absolute**.

## Development

```bash
.venv/bin/python -m pytest
```

## License

BSD-3-Clause.
