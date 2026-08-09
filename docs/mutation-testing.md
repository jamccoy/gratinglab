# Mutation testing

Coverage says a line ran. Mutation testing asks whether anything would have
**noticed if it were wrong** — it corrupts the source one edit at a time and
reports which corruptions the suite fails to catch. For a project whose entire
value is defensible numbers, that is the more useful question, and it has now
found real gaps twice (see [`findings.md`](findings.md#mutation-testing-found-a-real-gap)).

## Running it

```bash
.venv/bin/python -m pip install -e ".[mut]"
.venv/bin/mutmut run
```

On demand, never in CI. A full sweep is orders of magnitude slower than the
suite, and a surviving mutant is a prompt to think rather than a build failure.

Useful subsets and follow-ups:

```bash
.venv/bin/mutmut run "gratinglab.geometry.*"   # one module
.venv/bin/mutmut results                        # survivors, by name
.venv/bin/mutmut show gratinglab.geometry.x_beta__mutmut_13   # the diff
.venv/bin/mutmut browse                         # interactive
```

Note the mutant-name prefix is `gratinglab.`, **not** `src.gratinglab.` — the
src layout does not appear in the names.

## What is mutated, and what is not

`[tool.mutmut]` in `pyproject.toml` is the authority. In summary:

| Excluded | Why |
|---|---|
| `gui/` | Computes nothing. A surviving mutant there means a misplaced button, not a wrong number — which is the whole reason that layer was arranged to be thin. |
| `io/` | Pinned by round-trip and corpus-parsing tests, which a mutation operator does not reach usefully. |
| `illumination.py`, `problem.py`, `profiles.py` | **A tool limitation, not a judgement.** See below. |

Everything else — `geometry`, `checks`, `convergence`, `compare`, `result`,
`solvers/` — is mutated. That is where the physics is.

### The pydantic exclusion

mutmut rewrites every function into a trampoline plus a class attribute holding
the mutant table. On a `BaseModel` subclass, pydantic rejects the un-annotated
attribute outright and the module fails to import — so the choice is not
"mutate them badly", it is "mutate them or import them". Their validators are
covered by `tests/test_illumination.py` and `tests/test_profiles.py`.

If this matters more later, the fix is upstream (`model_config['ignored_types']`
would have to be set on every model, which is production code changed for a test
tool — not worth it today).

## Two things that make this work at all

Both are load-bearing and easy to break by accident.

**`pythonpath = ["src"]`** in `[tool.pytest.ini_options]`. mutmut copies the
project to `mutants/` and runs pytest there. Without this, `import gratinglab`
resolves to the *real* source via the editable install, every mutant "survives",
and the run reports a catastrophe that never happened. It is a no-op for an
ordinary `pytest`, where it names the same path the editable install already
provides.

**No `tests/__init__.py`.** With one, pytest imports conftest as
`tests.conftest`, which cannot resolve when pytest is driven in-process from a
console script — the current directory is not on `sys.path` there. The package
layout was not buying anything: the one thing it enabled, a relative
`from .conftest import reference_dir`, is better served by
[`tests/corpus.py`](../tests/corpus.py), since `conftest.py` is not an ordinary
module and importing from it is fragile by construction.

## Version pin

`mut = ["mutmut>=3.2,<3.3"]`, deliberately not in `dev`.

mutmut 3.3+ requires `libcst`, which ships no macOS x86_64 wheel and must be
compiled with a Rust toolchain — mutmut's own README names this case ("known for
at least the `x86_64-darwin` architecture"). 3.2.3 uses `parso` and is pure
Python. Lift the pin when either libcst ships the wheel or the machine grows a
Rust toolchain; nothing else depends on it.

## Where it stands

796 mutants across the physics core, 763 killed — **95.9%**. The per-module
breakdown, and the reason `check_reciprocity` accounts for 18 of the 33
survivors on its own, are in
[`findings.md`](findings.md#the-full-sweep-959-and-where-the-4-is).

A full sweep takes roughly half an hour on this machine (~6 mutants/second),
which is why it is a deliberate act rather than a CI step.

## Reading a survivor

Not every survivor is a missing test. An **equivalent mutant** changes the source
without changing behaviour, and no test can kill it. Two of the four survivors in
`geometry.beta` were equivalent, verified by evaluating both versions rather than
by assuming:

- `|s| <= 1.0` → `|s| <= 2.0`: `arcsin` of anything outside `[-1, 1]` is NaN
  regardless, so the guard is doing less than it looks.
- `np.asarray(s, dtype=np.float64)` → `np.asarray(s)`: NumPy promotes anyway.

So the workflow is: read the diff, work out whether the change is observable at
all, and only then decide whether it is a missing test. Record the equivalent
ones — in `findings.md` — so the next reader does not re-derive them.
