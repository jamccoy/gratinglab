# Mutation testing

Coverage reports whether a line ran. Mutation testing asks whether the test suite
would have detected the line being wrong: it corrupts the source one edit at a time
and reports which corruptions the suite fails to catch. For a project whose value
lies in producing defensible numbers, this is the more useful question, and it has so
far found real gaps twice (see [`findings.md`](findings.md#mutation-testing-found-a-real-gap)).

## Running it

```bash
.venv/bin/python -m pip install -e ".[mut]"
.venv/bin/mutmut run
```

Mutation testing is run on demand and never in CI. A full sweep is orders of
magnitude slower than the suite, and a surviving mutant calls for judgement instead of
constituting a build failure.

Useful subsets and follow-up commands:

```bash
.venv/bin/mutmut run "gratinglab.geometry.*"   # one module
.venv/bin/mutmut results                        # survivors, by name
.venv/bin/mutmut show gratinglab.geometry.x_beta__mutmut_13   # the diff
.venv/bin/mutmut browse                         # interactive
```

The mutant-name prefix is `gratinglab.` and not `src.gratinglab.`, because the src
layout does not appear in the names.

## Scope of mutation

`[tool.mutmut]` in `pyproject.toml` is the authority. In summary:

| Excluded | Why |
|---|---|
| `gui/` | Computes nothing. A surviving mutant there means a misplaced button and not a wrong number, which is the reason that layer was arranged to be thin. |
| `io/` | Pinned by round-trip and corpus-parsing tests, which a mutation operator does not reach usefully. |
| `illumination.py`, `problem.py`, `profiles.py` | Excluded because of a tool limitation, not as a judgement about these modules. See below. |

Everything else (`geometry`, `checks`, `convergence`, `compare`, `result`,
`solvers/`) is mutated, since that is where the physics resides.

### The pydantic exclusion

mutmut rewrites every function into a trampoline plus a class attribute holding the
mutant table. On a `BaseModel` subclass, pydantic rejects the un-annotated attribute
outright and the module fails to import. The available choice is therefore between
mutating these modules and importing them, and mutating them less thoroughly is not an
option. Their validators are covered by `tests/test_illumination.py` and
`tests/test_profiles.py`.

If this becomes more important later, the fix lies upstream:
`model_config['ignored_types']` would have to be set on every model, which would mean
changing production code for the sake of a test tool, and this is not justified at
present.

## Two configuration requirements

Both of the following are necessary for the sweep to work, and both are easy to break
by accident.

**`pythonpath = ["src"]`** in `[tool.pytest.ini_options]`. mutmut copies the project
to `mutants/` and runs pytest there. Without this setting, `import gratinglab`
resolves to the real source through the editable install, every mutant "survives", and
the run reports a failure of the suite that did not occur. The setting has no effect
on an ordinary `pytest` run, where it names the same path that the editable install
already provides.

**No `tests/__init__.py`.** If this file exists, pytest imports conftest as
`tests.conftest`, which cannot be resolved when pytest is driven in-process from a
console script, because the current directory is not on `sys.path` in that case. The
package layout provided no benefit. The one thing it enabled, a relative
`from .conftest import reference_dir`, is better served by
[`tests/corpus.py`](../tests/corpus.py), since `conftest.py` is not an ordinary module
and importing from it is inherently fragile.

## Version pin

`mut = ["mutmut>=3.2,<3.3"]`, which is deliberately kept out of `dev`.

mutmut 3.3+ requires `libcst`, which ships no macOS x86_64 wheel and must be compiled
with a Rust toolchain; mutmut's own README names this case ("known for at least the
`x86_64-darwin` architecture"). 3.2.3 uses `parso` and is pure Python. The pin can be
lifted when either libcst ships the wheel or the machine acquires a Rust toolchain;
nothing else depends on it.

## Current status

The physics core yields 796 mutants. The first full sweep killed 763 (95.9%), and 18
of the 33 survivors were in `check_reciprocity` alone. Those survivors have since been
closed, so `checks.py` and `geometry.py` both stand at 100% (the two remaining
survivors in `geometry` were verified to be equivalent mutants). The breakdown, and
the work that closing them required, are given in
[`findings.md`](findings.md#the-full-sweep-959-and-where-the-4-is).

The main lesson from that round is that a strategy which leaves no trace in its own
output cannot be tested through the public API. `check_reciprocity` selects which
orders to test, and its report does not state which orders were selected, so verifying
it required a recording stand-in solver instead of a more elaborate assertion against
the real one.

A full sweep takes roughly half an hour on this machine (~6 mutants/second), which is
why it is run deliberately and not as a CI step.

## Reading a survivor

A survivor does not always indicate a missing test. An *equivalent mutant* changes the
source without changing behaviour, and no test can kill it. Two of the four survivors
in `geometry.beta` were equivalent, which was verified by evaluating both versions
instead of being assumed:

- `|s| <= 1.0` → `|s| <= 2.0`: `arcsin` of any value outside `[-1, 1]` is NaN in
  either case, so the guard has less effect than it appears to.
- `np.asarray(s, dtype=np.float64)` → `np.asarray(s)`: NumPy promotes the type in
  either case.

The workflow is therefore to read the diff, determine whether the change is observable
at all, and only then decide whether it indicates a missing test. Equivalent mutants
are recorded in `findings.md`, so that later readers do not have to re-derive them.

### A third category: mutants observable only where the model is already invalid

M16-D left one survivor that is neither a missing test nor strictly an equivalent
mutant. In `_local_reflected_efficiency`, replacing $\sqrt{r_i}\sqrt{r_d}$ with
$\sqrt{r_i r_d}$ does change the output, but only where the branch wrap is non-uniform
across the groove, which is the Brewster region near normal incidence on a deep
groove. Everywhere else the two forms differ by a global sign that cancels out of
$|\int\cdots|^2$, and on Au at 1–6 nm they agree to $10^{-15}$.

Killing this mutant would require asserting on orders carrying $10^{-9}$, in a regime
that the solver already reports as outside what the model can represent. Such a test
would pin numerical noise only to raise the kill count. The mutant is therefore left
alive deliberately. The reasoning is recorded in the `_local_reflected_efficiency`
docstring and in `test_the_two_square_roots_are_taken_separately`, which asserts the
property that makes the mutation harmless instead of the arbitrary choice between the
two forms.

This case is treated as a separate category because the rule "survivor ⇒ missing
test" would have produced a poor test here, and the rule "survivor ⇒ equivalent
mutant" would have been a false claim.
