"""The quick, cool pass over the suite: no Qt windows, no BLAS stampede.

Why this exists. A full `pytest` here is about 2750 tests and, on the macOS
CI runner, just under seven minutes; on an older laptop it is appreciably
longer and holds every core at full tilt for all of it. The ten slowest
tests account for only a quarter of that, so there is no hot spot to carve
out -- the cost is thousands of small solves, each of which asks LAPACK for
as many threads as the machine will give. That is the right default for a
real scan and the wrong one for a test run, where the work is already
plentiful and the threads only contend.

So this does three things a bare `pytest` does not:

* **Deselects the Qt tests** (`-m "not gui"`). Those build real windows,
  which is the slowest part of a run and the part a tired machine feels
  first.
* **Caps the BLAS pools at one thread** through the standard environment
  variables, which must be set before NumPy is imported -- hence a wrapper
  process rather than a fixture. The solver library already defaults to
  ``blas_threads=1`` per solve, but the tests reach plenty of NumPy that
  never goes through a solver.
* **Asks for a lower scheduling priority**, so a long run stays in the
  background of whatever else is open.

None of this changes what the tests assert, and none of it is a substitute
for the full suite: CI runs everything, on four platform/version legs, and
that is what a push is judged by. This is for the loop before the push.

Usage::

    python tools/test_fast.py                     # the quick pass
    python tools/test_fast.py tests/test_ggp.py   # narrower still
    python tools/test_fast.py --gui               # the Qt tests as well

Anything else on the command line is handed to pytest unchanged, so
``-k``, ``-x``, ``--durations`` and friends all work.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

#: Every library that might start a thread pool behind NumPy. Set together,
#: because which one is live depends on how NumPy and SciPy were built, and a
#: wheel from PyPI does not always use the one you would guess.
THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)

#: Unix only, and advisory: 10 is well below normal without being so low that
#: the run never finishes. Windows has no `os.nice`, so it simply goes without.
NICENESS = 10


def main(argv: list[str]) -> int:
    args = list(argv)
    include_gui = "--gui" in args
    if include_gui:
        args.remove("--gui")

    env = dict(os.environ)
    # setdefault, so `OMP_NUM_THREADS=4 python tools/test_fast.py` still means
    # what it says: this is a better default, not a decision taken from you.
    for var in THREAD_VARS:
        env.setdefault(var, "1")

    command = [sys.executable, "-m", "pytest"]
    if not include_gui:
        command += ["-m", "not gui"]
    command += args

    # Lowered here rather than in a `preexec_fn`: niceness is inherited, so
    # the pytest child starts where this process is, and `preexec_fn` is
    # documented as unsafe and does not exist on Windows anyway.
    try:
        os.nice(NICENESS)
    except (AttributeError, OSError):
        # No `nice` on this platform, or not permitted. Not worth failing a
        # test run over.
        pass

    print(
        "[test_fast] "
        + ("with" if include_gui else "without")
        + " the Qt tests, BLAS capped at "
        + env["OMP_NUM_THREADS"]
        + " thread(s)",
        file=sys.stderr,
    )
    return subprocess.call(command, cwd=REPO, env=env)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
