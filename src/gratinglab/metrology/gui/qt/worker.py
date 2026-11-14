"""
The analysis the window runs off the UI thread.

The thread is endstation's (`endstation.qt.jobs.JobRunner`, owned by the
window); this module is only the work. A single scan takes about 0.3 s, so the
thread gives little benefit a user would notice today, and matters once a
compare mode runs eight scans.

There is nothing to cancel inside the work: `analyze_single_file` reaches
NumPy-bound loops with no check point and reports no progress. When a newer
run supersedes this one, the window cancels it and the runner drops its
result on arrival, while the CPU finishes the computation.

The work returns plain data and never touches a figure; see the rule in
`gratinglab.metrology.gui`.
"""
from __future__ import annotations

import contextlib
import io

__all__ = ["analyse"]


def analyse(filename, settings):
    """Analyse one file, swallowing what the analysis prints.

    Returns ``(result, settings)``, so the window draws the result with the
    settings that produced it rather than whatever the form says by then.
    """
    from ...analyzer import analyze_single_file

    with contextlib.redirect_stdout(io.StringIO()):
        result = analyze_single_file(filename, show_plots=False, settings=settings)
    return result, settings
