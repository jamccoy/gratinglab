"""Layout floors that cannot promise less than their contents need.

The implementation moved to :mod:`endstation.qt.sizing`, shared with the other
projects here -- two of which still hardcode a pixel count and will inherit the
derived version when they adopt it.

Re-exported rather than replaced at the call sites so `scalar_tab.py` and
`integral_tab.py` keep importing from where they always have, and so this
docstring stays next to the tabs it constrains.

Why it is derived rather than a literal: a ``QScrollArea`` with
``setWidgetResizable(True)`` will not shrink its widget below that widget's own
``minimumSizeHint``. A ``minimumWidth`` under that figure is not a tight budget,
it is a constraint Qt cannot satisfy, and it answers with a horizontal scrollbar
and a clipped control. That shipped twice here unnoticed, because nothing in the
suite asserted on a pixel. The number depends on the platform's push-button
metrics and scrollbar width, so a literal that is right on one machine is wrong
on another -- and CI runs three.
"""

from __future__ import annotations

from endstation.qt.sizing import fit_width_to_contents

__all__ = ["fit_width_to_contents"]
