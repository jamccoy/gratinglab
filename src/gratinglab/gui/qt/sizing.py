"""Layout floors that cannot promise less than their contents need.

The implementation moved to :mod:`endstation.qt.sizing`, where it is shared
with the other projects here. Two of those still hardcode a pixel count and
will inherit the derived version when they adopt it.

It is re-exported instead of being replaced at the call sites, so that
`scalar_tab.py` and `integral_tab.py` keep importing it from the same place
and this docstring stays next to the tabs it constrains.

The width is derived instead of given as a literal for the following reason.
A ``QScrollArea`` with ``setWidgetResizable(True)`` will not shrink its
widget below that widget's own ``minimumSizeHint``. A ``minimumWidth`` below
that figure is therefore a constraint Qt cannot satisfy, and Qt responds with
a horizontal scrollbar and a clipped control. This was shipped twice here
without being noticed, because nothing in the suite asserted on a pixel. The
figure depends on the platform's push-button metrics and scrollbar width, so
a literal that is correct on one machine is wrong on another, and CI runs on
three.
"""

from __future__ import annotations

from endstation.qt.sizing import fit_width_to_contents

__all__ = ["fit_width_to_contents"]
