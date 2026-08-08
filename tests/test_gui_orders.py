"""Which orders the efficiency plot draws.

Headless. This replaces a single untested line in the widget layer
(`if values.max() < 1e-4: continue`), so the tests that matter most are the
ones showing the rule is now visible and overridable rather than silent.
"""

import numpy as np
import pytest

from gratinglab.gui.orders import (
    DEFAULT_LIMIT,
    DEFAULT_THRESHOLD,
    OrderSummary,
    carry_over,
    default_visible,
    describe,
    summarize,
)
from gratinglab.gui.scalar_options import ScalarOptionsState, build_options
from gratinglab.gui.state import FormState, build
from gratinglab.result import EfficiencyScan, Provenance
from gratinglab.solvers import scalar


def make_scan(peaks, propagating=None):
    """A scan whose orders peak at the given values, m = 0, 1, 2, ..."""
    peaks = np.asarray(peaks, dtype=float)
    efficiency = np.vstack([peaks * 0.5, peaks])  # two wavelengths; max is `peaks`
    if propagating is None:
        propagating = np.ones_like(efficiency, dtype=bool)
    return EfficiencyScan(
        wavelengths=np.array([1.0, 2.0]),
        orders=np.arange(len(peaks)),
        efficiency=np.where(propagating, efficiency, 0.0),
        propagating=np.asarray(propagating, dtype=bool),
        provenance=Provenance("test"),
    )


@pytest.fixture(scope="module")
def default_scan():
    """What the app shows on open: off-plane, 315.15 nm, λ 1-5 nm."""
    parsed = build(FormState())
    options = build_options(
        parsed.problem, parsed.illumination, parsed.wavelengths, ScalarOptionsState()
    )
    return scalar.solve(
        parsed.problem, parsed.illumination, parsed.wavelengths, **options
    )


class TestSummarize:
    def test_reports_one_entry_per_order_in_scan_order(self, default_scan):
        summaries = summarize(default_scan)
        assert [s.order for s in summaries] == [int(m) for m in default_scan.orders]

    def test_column_indexes_the_right_efficiency_column(self, default_scan):
        """The mapping a caller would otherwise have to rebuild, and could get
        wrong -- silently plotting m=+2's curve under m=+3's label."""
        for summary in summarize(default_scan):
            column = default_scan.efficiency[:, summary.column]
            assert summary.peak == pytest.approx(column.max())

    def test_peak_is_the_maximum_over_wavelength(self):
        assert [s.peak for s in summarize(make_scan([0.4, 0.1]))] == [0.4, 0.1]

    def test_an_order_evanescent_everywhere_is_marked(self):
        propagating = np.array([[True, False], [True, False]])
        summaries = summarize(make_scan([0.4, 0.3], propagating))
        assert summaries[0].ever_propagating
        assert not summaries[1].ever_propagating

    def test_an_order_propagating_only_sometimes_still_counts(self):
        """Orders open and close across a wavelength scan; one that carries
        power at any wavelength is worth drawing."""
        propagating = np.array([[True, True], [True, False]])
        assert summarize(make_scan([0.4, 0.3], propagating))[1].ever_propagating


class TestDefaultVisible:
    def test_keeps_what_the_old_cutoff_kept(self, default_scan):
        """The whole point of retaining 1e-4 as the default: the plot looks
        the same on first open as it did before the rule became visible."""
        summaries = summarize(default_scan)
        old_rule = {s.order for s in summaries if s.peak >= 1e-4}
        assert default_visible(summaries) == old_rule

    def test_the_default_case_hides_nothing(self, default_scan):
        """Measured, not assumed: all 16 orders clear the threshold here, so
        this milestone changes no pixels on the geometry the app opens with."""
        summaries = summarize(default_scan)
        assert len(summaries) == 16
        assert default_visible(summaries) == {s.order for s in summaries}

    def test_drops_orders_below_the_threshold(self):
        summaries = summarize(make_scan([0.5, 1e-6, 0.2]))
        assert default_visible(summaries) == {0, 2}

    def test_the_threshold_is_inclusive(self):
        """An order exactly at the bound is kept; a strict `>` would make the
        documented number off by one epsilon."""
        assert 0 in default_visible(summarize(make_scan([DEFAULT_THRESHOLD])))

    def test_caps_the_count_at_the_limit(self):
        """The real protection. An X-ray geometry can propagate hundreds of
        orders, and the old code's only defence was a wider legend."""
        summaries = summarize(make_scan(np.linspace(0.9, 0.1, 60)))
        assert len(default_visible(summaries)) == DEFAULT_LIMIT

    def test_the_cap_keeps_the_largest_peaks(self):
        summaries = summarize(make_scan([0.1, 0.9, 0.5, 0.7]))
        assert default_visible(summaries, limit=2) == {1, 3}

    def test_ties_resolve_by_order_so_redraws_are_stable(self):
        """A plot that redrew differently for the same scan would be its own
        small bug."""
        summaries = summarize(make_scan([0.5, 0.5, 0.5]))
        first = default_visible(summaries, limit=2)
        assert first == default_visible(summaries, limit=2) == {0, 1}

    def test_an_empty_scan_is_not_an_error(self):
        assert default_visible(()) == frozenset()

    def test_returns_orders_not_columns(self):
        """Orders survive a re-solve that changes how many propagate; column
        indices do not. carry_over depends on this."""
        summaries = (OrderSummary(order=-4, column=0, peak=0.5, ever_propagating=True),)
        assert default_visible(summaries) == {-4}


class TestCarryOver:
    def test_an_unchecked_order_stays_unchecked(self):
        """The behaviour that makes the panel feel right: nudge the blaze
        angle, re-solve, and your choices are still there.

        Order 1 is strong enough that the default rule would show it. It was
        unchecked, and it existed last time, so it stays hidden.
        """
        summaries = summarize(make_scan([0.5, 0.5, 0.5]))
        assert carry_over({0, 2}, {0, 1, 2}, summaries) == {0, 2}

    def test_a_checked_order_stays_checked_even_below_threshold(self):
        """An explicit choice outranks the default rule -- otherwise looking
        at a weak order would be impossible across a re-solve."""
        summaries = summarize(make_scan([0.5, 1e-9]))
        assert 1 in carry_over({0, 1}, {0, 1}, summaries)

    def test_a_newly_appeared_order_gets_the_default_rule(self):
        summaries = summarize(make_scan([0.5, 0.4, 1e-9]))
        result = carry_over({0}, {0}, summaries)
        assert 1 in result  # new and strong: shown
        assert 2 not in result  # new and negligible: hidden

    def test_hidden_on_purpose_is_distinguished_from_never_seen(self):
        """The whole reason carry_over takes the previous *orders* and not
        just the previous *selection*. Both orders are absent from the visible
        set; only one of them is new, and only that one gets the default."""
        summaries = summarize(make_scan([0.5, 0.5, 0.5]))
        result = carry_over({0}, {0, 1}, summaries)
        assert 1 not in result  # existed, unchecked -> stays hidden
        assert 2 in result  # brand new, strong -> shown

    def test_an_order_that_no_longer_exists_is_forgotten(self):
        """A remembered selection would resurrect on an unrelated re-solve."""
        assert 99 not in carry_over({0, 99}, {0, 99}, summarize(make_scan([0.5])))

    def test_deselecting_everything_survives_a_re_solve(self):
        """A legitimate state, and it must stick, or the None button would
        undo itself the moment anything was re-solved."""
        summaries = summarize(make_scan([0.5, 0.4]))
        assert carry_over(set(), {0, 1}, summaries) == frozenset()

    def test_first_solve_is_just_the_default(self):
        summaries = summarize(make_scan([0.5, 1e-9]))
        assert carry_over(frozenset(), frozenset(), summaries) == default_visible(
            summaries
        )


class TestDescribe:
    def test_signs_the_order(self):
        """`+3` and `-3` diffract to opposite sides; an unsigned `3` next to
        `-3` in a list reads as ambiguous."""
        summary = OrderSummary(order=3, column=0, peak=0.4127, ever_propagating=True)
        assert describe(summary).startswith("m=+3")

    def test_shows_the_peak_that_the_rule_judged(self):
        summary = OrderSummary(order=3, column=0, peak=0.4127, ever_propagating=True)
        assert "0.4127" in describe(summary)

    def test_an_evanescent_order_says_so_instead_of_showing_zero(self):
        """`peak 0.0000` would read as a numerical result rather than as
        'this order does not propagate here'."""
        summary = OrderSummary(order=9, column=0, peak=0.0, ever_propagating=False)
        assert "evanescent" in describe(summary)
        assert "peak" not in describe(summary)


class TestPurity:
    def test_imports_no_toolkit_and_no_plotting(self):
        import ast
        import inspect

        import gratinglab.gui.orders as module

        imported: set[str] = set()
        for node in ast.walk(ast.parse(inspect.getsource(module))):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])

        assert not {"tkinter", "PySide6", "PyQt6", "matplotlib"} & imported
