r"""What a grating geometry looks like.

Headless -- no Qt, no matplotlib, no display. The whole point of the split is
that "does the beam strike the active facet" and "is an evanescent order still
in the picture" are checkable without opening a window.

The reference geometry throughout is the project's primary application: a
315.15 nm blazed grating at :math:`\delta = 29.5°`, off-plane at
:math:`\alpha = 25°`, :math:`\gamma = 1.5°`, where
:math:`\beta_b = 2\delta - \alpha = +34.0°`.
"""

import numpy as np
import pytest

from gratinglab.geometry import (
    blaze_direction,
    is_propagating,
    order_range,
    sin_beta,
)
from gratinglab.gui import diagram
from gratinglab.gui.diagram import (
    MAX_ORDER_LABELS,
    Diagram,
    blaze_jump,
    blaze_targets,
    build,
    facet_normal,
    label_orders,
    nearest_index,
    order_marks,
    order_span,
    x_nm,
)
from gratinglab.gui.provenance import ALARM_TAGS
from gratinglab.illumination import Illumination
from gratinglab.problem import Problem
from gratinglab.profiles import Blazed, Lamellar, Sinusoidal

UNPOL = "unpolarized"

PERIOD = 315.15
DELTA = 29.5
ALPHA = 25.0
GAMMA = 1.5


def blazed(period=PERIOD, delta=DELTA, antiblaze=70.5):
    return Problem(
        period=period, profile=Blazed(blaze_angle=delta, antiblaze_angle=antiblaze)
    )


def offplane(alpha=ALPHA, gamma=GAMMA):
    return Illumination.offplane(graze=gamma, azimuth=alpha, polarization=UNPOL)


@pytest.fixture(scope="module")
def reference():
    return build(blazed(), offplane(), 3.0)


def arrow(d: Diagram, tag: str, panel: str = "main"):
    """The single arrow with this tag on this panel."""
    found = [a for a in d.arrows if a.tag == tag and a.panel == panel]
    assert len(found) == 1, f"expected one {tag!r} arrow, found {len(found)}"
    return found[0]


def unit(a):
    v = np.array([a.x1 - a.x0, a.y1 - a.y0])
    return v / np.linalg.norm(v)


class TestOrderSpan:
    def test_the_propagating_subset_is_exactly_order_range(self):
        """The picture must agree with the solver about what propagates."""
        marks = order_marks(blazed(), offplane(), 3.0)
        live = {m.order for m in marks if m.propagating}
        expected = set(
            int(m)
            for m in order_range(3.0, PERIOD, np.sin(np.radians(ALPHA)), np.sin(np.radians(GAMMA)))
        )
        assert live == expected

    def test_the_span_pads_beyond_the_propagating_set(self):
        """An order that has just passed off should leave a trace; order_range
        alone would erase it."""
        live = order_range(3.0, PERIOD, np.sin(np.radians(ALPHA)), np.sin(np.radians(GAMMA)))
        span = order_span(3.0, PERIOD, np.sin(np.radians(ALPHA)), np.sin(np.radians(GAMMA)), pad=2)
        assert span.min() == live.min() - 2
        assert span.max() == live.max() + 2

    def test_every_sin_beta_matches_the_geometry_module(self):
        """This module arranges; it never re-derives."""
        marks = order_marks(blazed(), offplane(), 3.0)
        for mark in marks:
            expected = float(
                sin_beta(
                    mark.order, 3.0, PERIOD,
                    np.sin(np.radians(ALPHA)), np.sin(np.radians(GAMMA)),
                )
            )
            assert mark.sin_beta == pytest.approx(expected)

    def test_beta_is_none_exactly_when_the_order_is_evanescent(self):
        """Not nan, not 0.0 -- an evanescent order has no direction, and a
        number there would invite one."""
        marks = order_marks(blazed(), offplane(), 3.0)
        assert any(m.propagating for m in marks) and any(not m.propagating for m in marks)
        for mark in marks:
            assert (mark.beta is None) == (not mark.propagating)
            if mark.beta is not None:
                assert np.isfinite(mark.beta)


class TestEvanescentAreNeverDropped:
    """conventions.md §4: retained in the record, never silently dropped,
    never NaN. A picture is a record."""

    #: Order 3 passes off above 3.912 nm at this geometry, so it propagates at
    #: LIVE and does not at PAST. Pinned by test_the_case_really_is_a_cutoff --
    #: which caught the first attempt at these numbers, where the order was
    #: already evanescent at both wavelengths and the whole class was vacuous.
    M = 3
    LIVE = 3.0
    PAST = 4.5

    def test_the_case_really_is_a_cutoff(self):
        """Non-vacuity. Without this, the tests below pass on a geometry where
        nothing ever passes off, and check nothing at all."""
        sa, sg = np.sin(np.radians(ALPHA)), np.sin(np.radians(GAMMA))
        assert is_propagating(sin_beta(self.M, self.LIVE, PERIOD, sa, sg))
        assert not is_propagating(sin_beta(self.M, self.PAST, PERIOD, sa, sg))

    def test_an_order_that_passes_off_stays_in_the_diagram(self):
        before = build(blazed(), offplane(), self.LIVE)
        after = build(blazed(), offplane(), self.PAST)

        assert self.M in {m.order for m in before.orders}
        assert self.M in {m.order for m in after.orders}, "silently dropped"

        assert before.mark(self.M).propagating
        assert not after.mark(self.M).propagating
        assert after.mark(self.M).beta is None
        assert abs(after.mark(self.M).sin_beta) > 1.0, "clipped, not unclipped"

    def test_an_evanescent_order_gets_no_ray_and_is_not_faked_into_one(self):
        after = build(blazed(), offplane(), self.PAST)
        assert not any(a.order == self.M for a in after.arrows)
        assert any(
            m.order == self.M and not m.propagating for m in after.orders
        ), "must still appear on the ladder"

    def test_the_evanescent_order_is_drawn_outside_the_propagating_window(self):
        """Present *and* visibly outside, not present and indistinguishable."""
        after = build(blazed(), offplane(), self.PAST)
        ladder = [m for m in after.markers if m.panel == "ladder"]
        evanescent = [m for m in ladder if m.tag == "evanescent"]
        assert evanescent, "no evanescent marker on the ladder"
        assert all(abs(m.x) > 1.0 for m in evanescent)

    def test_the_caption_counts_both_kinds(self):
        after = build(blazed(), offplane(), self.PAST)
        text = " ".join(c.text for c in after.captions)
        assert "evanescent" in text
        dark = sum(1 for m in after.orders if not m.propagating)
        assert str(dark) in text


class TestHandedness:
    """Ties `Blazed`'s drawn orientation to `geometry.blaze_direction`.

    Nothing before this module drew a profile and a ray in the same frame, so
    nothing before it could catch a mirror flip between them -- and there was
    one to catch. See docs/conventions.md §3 and docs/findings.md.
    """

    def test_the_active_facet_reflects_the_incident_ray_into_the_blaze_direction(
        self, reference
    ):
        incident = unit(arrow(reference, "incident"))
        normal = unit(arrow(reference, "normal"))
        reflected = incident - 2.0 * (incident @ normal) * normal
        assert np.allclose(reflected, unit(arrow(reference, "blaze")))

    def test_the_incident_ray_strikes_the_active_facet_from_the_front(self, reference):
        """A mirrored drawing puts the beam on the anti-blaze facet, arriving
        from behind the one it is supposed to illuminate."""
        incident = unit(arrow(reference, "incident"))
        normal = unit(arrow(reference, "normal"))
        assert incident @ normal < 0.0

    def test_the_facet_normal_sits_at_the_blaze_angle(self):
        nx, ny = facet_normal(blazed())
        assert np.degrees(np.arctan2(nx, ny)) == pytest.approx(DELTA)

    def test_the_blaze_arrow_azimuth_equals_blaze_direction(self, reference):
        expected = blaze_direction(np.radians(DELTA), np.radians(ALPHA))
        assert arrow(reference, "blaze").azimuth == pytest.approx(expected)
        assert np.degrees(expected) == pytest.approx(34.0)

    def test_at_littrow_the_blaze_arrow_retroreflects(self):
        """alpha = delta sends the blazed order straight back along the
        incident ray -- an independent check of the whole angle chain."""
        d = build(blazed(delta=20.0), offplane(alpha=20.0), 3.0)
        assert np.allclose(unit(arrow(d, "blaze")), -unit(arrow(d, "incident")))

    def test_x_nm_runs_against_the_dispersion_direction(self):
        """The convention itself, as a one-line assertion."""
        assert x_nm(1.0, PERIOD) == pytest.approx(-PERIOD)


class TestKnownGeometry:
    @pytest.mark.parametrize("alpha", [0.0, 10.0, 25.0, -30.0])
    @pytest.mark.parametrize("gamma", [90.0, 30.0, 1.5])
    def test_the_zeroth_order_is_the_specular_reflection(self, alpha, gamma):
        r""":math:`\beta_0 = -\alpha`, so the m=0 arrow mirrors the incident
        one about :math:`\hat{n}`."""
        d = build(blazed(), offplane(alpha=alpha, gamma=gamma), 3.0)
        zero = arrow(d, "zero")
        assert np.degrees(zero.azimuth) == pytest.approx(-alpha, abs=1e-9)

    def test_every_order_tip_lies_on_the_cone(self, reference):
        r"""":math:`k_z` is conserved", drawn. All propagating orders share the
        same projected length, so their tips are equidistant from the strike
        point -- which is what makes the cone circle meaningful rather than
        decorative."""
        cone = next(p for p in reference.paths if p.tag == "cone")
        centre = np.array([arrow(reference, "incident").x1, arrow(reference, "incident").y1])
        radius = np.hypot(cone.x - centre[0], cone.y - centre[1])
        assert np.allclose(radius, radius[0])

        for a in reference.arrows:
            if a.order is None:
                continue
            assert np.hypot(a.x1 - centre[0], a.y1 - centre[1]) == pytest.approx(radius[0])

    def test_orders_are_symmetric_about_the_normal_at_normal_incidence(self):
        d = build(blazed(), offplane(alpha=0.0, gamma=90.0), 60.0)
        by_order = {a.order: a for a in d.arrows if a.order is not None}
        for m in (1, 2, 3):
            if m in by_order and -m in by_order:
                assert by_order[m].azimuth == pytest.approx(-by_order[-m].azimuth)

    def test_the_surface_spans_the_requested_number_of_periods(self):
        d = build(blazed(), offplane(), 3.0, periods=3)
        surface = next(p for p in d.paths if p.tag == "surface")
        assert surface.x.max() - surface.x.min() == pytest.approx(3 * PERIOD)


class TestBlazeArrowOnlyWhenThereIsOne:
    @pytest.mark.parametrize(
        "profile",
        [Sinusoidal(depth_fraction=0.15), Lamellar(depth_fraction=0.2, duty_cycle=0.5)],
        ids=["sinusoid", "lamellar"],
    )
    def test_no_blaze_arrow_for_a_profile_without_a_blaze_angle(self, profile):
        d = build(Problem(period=PERIOD, profile=profile), offplane(), 3.0)
        assert not any(a.tag == "blaze" for a in d.arrows)

    def test_but_a_blazed_profile_does_get_one(self, reference):
        """Non-vacuity for the test above."""
        assert any(a.tag == "blaze" for a in reference.arrows)

    def test_and_the_caption_says_why_it_is_missing(self):
        d = build(Problem(period=PERIOD, profile=Sinusoidal(depth_fraction=0.15)), offplane(), 3.0)
        text = " ".join(c.text for c in d.captions)
        assert "no blaze angle" in text


class TestCaptionsDoNotLie:
    def test_a_correct_diagram_is_not_styled_as_a_problem(self, reference):
        """Nothing here is wrong, so nothing here may look wrong -- the same
        rule gui/provenance.py exists to enforce."""
        assert all(c.tag not in ALARM_TAGS for c in reference.captions)

    def test_a_caption_states_gamma_and_that_k_z_is_out_of_the_plane(self, reference):
        text = " ".join(c.text for c in reference.captions)
        assert "k_z" in text
        assert "1.5" in text

    def test_a_caption_says_zeta_is_not_drawn(self, reference):
        """zeta is a 3D angle to the facet plane, not an angle in this view --
        1.495 deg against 85.5 deg. An arc labelled zeta would be the worst
        lie available here."""
        text = " ".join(c.text for c in reference.captions)
        assert "ζ" in text and "not drawn" in text

    def test_a_caption_says_ray_lengths_are_arbitrary(self, reference):
        assert "arbitrary" in " ".join(c.text for c in reference.captions)

    def test_in_plane_says_the_diagram_is_the_whole_geometry(self):
        """At gamma = 90 the projection *is* the truth, and hedging would be
        its own small dishonesty."""
        d = build(blazed(), Illumination.classical(alpha=10.0, polarization=UNPOL), 600.0)
        text = " ".join(c.text for c in d.captions)
        assert "in-plane" in text
        assert "not a projection" in text
        assert "out of the page" not in text


class TestLabelThinning:
    def test_labels_are_capped(self):
        d = build(blazed(), offplane(), 1.0)  # short wavelength: many orders
        live = sum(1 for m in d.orders if m.propagating)
        assert live > MAX_ORDER_LABELS, "vacuous unless the case is crowded"
        assert sum(1 for a in d.arrows if a.order is not None and a.label) <= MAX_ORDER_LABELS

    def test_every_propagating_order_still_gets_an_arrow(self):
        """Only the *label* is dropped. Dropping the arrow would be the silent
        omission this module exists to avoid."""
        d = build(blazed(), offplane(), 1.0)
        drawn = {a.order for a in d.arrows if a.order is not None}
        assert drawn == {m.order for m in d.orders if m.propagating}

    def test_the_zeroth_order_is_always_labelled(self):
        assert 0 in label_orders(order_marks(blazed(), offplane(), 1.0))

    def test_the_caption_reports_the_true_total_when_labels_are_thinned(self):
        d = build(blazed(), offplane(), 1.0)
        assert "labelled" in " ".join(c.text for c in d.captions)


class TestWavelengthChoice:
    def test_nearest_index_snaps_to_the_grid(self):
        grid = np.linspace(1.0, 5.0, 200)
        assert nearest_index(grid, 3.0) == pytest.approx(np.argmin(np.abs(grid - 3.0)))

    def test_blaze_jump_finds_an_order_inside_the_default_scan(self):
        """At the shipped defaults lambda_b(m=2) = 4.05 nm, inside 1-5 nm."""
        grid = np.linspace(1.0, 5.0, 200)
        index, reason = blaze_jump(grid, blazed(), offplane())
        assert index is not None
        assert reason == "", "a usable jump must carry no excuse"
        assert 4.0 < grid[index] < 4.1

    def test_a_profile_without_a_blaze_angle_says_so(self):
        grid = np.linspace(1.0, 5.0, 200)
        index, reason = blaze_jump(
            grid, Problem(period=PERIOD, profile=Sinusoidal(depth_fraction=0.15)), offplane()
        )
        assert index is None
        assert "no blaze angle" in reason

    def test_a_scan_that_misses_every_blaze_wavelength_says_so(self):
        """Both branches covered, so neither assertion is vacuous."""
        grid = np.linspace(10.0, 12.0, 50)
        index, reason = blaze_jump(grid, blazed(), offplane())
        assert index is None
        assert "10.00-12.00 nm scan" in reason

    def test_the_reason_is_non_empty_exactly_when_there_is_no_index(self):
        grid = np.linspace(1.0, 5.0, 200)
        for problem in (blazed(), Problem(period=PERIOD, profile=Sinusoidal(depth_fraction=0.15))):
            index, reason = blaze_jump(grid, problem, offplane())
            assert (index is None) == bool(reason)

    def test_blaze_targets_are_empty_for_an_unblazed_profile(self):
        assert blaze_targets(Problem(period=PERIOD, profile=Sinusoidal(depth_fraction=0.1)), offplane()) == ()


class TestDiagramShape:
    def test_every_panel_has_limits(self, reference):
        for panel in ("main", "ladder"):
            assert panel in reference.limits
            (x0, x1), (y0, y1) = reference.limits[panel]
            assert x1 > x0 and y1 > y0

    def test_on_partitions_the_primitives_by_panel(self, reference):
        total = len(reference.paths) + len(reference.arrows) + len(reference.markers)
        assert sum(len(reference.on(p)) for p in ("main", "ladder")) == total

    def test_every_tag_used_has_a_colour(self, reference):
        used = {p.tag for p in reference.paths} | {a.tag for a in reference.arrows}
        used |= {m.tag for m in reference.markers}
        assert used <= set(diagram.TAG_COLORS)

    def test_mark_raises_for_an_order_that_is_not_there(self, reference):
        """An order that should be in the record and isn't is a bug, not a
        default."""
        with pytest.raises(KeyError):
            reference.mark(999)

    def test_there_is_no_gamma_panel_any_more(self, reference):
        """Retired in M13-I. `diagram3d` draws the real cone, with γ as an
        angle in a real scene rather than a sliver standing in beside one, and
        two drawings of one angle are two answers to one question.

        γ did not lose its evidence with it: see
        `test_gui_diagram3d.py::TestTheConeProperty`, which asserts every drawn
        ray -- including the incident one -- sits at exactly γ from the cone
        axis, and that the transverse extent is `sin γ` rather than a
        normalised stand-in.
        """
        assert "cone" not in reference.limits
        assert not [p for p in reference.paths if p.panel == "cone"]
        assert not [a for a in reference.arrows if a.panel == "cone"]

    def test_but_the_cone_itself_is_still_drawn_in_the_main_panel(self, reference):
        """Non-vacuity, and the distinction the removal turns on: `"cone"` as
        a *panel* is gone; `"cone"` as a *tag* -- the end-on circle every
        propagating order lands on -- is what made the panel redundant, and
        stays."""
        circle = next(p for p in reference.paths if p.tag == "cone")
        assert circle.panel == "main"

    def test_every_ray_carries_the_same_out_of_page_marker(self):
        """What 'k_z is conserved' looks like: the same glyph on every ray."""
        d = build(blazed(), offplane(), 3.0)
        tips = [m for m in d.markers if m.glyph == "out-of-page"]
        assert len(tips) == sum(1 for m in d.orders if m.propagating)


class TestPurity:
    def test_imports_no_toolkit_and_no_plotting(self):
        """A diagram module is the single most tempting place to reach for
        'just a QColor'."""
        import ast
        import inspect

        imported: set[str] = set()
        for node in ast.walk(ast.parse(inspect.getsource(diagram))):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module.split(".")[0])

        assert not {"tkinter", "PySide6", "PyQt6", "matplotlib"} & imported
