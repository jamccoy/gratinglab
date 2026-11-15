"""The column a crowd of labels is set in, without a toolkit.

The 3D view at grazing incidence lands every order's tip in one knot, and
labels placed round each tip filled it. These hold the pure layout to what
it is for: crowds are found, a column is in tip order, and it keeps off
every tip and every ray it is given.
"""

import pytest

from gratinglab.gui import label_layout as L

BOUNDS = (0.0, 0.0, 600.0, 400.0)
SIZE = (30.0, 12.0)


class TestCrowds:
    def test_near_tips_group_and_a_far_one_does_not(self):
        tips = [(100, 100), (110, 104), (120, 100), (400, 300)]
        assert [sorted(c) for c in L.crowds(tips)] == [[0, 1, 2]]

    def test_a_chain_is_one_crowd(self):
        """Single linkage: each within reach of the next is one crowd, even
        when the ends are further apart than the reach."""
        tips = [(0, 0), (30, 0), (60, 0), (90, 0)]
        assert [sorted(c) for c in L.crowds(tips)] == [[0, 1, 2, 3]]

    def test_lone_tips_make_no_crowd(self):
        assert L.crowds([(0, 0), (200, 0), (0, 200)]) == []


class TestColumn:
    def test_rows_follow_the_tips_down_the_screen(self):
        tips = [(100, 110), (104, 100), (108, 120)]
        rects = L.column([0, 1, 2], tips, [SIZE] * 3, BOUNDS)
        tops = {i: r[1] for i, r in rects.items()}
        assert tops[1] < tops[0] < tops[2]

    def test_it_covers_no_tip_and_no_label_overlaps(self):
        tips = [(100, 110), (104, 100), (108, 120)]
        rects = L.column([0, 1, 2], tips, [SIZE] * 3, BOUNDS, keep_clear=tips)
        boxes = list(rects.values())
        assert not any(L.covers(b, t) for b in boxes for t in tips)
        assert not any(L.overlaps(a, b) for k, a in enumerate(boxes) for b in boxes[k + 1:])

    def test_it_keeps_off_a_ray_it_is_given(self):
        """A ray leaving the crowd to the right pushes the column off it."""
        tips = [(100, 200), (104, 204)]
        ray = ((100, 200), (500, 200))
        rects = L.column([0, 1], tips, [SIZE] * 2, BOUNDS, lines=[ray])
        assert not any(L.crosses(b, *ray) for b in rects.values())
        unconstrained = L.column([0, 1], tips, [SIZE] * 2, BOUNDS)
        assert any(L.crosses(b, *ray) for b in unconstrained.values()), (
            "the test proves nothing unless the ray was in the way to begin with")

    def test_it_goes_left_when_the_right_is_out_of_frame(self):
        tips = [(590, 200), (594, 206)]
        rects = L.column([0, 1], tips, [SIZE] * 2, BOUNDS)
        assert all(r[0] + r[2] <= 590 for r in rects.values())

    def test_no_room_is_none(self):
        tips = [(10, 10), (14, 12)]
        assert L.column([0, 1], tips, [(500.0, 300.0)] * 2, (0, 0, 520, 320)) is None


class TestCrosses:
    @pytest.mark.parametrize("a, b, expected", [
        ((0, 5), (20, 5), True),      # straight through
        ((0, 0), (20, 20), True),     # corner to corner
        ((0, 20), (20, 20), False),   # passes below
        ((5, 5), (6, 6), True),       # wholly inside
        ((-10, 5), (-1, 5), False),   # stops short
    ])
    def test_cases(self, a, b, expected):
        assert L.crosses((0, 0, 10, 10), a, b) is expected
