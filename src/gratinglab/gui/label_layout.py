"""Labels whose anchors crowd together, set in a column beside the crowd.

At grazing incidence every ray of the 3D view leaves along nearly the same
line, so their tips land within a few pixels of one another. Each label
placed beside its own tip then fills that knot, where the rays and tip
markers cross it. This module finds such crowds and lays their labels out as
one column beside the crowd, in the order the tips run down the screen, so
the leader lines joining each label to its tip do not cross.

Pure: screen points and ``(x, y, w, h)`` rectangles in, rectangles out, and
no toolkit. :class:`gratinglab.gui.qt.scene3d_panel.SceneCanvas` calls it and
paints the result.
"""

from __future__ import annotations

from collections.abc import Sequence

__all__ = ["CROWD_PX", "GAP_PX", "column", "covers", "crosses", "crowds"]

#: A tip nearer than this to another tip joins that tip's crowd.
CROWD_PX = 36.0
#: Between the crowd's outermost tip and the column.
GAP_PX = 14.0
#: Between rows of the column.
ROW_GAP = 2.0
#: How far a label keeps from any point it must not cover.
CLEAR_PX = 3.0

Point = tuple[float, float]
Rect = tuple[float, float, float, float]


def crowds(tips: Sequence[Point], reach: float = CROWD_PX) -> list[list[int]]:
    """Groups of two or more tips, each tip within ``reach`` of another in its group."""
    n = len(tips)
    group = list(range(n))

    def root(i: int) -> int:
        while group[i] != i:
            group[i] = group[group[i]]
            i = group[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            dx, dy = tips[i][0] - tips[j][0], tips[i][1] - tips[j][1]
            if dx * dx + dy * dy <= reach * reach:
                group[root(i)] = root(j)
    found: dict[int, list[int]] = {}
    for i in range(n):
        found.setdefault(root(i), []).append(i)
    return [members for members in found.values() if len(members) > 1]


def overlaps(a: Rect, b: Rect) -> bool:
    return a[0] < b[0] + b[2] and b[0] < a[0] + a[2] and a[1] < b[1] + b[3] and b[1] < a[1] + a[3]


def covers(rect: Rect, point: Point, pad: float = CLEAR_PX) -> bool:
    x, y, w, h = rect
    return x - pad <= point[0] <= x + w + pad and y - pad <= point[1] <= y + h + pad


def inside(rect: Rect, bounds: Rect) -> bool:
    x, y, w, h = rect
    bx, by, bw, bh = bounds
    return x >= bx and y >= by and x + w <= bx + bw and y + h <= by + bh


def crosses(rect: Rect, a: Point, b: Point) -> bool:
    """Whether the segment from ``a`` to ``b`` passes through ``rect``.

    Liang and Barsky's clip: the segment is cut to each of the rectangle's
    four slabs in turn, and crosses it if anything of it is left."""
    x, y, w, h = rect
    lo, hi = 0.0, 1.0
    dx, dy = b[0] - a[0], b[1] - a[1]
    for step, gap in ((-dx, a[0] - x), (dx, x + w - a[0]), (-dy, a[1] - y), (dy, y + h - a[1])):
        if step == 0:
            if gap < 0:
                return False
            continue
        t = gap / step
        if step < 0:
            lo = max(lo, t)
        else:
            hi = min(hi, t)
        if lo > hi:
            return False
    return True


def column(
    members: Sequence[int],
    tips: Sequence[Point],
    sizes: Sequence[tuple[float, float]],
    bounds: Rect,
    keep_clear: Sequence[Point] = (),
    taken: Sequence[Rect] = (),
    lines: Sequence[tuple[Point, Point]] = (),
) -> dict[int, Rect] | None:
    """A rectangle for each member's label, stacked beside the crowd, or None.

    The column goes to the right of the crowd if it fits, else to the left,
    first level with the crowd and then stepped up and down. A column fits
    when it is in ``bounds``, overlaps nothing in ``taken``, covers no point
    in ``keep_clear`` and is crossed by no segment in ``lines``: the rays
    leave the crowd as a fan, and a column across the fan is read through
    it. On the right the labels are left-aligned, on the left right-aligned,
    so each sits against the side its leader comes from.
    """
    rows = sorted(members, key=lambda i: tips[i][1])
    width = max(sizes[i][0] for i in rows)
    height = sum(sizes[i][1] for i in rows) + ROW_GAP * (len(rows) - 1)
    xs = [tips[i][0] for i in rows]
    middle = sum(tips[i][1] for i in rows) / len(rows)
    shifts = [0.0]
    for step in (0.5, 1.0, 1.5, 2.0, 3.0):
        shifts += [-step * height, step * height]
    for side in (1, -1):
        left = max(xs) + GAP_PX if side > 0 else min(xs) - GAP_PX - width
        for shift in shifts:
            top = middle - height / 2 + shift
            top = min(max(top, bounds[1]), bounds[1] + bounds[3] - height)
            placed: dict[int, Rect] = {}
            y = top
            for i in rows:
                w, h = sizes[i]
                x = left if side > 0 else left + width - w
                placed[i] = (x, y, w, h)
                y += h + ROW_GAP
            boxes = list(placed.values())
            if all(inside(b, bounds) for b in boxes) \
                    and not any(overlaps(b, t) for b in boxes for t in taken) \
                    and not any(covers(b, p) for b in boxes for p in keep_clear) \
                    and not any(crosses(b, *line) for b in boxes for line in lines):
                return placed
    return None
