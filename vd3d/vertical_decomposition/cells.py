"""Build trapezoid cells of a 2D vertical decomposition.

Walls are y-parallel (``x = const``). A cell is a (possibly unbounded)
trapezoid: left/right ``x`` or unbounded, lower/upper supporting line or
unbounded.
"""

from __future__ import annotations

from vd3d.arrangement2d.invariants import point_in_face
from vd3d.arrangement2d.types import Arrangement2D, Vertex
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar
from vd3d.vertical_decomposition.geom import (
    line_is_vertical,
    x_of_vertical_line,
    y_at_x,
)
from vd3d.vertical_decomposition.rays import first_hit, vertical_rays_from
from vd3d.vertical_decomposition.types import (
    Hit,
    VDCell2D,
    VerticalDecomposition,
    VerticalRay,
    VerticalWall,
)


def compute_vertical_decomposition(
    arrangement: Arrangement2D,
) -> VerticalDecomposition:
    """``COMPUTE_VERTICAL_DECOMPOSITION``: walls from vertices, then trapezoids."""
    walls = collect_decomposition_walls(arrangement)
    cells = build_cells(arrangement, walls)
    return VerticalDecomposition(
        arrangement=arrangement, cells=tuple(cells), walls=tuple(walls)
    )


def collect_decomposition_walls(
    arrangement: Arrangement2D,
) -> list[VerticalWall]:
    """Shoot ±y rays from every vertex and insert the resulting walls."""
    walls: list[VerticalWall] = []
    for vertex in arrangement.vertices:
        for ray in vertical_rays_from(vertex):
            hit = first_hit(arrangement, ray)
            insert_decomposition_segment(walls, vertex, ray, hit)
    return finalize_decomposition_walls(arrangement, walls)


def finalize_decomposition_walls(
    arrangement: Arrangement2D, walls: list[VerticalWall]
) -> list[VerticalWall]:
    """Add vertical-input-line walls, sort, and assign dense ids."""
    _ensure_vertical_line_walls(arrangement, walls)
    walls.sort(key=_wall_sort_key)
    for i, wall in enumerate(walls):
        wall.id = i
    return walls


def insert_decomposition_segment(
    walls: list[VerticalWall],
    vertex: Vertex,
    ray: VerticalRay,
    hit: Hit,
) -> VerticalWall:
    """``INSERT_DECOMPOSITION_SEGMENT``: the wall from ``vertex`` along ``ray``."""
    x = vertex.point.x
    if ray.direction > 0:
        y_min = vertex.point.y
        y_max = None if hit.unbounded else hit.point.y  # type: ignore[union-attr]
        bottom_vertex_id = vertex.id
        top_vertex_id = None if hit.unbounded else hit.vertex_id
        bottom_line_index = None
        top_line_index = None if hit.unbounded else hit.line_index
    else:
        y_max = vertex.point.y
        y_min = None if hit.unbounded else hit.point.y  # type: ignore[union-attr]
        top_vertex_id = vertex.id
        bottom_vertex_id = None if hit.unbounded else hit.vertex_id
        top_line_index = None
        bottom_line_index = None if hit.unbounded else hit.line_index

    for existing in walls:
        if (
            existing.x == x
            and existing.y_min == y_min
            and existing.y_max == y_max
        ):
            return existing

    wall = VerticalWall(
        id=len(walls),
        x=x,
        y_min=y_min,
        y_max=y_max,
        bottom_vertex_id=bottom_vertex_id,
        top_vertex_id=top_vertex_id,
        bottom_line_index=bottom_line_index,
        top_line_index=top_line_index,
    )
    walls.append(wall)
    return wall


def build_cells(
    arrangement: Arrangement2D, walls: list[VerticalWall]
) -> list[VDCell2D]:
    """``BUILD_CELLS`` / ``MERGE_OR_CREATE_DECOMPOSITION_CELLS``.

    Split the plane into open vertical strips at vertex ``x``-coordinates
    (and vertical input lines). Inside a strip the non-vertical lines do
    not cross, so they sort by ``y`` into trapezoids.
    """
    xs = _strip_xs(arrangement)
    strips = _strips(xs)
    lines = arrangement.lines
    cells: list[VDCell2D] = []

    for left_x, right_x in strips:
        sample = _sample_x(left_x, right_x)
        ordered = _non_vertical_sorted_at(lines, sample)
        bands: list[tuple[int | None, int | None]] = []
        if not ordered:
            bands.append((None, None))
        else:
            bands.append((None, ordered[0]))
            for lower, upper in zip(ordered, ordered[1:]):
                bands.append((lower, upper))
            bands.append((ordered[-1], None))

        for lower_line, upper_line in bands:
            representative = _representative(
                left_x, right_x, lower_line, upper_line, lines
            )
            source_face = _source_face(arrangement, representative)
            cells.append(
                VDCell2D(
                    id=len(cells),
                    left_x=left_x,
                    right_x=right_x,
                    lower_line=lower_line,
                    upper_line=upper_line,
                    representative=representative,
                    source_face=source_face,
                )
            )

    cells = _merge_unwalled_strips(cells, walls, arrangement)
    _attach_walls(cells, walls, lines)
    _attach_neighbors(cells, lines)
    return cells


def _ensure_vertical_line_walls(
    arrangement: Arrangement2D, walls: list[VerticalWall]
) -> None:
    """A vertical input line with no vertices is still a full-line wall."""
    covered_x = {wall.x for wall in walls}
    for line in arrangement.lines:
        if not line_is_vertical(line):
            continue
        x = x_of_vertical_line(line)
        if x in covered_x:
            continue
        if any(vertex.point.x == x for vertex in arrangement.vertices):
            continue
        walls.append(
            VerticalWall(
                id=len(walls),
                x=x,
                y_min=None,
                y_max=None,
            )
        )
        covered_x.add(x)


def _wall_sort_key(
    wall: VerticalWall,
) -> tuple[Scalar, int, Scalar, int, Scalar]:
    """``x`` ascending; unbounded ``y_min`` first; unbounded ``y_max`` last."""
    y_min_flag = 0 if wall.y_min is None else 1
    y_min_val = 0 if wall.y_min is None else wall.y_min
    y_max_flag = 1 if wall.y_max is None else 0
    y_max_val = 0 if wall.y_max is None else wall.y_max
    return (wall.x, y_min_flag, y_min_val, y_max_flag, y_max_val)


def _strip_xs(arrangement: Arrangement2D) -> list[Scalar]:
    xs = {vertex.point.x for vertex in arrangement.vertices}
    for line in arrangement.lines:
        if line_is_vertical(line):
            xs.add(x_of_vertical_line(line))
    return sorted(xs)


def _strips(
    xs: list[Scalar],
) -> list[tuple[Scalar | None, Scalar | None]]:
    if not xs:
        return [(None, None)]
    strips: list[tuple[Scalar | None, Scalar | None]] = [(None, xs[0])]
    for left, right in zip(xs, xs[1:]):
        strips.append((left, right))
    strips.append((xs[-1], None))
    return strips


def _sample_x(left: Scalar | None, right: Scalar | None) -> Scalar:
    if left is None and right is None:
        return Scalar(0)
    if left is None:
        return right - 1  # type: ignore[operator]
    if right is None:
        return left + 1
    return (left + right) / 2


def _non_vertical_sorted_at(
    lines: tuple[Line2D, ...], x: Scalar
) -> list[int]:
    scored: list[tuple[Scalar, int]] = []
    for i, line in enumerate(lines):
        if line_is_vertical(line):
            continue
        scored.append((y_at_x(line, x), i))
    scored.sort()
    ys = [y for y, _ in scored]
    if len(ys) != len(set(ys)):
        raise RuntimeError(f"non-vertical lines are not strictly ordered at x={x}")
    return [i for _, i in scored]


def _representative(
    left_x: Scalar | None,
    right_x: Scalar | None,
    lower_line: int | None,
    upper_line: int | None,
    lines: tuple[Line2D, ...],
) -> Point2D:
    x = _sample_x(left_x, right_x)
    if lower_line is not None and upper_line is not None:
        y = (y_at_x(lines[lower_line], x) + y_at_x(lines[upper_line], x)) / 2
    elif lower_line is not None:
        y = y_at_x(lines[lower_line], x) + 1
    elif upper_line is not None:
        y = y_at_x(lines[upper_line], x) - 1
    else:
        y = Scalar(0)
    return Point2D(x, y)


def _source_face(arrangement: Arrangement2D, point: Point2D) -> int:
    matches = [
        face.id
        for face in arrangement.faces
        if point_in_face(arrangement, face, point, closed=False)
    ]
    if len(matches) != 1:
        raise RuntimeError(
            f"representative {point} lies in {len(matches)} arrangement faces"
        )
    return matches[0]


def y_range_at(
    cell: VDCell2D, x: Scalar, lines: tuple[Line2D, ...]
) -> tuple[Scalar | None, Scalar | None]:
    """Open ``y``-interval of ``cell`` on the vertical line at ``x``."""
    lo = y_at_x(lines[cell.lower_line], x) if cell.lower_line is not None else None
    hi = y_at_x(lines[cell.upper_line], x) if cell.upper_line is not None else None
    return lo, hi


def interval_overlap_open(
    lo1: Scalar | None,
    hi1: Scalar | None,
    lo2: Scalar | None,
    hi2: Scalar | None,
) -> bool:
    """True iff the open intervals ``(lo1, hi1)`` and ``(lo2, hi2)`` overlap.

    ``None`` on the left is ``-∞``; ``None`` on the right is ``+∞``.
    """
    lo = _max_lo(lo1, lo2)
    hi = _min_hi(hi1, hi2)
    if lo is None or hi is None:
        return True
    return lo < hi


def interval_positive(lo: Scalar | None, hi: Scalar | None) -> bool:
    if lo is None or hi is None:
        return True
    return lo < hi


def _max_lo(a: Scalar | None, b: Scalar | None) -> Scalar | None:
    if a is None:
        return b
    if b is None:
        return a
    return a if a > b else b


def _min_hi(a: Scalar | None, b: Scalar | None) -> Scalar | None:
    if a is None:
        return b
    if b is None:
        return a
    return a if a < b else b


def _wall_covers(
    walls: list[VerticalWall],
    x: Scalar,
    lo: Scalar | None,
    hi: Scalar | None,
) -> bool:
    """True iff some wall at ``x`` overlaps the open y-interval ``(lo, hi)``."""
    if not interval_positive(lo, hi):
        return False
    for wall in walls:
        if wall.x != x:
            continue
        if interval_overlap_open(lo, hi, wall.y_min, wall.y_max):
            return True
    return False


def _merge_unwalled_strips(
    cells: list[VDCell2D],
    walls: list[VerticalWall],
    arrangement: Arrangement2D,
) -> list[VDCell2D]:
    """Merge left/right strip cells in the same band when no wall sits between them.

    Vertex ``x``-coordinates over-split bands that a ray never reaches.
    ``MERGE_OR_CREATE_DECOMPOSITION_CELLS`` glues those pieces back together.
    """
    if not cells:
        return cells

    lines = arrangement.lines
    parent = list(range(len(cells)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    by_band: dict[tuple[int | None, int | None], list[VDCell2D]] = {}
    for cell in cells:
        by_band.setdefault((cell.lower_line, cell.upper_line), []).append(cell)

    for group in by_band.values():
        group.sort(key=lambda c: _strip_sort_key((c.left_x, c.right_x)))
        for left_cell, right_cell in zip(group, group[1:]):
            if left_cell.right_x is None or left_cell.right_x != right_cell.left_x:
                continue
            wall_x = left_cell.right_x
            lo, hi = y_range_at(left_cell, wall_x, lines)
            if _wall_covers(walls, wall_x, lo, hi):
                continue
            union(left_cell.id, right_cell.id)

    groups: dict[int, list[VDCell2D]] = {}
    for cell in cells:
        groups.setdefault(find(cell.id), []).append(cell)

    merged: list[VDCell2D] = []
    for group in groups.values():
        group.sort(key=lambda c: _strip_sort_key((c.left_x, c.right_x)))
        left_x = group[0].left_x
        right_x = group[-1].right_x
        lower_line = group[0].lower_line
        upper_line = group[0].upper_line
        representative = _representative(left_x, right_x, lower_line, upper_line, lines)
        merged.append(
            VDCell2D(
                id=len(merged),
                left_x=left_x,
                right_x=right_x,
                lower_line=lower_line,
                upper_line=upper_line,
                representative=representative,
                source_face=_source_face(arrangement, representative),
            )
        )

    merged.sort(key=lambda c: (_strip_sort_key((c.left_x, c.right_x)), _band_sort_key(c, lines)))
    for i, cell in enumerate(merged):
        cell.id = i
    return merged


def _attach_walls(
    cells: list[VDCell2D],
    walls: list[VerticalWall],
    lines: tuple[Line2D, ...],
) -> None:
    for cell in cells:
        ids: list[int] = []
        for wall in walls:
            if cell.left_x is not None and wall.x == cell.left_x:
                lo, hi = y_range_at(cell, wall.x, lines)
                if interval_positive(lo, hi) and interval_overlap_open(
                    lo, hi, wall.y_min, wall.y_max
                ):
                    ids.append(wall.id)
            elif cell.right_x is not None and wall.x == cell.right_x:
                lo, hi = y_range_at(cell, wall.x, lines)
                if interval_positive(lo, hi) and interval_overlap_open(
                    lo, hi, wall.y_min, wall.y_max
                ):
                    ids.append(wall.id)
        cell.vertical_walls = tuple(ids)


def _attach_neighbors(cells: list[VDCell2D], lines: tuple[Line2D, ...]) -> None:
    n = len(cells)
    adj: list[set[int]] = [set() for _ in range(n)]

    for i, a in enumerate(cells):
        for b in cells[i + 1 :]:
            if _cells_adjacent(a, b, lines):
                adj[a.id].add(b.id)
                adj[b.id].add(a.id)

    for cell in cells:
        cell.neighbors = tuple(sorted(adj[cell.id]))


def _cells_adjacent(
    a: VDCell2D, b: VDCell2D, lines: tuple[Line2D, ...]
) -> bool:
    if a.right_x is not None and a.right_x == b.left_x:
        lo_a, hi_a = y_range_at(a, a.right_x, lines)
        lo_b, hi_b = y_range_at(b, b.left_x, lines)
        if interval_positive(lo_a, hi_a) and interval_positive(lo_b, hi_b):
            if interval_overlap_open(lo_a, hi_a, lo_b, hi_b):
                return True
    if b.right_x is not None and b.right_x == a.left_x:
        lo_a, hi_a = y_range_at(a, a.left_x, lines)
        lo_b, hi_b = y_range_at(b, b.right_x, lines)
        if interval_positive(lo_a, hi_a) and interval_positive(lo_b, hi_b):
            if interval_overlap_open(lo_a, hi_a, lo_b, hi_b):
                return True
    share_line = (
        (a.upper_line is not None and a.upper_line == b.lower_line)
        or (a.lower_line is not None and a.lower_line == b.upper_line)
    )
    if share_line and interval_overlap_open(a.left_x, a.right_x, b.left_x, b.right_x):
        return True
    return False


def _band_sort_key(
    cell: VDCell2D, lines: tuple[Line2D, ...]
) -> tuple[int, Scalar]:
    x = _sample_x(cell.left_x, cell.right_x)
    if cell.lower_line is not None:
        return (1, y_at_x(lines[cell.lower_line], x))
    return (0, Scalar(0))


def _strip_sort_key(
    key: tuple[Scalar | None, Scalar | None],
) -> tuple[int, Scalar]:
    left, _right = key
    if left is None:
        return (0, Scalar(0))
    return (1, left)
