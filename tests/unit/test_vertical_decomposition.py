"""Phase 3 — 2D vertical decomposition. Design tests 9–12."""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest

from tests.oracles.simple_arrangement import random_simple_lines
from tests.oracles.vd2d import assert_grid_partition, assert_point_in_exactly_one_open_cell
from vd3d.arrangement2d import build_line_arrangement, verify_arrangement_invariants
from vd3d.geometry import Line2D, Point2D, as_scalar
from vd3d.vertical_decomposition import (
    UNBOUNDED,
    collect_decomposition_walls,
    compute_vertical_decomposition,
    first_hit,
    verify_vd_invariants,
    vertical_rays_from,
)
from vd3d.vertical_decomposition.invariants import point_in_cell

_FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "triangle_vd.json"


def _x_eq_0() -> Line2D:
    return Line2D(a=1, b=0, c=0, id=0)


def _y_eq_0() -> Line2D:
    return Line2D(a=0, b=1, c=0, id=1)


def _x_plus_y_eq_1() -> Line2D:
    return Line2D(a=1, b=1, c=-1, id=2)


def _two_axes():
    return [_x_eq_0(), _y_eq_0()]


def _two_diagonals():
    return [Line2D(a=1, b=-1, c=0, id=0), Line2D(a=1, b=1, c=0, id=1)]  # y=x, y=-x


def _triangle():
    return [_x_eq_0(), _y_eq_0(), _x_plus_y_eq_1()]


def _vertex_at(arrangement, point: Point2D):
    for vertex in arrangement.vertices:
        if vertex.point == point:
            return vertex
    raise AssertionError(f"no vertex at {point}")


# ---------------------------------------------------------------------------
# Step 3.1 — vertical rays and first hit
# ---------------------------------------------------------------------------


def test_single_horizontal_line_no_rays():
    """Design Test 9: one line, no vertices, trivial two-cell split."""
    arr = build_line_arrangement([Line2D(a=0, b=1, c=0)])
    assert arr.vertices == ()
    walls = collect_decomposition_walls(arr)
    assert walls == []
    vd = compute_vertical_decomposition(arr)
    verify_vd_invariants(vd)
    assert len(vd.cells) == 2
    assert all(cell.left_x is None and cell.right_x is None for cell in vd.cells)
    bands = {(cell.lower_line, cell.upper_line) for cell in vd.cells}
    assert bands == {(None, 0), (0, None)}
    assert_grid_partition(vd, half=3)


def test_single_vertical_line_two_halfplanes():
    """Design Test 9: vertical line is itself the wall."""
    arr = build_line_arrangement([_x_eq_0()])
    vd = compute_vertical_decomposition(arr)
    verify_vd_invariants(vd)
    assert len(vd.cells) == 2
    assert len(vd.walls) == 1
    assert vd.walls[0].x == 0
    assert vd.walls[0].y_min is None and vd.walls[0].y_max is None
    xs = {cell.left_x for cell in vd.cells} | {cell.right_x for cell in vd.cells}
    assert xs == {None, 0}
    assert_grid_partition(vd, half=3)


def test_two_axes_rays_from_origin_miss():
    """Design Test 10: from the crossing, ±y rays miss."""
    arr = build_line_arrangement(_two_axes())
    origin = _vertex_at(arr, Point2D(0, 0))
    plus, minus = vertical_rays_from(origin)
    assert plus.direction == 1 and minus.direction == -1
    assert plus.origin == Point2D(0, 0)
    assert first_hit(arr, plus) is UNBOUNDED
    assert first_hit(arr, minus) is UNBOUNDED


def test_two_diagonals_rays_from_origin_miss():
    """Design Test 10: incident lines are not hits above/below."""
    arr = build_line_arrangement(_two_diagonals())
    origin = _vertex_at(arr, Point2D(0, 0))
    plus, minus = vertical_rays_from(origin)
    assert first_hit(arr, plus) is UNBOUNDED
    assert first_hit(arr, minus) is UNBOUNDED
    walls = collect_decomposition_walls(arr)
    assert len(walls) == 2
    xs = {w.x for w in walls}
    assert xs == {0}


def test_triangle_hits_along_vertical_side():
    arr = build_line_arrangement(_triangle())
    v00 = _vertex_at(arr, Point2D(0, 0))
    v01 = _vertex_at(arr, Point2D(0, 1))
    v10 = _vertex_at(arr, Point2D(1, 0))
    plus00, minus00 = vertical_rays_from(v00)
    hit_up = first_hit(arr, plus00)
    assert not hit_up.unbounded
    assert hit_up.point == Point2D(0, 1)
    assert hit_up.vertex_id == v01.id
    assert first_hit(arr, minus00) is UNBOUNDED

    plus01, minus01 = vertical_rays_from(v01)
    assert first_hit(arr, plus01) is UNBOUNDED
    hit_down = first_hit(arr, minus01)
    assert not hit_down.unbounded
    assert hit_down.point == Point2D(0, 0)
    assert hit_down.vertex_id == v00.id

    plus10, minus10 = vertical_rays_from(v10)
    assert first_hit(arr, plus10) is UNBOUNDED
    assert first_hit(arr, minus10) is UNBOUNDED


def test_first_hit_interior_of_edge():
    """From (1,1) in y=0, y=x, y=-x+2, -y hits y=0 at (1,0), not a vertex."""
    lines = [
        Line2D(a=0, b=1, c=0, id=0),
        Line2D(a=1, b=-1, c=0, id=1),
        Line2D(a=1, b=1, c=-2, id=2),
    ]
    arr = build_line_arrangement(lines)
    apex = _vertex_at(arr, Point2D(1, 1))
    plus, minus = vertical_rays_from(apex)
    assert first_hit(arr, plus) is UNBOUNDED
    hit = first_hit(arr, minus)
    assert not hit.unbounded
    assert hit.point == Point2D(1, 0)
    assert hit.vertex_id is None
    assert hit.line_index == 0
    assert hit.edge_id is not None


# ---------------------------------------------------------------------------
# Step 3.2 — trapezoid cells
# ---------------------------------------------------------------------------


def test_two_axes_four_cells():
    """Design Test 10: two axes → 4 unbounded cells, one vertical wall each."""
    arr = build_line_arrangement(_two_axes())
    vd = compute_vertical_decomposition(arr)
    verify_vd_invariants(vd)
    assert len(vd.cells) == 4
    assert all(cell.unbounded for cell in vd.cells)
    for cell in vd.cells:
        assert len(cell.vertical_walls) <= 4
        assert len(cell.vertical_walls) == 1
    assert_grid_partition(vd, half=4)
    assert not any(
        point_in_cell(vd, a, b.representative, closed=False)
        for a in vd.cells
        for b in vd.cells
        if a.id != b.id
    )


def test_two_diagonals_six_cells():
    """Generic two lines: the ±y walls through the crossing split 4 faces into 6."""
    arr = build_line_arrangement(_two_diagonals())
    vd = compute_vertical_decomposition(arr)
    verify_vd_invariants(vd)
    assert len(vd.cells) == 6
    assert len(vd.walls) == 2
    assert_grid_partition(vd, half=4)


def test_empty_arrangement_one_cell():
    arr = build_line_arrangement([])
    vd = compute_vertical_decomposition(arr)
    verify_vd_invariants(vd)
    assert len(vd.cells) == 1
    assert vd.cells[0].left_x is None
    assert vd.cells[0].right_x is None
    assert vd.cells[0].lower_line is None
    assert vd.cells[0].upper_line is None
    assert vd.walls == ()


def test_two_parallels_three_cells():
    arr = build_line_arrangement([Line2D(a=0, b=1, c=0), Line2D(a=0, b=1, c=-2)])
    vd = compute_vertical_decomposition(arr)
    verify_vd_invariants(vd)
    assert len(vd.cells) == 3
    assert all(cell.left_x is None and cell.right_x is None for cell in vd.cells)
    assert_grid_partition(vd, half=4)


def test_cell_union_covers_representatives_and_grid():
    arr = build_line_arrangement(_triangle())
    vd = compute_vertical_decomposition(arr)
    verify_vd_invariants(vd)
    assert_grid_partition(vd, half=4)
    for face in arr.faces:
        assert_point_in_exactly_one_open_cell(vd, face.representative)


# ---------------------------------------------------------------------------
# Step 3.3 — triangle fixture + random
# ---------------------------------------------------------------------------


def _scalar_json(value):
    if value is None:
        return None
    if value.denominator == 1:
        return int(value)
    return f"{value.numerator}/{value.denominator}"


def _load_triangle_fixture() -> dict:
    return json.loads(_FIXTURE.read_text(encoding="utf-8"))


def test_triangle_golden_fixture():
    """Design Test 11: frozen cell count, walls, adjacency."""
    arr = build_line_arrangement(_triangle())
    verify_arrangement_invariants(arr)
    vd = compute_vertical_decomposition(arr)
    verify_vd_invariants(vd)
    gold = _load_triangle_fixture()
    assert len(vd.cells) == gold["cell_count"]
    assert len(vd.walls) == gold["wall_count"]

    for expected, cell in zip(gold["cells"], vd.cells, strict=True):
        assert cell.id == expected["id"]
        assert _scalar_json(cell.left_x) == expected["left_x"]
        assert _scalar_json(cell.right_x) == expected["right_x"]
        assert cell.lower_line == expected["lower_line"]
        assert cell.upper_line == expected["upper_line"]
        assert list(cell.neighbors) == expected["neighbors"]
        assert list(cell.vertical_walls) == expected["vertical_walls"]

    for expected, wall in zip(gold["walls"], vd.walls, strict=True):
        assert wall.id == expected["id"]
        assert _scalar_json(wall.x) == expected["x"]
        assert _scalar_json(wall.y_min) == expected["y_min"]
        assert _scalar_json(wall.y_max) == expected["y_max"]

    triangle_cells = [
        cell
        for cell in vd.cells
        if cell.lower_line is not None
        and cell.upper_line is not None
        and cell.left_x == 0
        and cell.right_x == 1
    ]
    assert len(triangle_cells) == 1
    interior = triangle_cells[0]
    assert point_in_cell(vd, interior, Point2D("1/3", "1/3"), closed=False)
    assert_grid_partition(vd, half=4, step=as_scalar("1/2"))


@pytest.mark.parametrize("n", [3, 4, 5, 6])
@pytest.mark.parametrize("seed", range(5))
def test_random_simple_vd_partition(n: int, seed: int):
    """Design Test 12: random simple arrangements, invariants + unique cell."""
    rng = random.Random(seed)
    lines = random_simple_lines(rng, n)
    arr = build_line_arrangement(lines)
    verify_arrangement_invariants(arr)
    vd = compute_vertical_decomposition(arr)
    verify_vd_invariants(vd)
    assert all(len(cell.vertical_walls) <= 4 for cell in vd.cells)
    assert_grid_partition(vd, half=5, step=1)
    for cell in vd.cells:
        others = [
            other
            for other in vd.cells
            if other.id != cell.id
            and point_in_cell(vd, other, cell.representative, closed=False)
        ]
        assert others == []
