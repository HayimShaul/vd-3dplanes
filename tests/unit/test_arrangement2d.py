"""Phase 2 — 2D line arrangement. Design tests 5–8."""

from __future__ import annotations

import random

import pytest

from tests.oracles.simple_arrangement import random_simple_lines, simple_arrangement_counts
from vd3d.arrangement2d import (
    COINCIDENT,
    PARALLEL,
    build_line_arrangement,
    cmp_direction_ccw,
    compute_edge_pieces,
    compute_vertices,
    intersect_lines_2d,
    line_direction,
    sort_directions_ccw,
    verify_arrangement_invariants,
)
from vd3d.arrangement2d.invariants import (
    euler_characteristic_plane,
    point_in_face,
    twin_involution,
)
from vd3d.arrangement2d.predicates import quadrant
from vd3d.geometry import Line2D, Point2D


def _x_eq_0() -> Line2D:
    return Line2D(a=1, b=0, c=0, id=0)


def _y_eq_0() -> Line2D:
    return Line2D(a=0, b=1, c=0, id=1)


def _x_plus_y_eq_0() -> Line2D:
    return Line2D(a=1, b=1, c=0, id=2)


def _x_plus_y_eq_1() -> Line2D:
    return Line2D(a=1, b=1, c=-1, id=2)


def _two_axes():
    return [_x_eq_0(), _y_eq_0()]


def _three_concurrent():
    return [_x_eq_0(), _y_eq_0(), _x_plus_y_eq_0()]


def _triangle():
    return [_x_eq_0(), _y_eq_0(), _x_plus_y_eq_1()]


def _piece_counts(arrangement):
    rays = sum(1 for e in arrangement.edges if e.kind == "ray")
    segments = sum(1 for e in arrangement.edges if e.kind == "segment")
    whole = sum(1 for e in arrangement.edges if e.kind == "line")
    return rays, segments, whole


# ---------------------------------------------------------------------------
# Step 2.1 — line–line intersection
# ---------------------------------------------------------------------------


def test_intersect_two_crossing_axes():
    point = intersect_lines_2d(_x_eq_0(), _y_eq_0())
    assert point == Point2D(0, 0)
    assert _x_eq_0().contains(point)
    assert _y_eq_0().contains(point)


def test_intersect_generic_crossing_exact():
    left = Line2D(a=1, b=1, c=-3)  # x + y = 3
    right = Line2D(a=1, b=-1, c=-1)  # x - y = 1
    point = intersect_lines_2d(left, right)
    assert point == Point2D(2, 1)
    assert left.contains(point)
    assert right.contains(point)
    assert left.eval(point) == 0
    assert right.eval(point) == 0


def test_intersect_parallel_distinct():
    left = Line2D(a=1, b=0, c=0)  # x = 0
    right = Line2D(a=1, b=0, c=-2)  # x = 2
    assert intersect_lines_2d(left, right) is PARALLEL


def test_intersect_coincident_scaled():
    left = Line2D(a=1, b=0, c=0)
    right = Line2D(a=2, b=0, c=0)
    assert intersect_lines_2d(left, right) is COINCIDENT


def test_intersect_coincident_opposite_normal():
    left = Line2D(a=1, b=1, c=-1)
    right = Line2D(a=-2, b=-2, c=2)
    assert intersect_lines_2d(left, right) is COINCIDENT


def test_intersect_rejects_float_coefficients_via_line2d():
    with pytest.raises(TypeError):
        Line2D(a=1.0, b=0, c=0)


# ---------------------------------------------------------------------------
# Step 2.2 — vertices and edge pieces
# ---------------------------------------------------------------------------


def test_two_lines_one_vertex_four_rays():
    """Design Test 5."""
    lines = tuple(_two_axes())
    vertices = compute_vertices(lines)
    pieces = compute_edge_pieces(lines, vertices)
    assert len(vertices) == 1
    assert vertices[0].point == Point2D(0, 0)
    assert set(vertices[0].line_indices) == {0, 1}
    assert len(pieces) == 4
    assert all(p.kind == "ray" for p in pieces)


def test_three_concurrent_one_shared_vertex():
    """Design Test 6: one vertex, not three copies."""
    lines = tuple(_three_concurrent())
    vertices = compute_vertices(lines)
    assert len(vertices) == 1
    assert vertices[0].point == Point2D(0, 0)
    assert set(vertices[0].line_indices) == {0, 1, 2}
    pieces = compute_edge_pieces(lines, vertices)
    assert len(pieces) == 6
    assert all(p.kind == "ray" for p in pieces)


def test_three_lines_general_vertices_and_pieces():
    """Design Test 7."""
    lines = tuple(_triangle())
    vertices = compute_vertices(lines)
    pieces = compute_edge_pieces(lines, vertices)
    assert len(vertices) == 3
    points = {v.point for v in vertices}
    assert points == {Point2D(0, 0), Point2D(1, 0), Point2D(0, 1)}
    rays = [p for p in pieces if p.kind == "ray"]
    segments = [p for p in pieces if p.kind == "segment"]
    assert len(rays) == 6
    assert len(segments) == 3


def test_coincident_lines_rejected_by_vertices():
    with pytest.raises(ValueError, match="coincident"):
        compute_vertices((Line2D(a=1, b=0, c=0), Line2D(a=2, b=0, c=0)))


def test_single_line_is_one_whole_piece():
    vertices = compute_vertices((Line2D(a=0, b=1, c=-1),))
    pieces = compute_edge_pieces((Line2D(a=0, b=1, c=-1),), vertices)
    assert vertices == ()
    assert len(pieces) == 1
    assert pieces[0].kind == "line"


def test_two_parallels_no_vertices_two_whole_lines():
    lines = (Line2D(a=0, b=1, c=0), Line2D(a=0, b=1, c=-1))
    vertices = compute_vertices(lines)
    pieces = compute_edge_pieces(lines, vertices)
    assert vertices == ()
    assert len(pieces) == 2
    assert all(p.kind == "line" for p in pieces)


# ---------------------------------------------------------------------------
# Step 2.3 — vertex circulation
# ---------------------------------------------------------------------------


def test_twin_involution_two_axes():
    arr = build_line_arrangement(_two_axes())
    assert twin_involution(arr)
    for he in arr.half_edges:
        assert arr.half_edges[arr.half_edges[he.twin_id].twin_id] is he or (
            arr.half_edges[he.twin_id].twin_id == he.id
        )


def test_outgoing_ccw_at_origin_two_axes():
    arr = build_line_arrangement(_two_axes())
    vertex = arr.vertices[0]
    assert vertex.point == Point2D(0, 0)
    dirs = [arr.half_edges[hid].direction for hid in vertex.outgoing]
    assert len(dirs) == 4
    assert dirs == sort_directions_ccw(list(dirs))
    x_dir = line_direction(_x_eq_0())  # (0, 1)
    y_dir = line_direction(_y_eq_0())  # (-1, 0)
    expected = sort_directions_ccw(
        [x_dir, (-x_dir[0], -x_dir[1]), y_dir, (-y_dir[0], -y_dir[1])]
    )
    assert dirs == expected
    assert [quadrant(d) for d in dirs] == [0, 1, 2, 3]


def test_outgoing_ccw_generic_angles():
    lines = [
        Line2D(a=1, b=0, c=0),  # x = 0, dir (0, 1)
        Line2D(a=1, b=1, c=0),  # x + y = 0, dir (-1, 1)
    ]
    arr = build_line_arrangement(lines)
    dirs = [arr.half_edges[hid].direction for hid in arr.vertices[0].outgoing]
    assert len(dirs) == 4
    assert dirs == sort_directions_ccw(list(dirs))
    for i in range(len(dirs) - 1):
        assert cmp_direction_ccw(dirs[i], dirs[i + 1]) < 0


def test_quadrant_order_from_positive_x():
    ordered = sort_directions_ccw([(1, 0), (0, -1), (-1, 0), (0, 1)])
    assert [quadrant(d) for d in ordered] == [0, 1, 2, 3]
    assert ordered[0] == (1, 0)
    assert ordered[1] == (0, 1)
    assert ordered[2] == (-1, 0)
    assert ordered[3] == (0, -1)


# ---------------------------------------------------------------------------
# Step 2.4 — faces
# ---------------------------------------------------------------------------


def test_two_lines_four_unbounded_faces():
    """Design Test 5: 4 unbounded faces / sectors."""
    arr = build_line_arrangement(_two_axes())
    verify_arrangement_invariants(arr)
    assert len(arr.vertices) == 1
    assert _piece_counts(arr) == (4, 0, 0)
    assert len(arr.faces) == 4
    assert len(arr.bounded_faces) == 0
    assert len(arr.unbounded_faces) == 4
    assert euler_characteristic_plane(arr)


def test_triangle_one_bounded_six_unbounded():
    """Design Test 7: 1 bounded triangle + 6 unbounded."""
    arr = build_line_arrangement(_triangle())
    verify_arrangement_invariants(arr)
    assert len(arr.vertices) == 3
    assert _piece_counts(arr) == (6, 3, 0)
    assert len(arr.bounded_faces) == 1
    assert len(arr.unbounded_faces) == 6
    assert len(arr.faces) == 7
    bounded = arr.bounded_faces[0]
    cycle = arr.cycle(bounded)
    assert len(cycle) == 3
    assert all(arr.edges[he.edge_id].kind == "segment" for he in cycle)
    # The triangle's representative sits inside x>0, y>0, x+y<1.
    p = bounded.representative
    assert p.x > 0 and p.y > 0 and p.x + p.y < 1
    assert point_in_face(arr, bounded, Point2D("1/3", "1/3"), closed=False)


def test_three_concurrent_six_unbounded():
    arr = build_line_arrangement(_three_concurrent())
    verify_arrangement_invariants(arr)
    assert len(arr.vertices) == 1
    assert len(arr.faces) == 6
    assert len(arr.bounded_faces) == 0


def test_single_line_two_halfplanes():
    arr = build_line_arrangement([Line2D(a=0, b=1, c=0)])
    verify_arrangement_invariants(arr)
    assert len(arr.vertices) == 0
    assert _piece_counts(arr) == (0, 0, 1)
    assert len(arr.faces) == 2
    assert all(f.unbounded for f in arr.faces)


def test_two_parallels_three_faces():
    arr = build_line_arrangement([Line2D(a=0, b=1, c=0), Line2D(a=0, b=1, c=-2)])
    verify_arrangement_invariants(arr)
    assert len(arr.vertices) == 0
    assert _piece_counts(arr) == (0, 0, 2)
    assert len(arr.faces) == 3
    assert all(f.unbounded for f in arr.faces)


def test_empty_arrangement_one_unbounded_face():
    arr = build_line_arrangement([])
    assert len(arr.faces) == 1
    assert arr.faces[0].unbounded
    assert euler_characteristic_plane(arr)


def test_vertex_on_incident_and_edge_on_line():
    arr = build_line_arrangement(_triangle())
    for vertex in arr.vertices:
        for hid in vertex.outgoing:
            he = arr.half_edges[hid]
            assert arr.lines[he.line_index].contains(vertex.point)
    for edge in arr.edges:
        line = arr.lines[edge.line_index]
        for vid in (edge.start_vertex, edge.end_vertex):
            if vid is not None:
                assert line.contains(arr.vertices[vid].point)


# ---------------------------------------------------------------------------
# Step 2.5 — randomized formula oracle
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", range(2, 13))
@pytest.mark.parametrize("seed", [0, 1, 7])
def test_random_simple_arrangement_matches_formula(n: int, seed: int):
    """Design Test 8."""
    rng = random.Random(seed)
    lines = random_simple_lines(rng, n)
    arr = build_line_arrangement(lines)
    verify_arrangement_invariants(arr)
    v, e, f = simple_arrangement_counts(n)
    assert len(arr.vertices) == v
    assert len(arr.edges) == e
    assert len(arr.faces) == f
    assert len(arr.bounded_faces) + len(arr.unbounded_faces) == f
