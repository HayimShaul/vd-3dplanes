"""Phase 4 — zone of a query line. Design tests 13–14."""

from __future__ import annotations

import random

import pytest

from tests.oracles.simple_arrangement import random_simple_lines
from tests.oracles.zone import (
    brute_force_crossings,
    brute_force_face_ids,
    random_query_missing_vertices,
)
from vd3d.arrangement2d import build_line_arrangement
from vd3d.arrangement2d.predicates import parameter_on_line
from vd3d.geometry import Line2D, Point2D
from vd3d.geometry.scalar import as_scalar
from vd3d.zone import (
    QueryOverlapsArrangement,
    QueryThroughVertex,
    compute_crossings,
    compute_supporting_line_zone,
    compute_zone,
    face_containing_point,
    face_on_other_side,
    point_on_line_at_parameter,
    verify_supporting_line_zone,
    verify_zone_invariants,
)


def _x_eq_0() -> Line2D:
    return Line2D(a=1, b=0, c=0, id=0)


def _y_eq_0() -> Line2D:
    return Line2D(a=0, b=1, c=0, id=1)


def _x_plus_y_eq_1() -> Line2D:
    return Line2D(a=1, b=1, c=-1, id=2)


def _triangle() -> list[Line2D]:
    return [_x_eq_0(), _y_eq_0(), _x_plus_y_eq_1()]


def _two_axes() -> list[Line2D]:
    return [_x_eq_0(), _y_eq_0()]


def _query_y_equals(value: str | int) -> Line2D:
    """Horizontal line ``y = value`` oriented left → right (direction ``(+1, 0)``)."""
    return Line2D(a=0, b=-1, c=value)


def _assert_crossings_match(left, right) -> None:
    assert len(left) == len(right)
    for a, b in zip(left, right):
        assert a.point == b.point
        assert a.t == b.t
        assert a.edge_id == b.edge_id
        assert a.vertex_id == b.vertex_id


# ---------------------------------------------------------------------------
# Step 4.1 — crossings
# ---------------------------------------------------------------------------


def test_point_on_line_at_parameter_roundtrip():
    line = Line2D(a=1, b=1, c=-1)
    for t in (0, 1, -3, "5/2"):
        point = point_on_line_at_parameter(line, t)
        assert line.contains(point)
        assert parameter_on_line(line, point) == as_scalar(t)


def test_point_on_horizontal_parameter_is_x():
    query = _query_y_equals("1/2")
    point = point_on_line_at_parameter(query, 3)
    assert point == Point2D(3, "1/2")
    assert parameter_on_line(query, point) == 3


def test_crossings_triangle_y_half_two_interior_hits():
    arr = build_line_arrangement(_triangle())
    query = _query_y_equals("1/2")
    crossings = compute_crossings(arr, query)
    assert [c.point for c in crossings] == [Point2D(0, "1/2"), Point2D("1/2", "1/2")]
    assert all(c.vertex_id is None for c in crossings)
    assert crossings[0].t < crossings[1].t
    for crossing in crossings:
        assert query.contains(crossing.point)
        edge = arr.edges[crossing.edge_id]
        assert arr.lines[edge.line_index].contains(crossing.point)
    _assert_crossings_match(crossings, brute_force_crossings(arr, query))


def test_crossings_skip_parallel_to_one_line():
    arr = build_line_arrangement(_triangle())
    query = _query_y_equals(-2)
    crossings = compute_crossings(arr, query)
    # Parallel to y=0, so only x=0 and x+y=1.
    assert len(crossings) == 2
    lines_hit = {arr.edges[c.edge_id].line_index for c in crossings}
    assert lines_hit == {0, 2}
    _assert_crossings_match(crossings, brute_force_crossings(arr, query))


def test_crossings_reject_overlap():
    arr = build_line_arrangement(_triangle())
    with pytest.raises(QueryOverlapsArrangement):
        compute_crossings(arr, _y_eq_0())


def test_crossings_empty_arrangement():
    arr = build_line_arrangement([])
    crossings = compute_crossings(arr, _query_y_equals(0))
    assert crossings == ()


def test_crossings_sorted_match_brute_force_random():
    rng = random.Random(4)
    for n in range(2, 7):
        arr = build_line_arrangement(random_simple_lines(rng, n))
        query = random_query_missing_vertices(rng, arr)
        got = compute_crossings(arr, query)
        assert [c.t for c in got] == sorted(c.t for c in got)
        _assert_crossings_match(got, brute_force_crossings(arr, query))
        assert all(c.vertex_id is None for c in got)


# ---------------------------------------------------------------------------
# Step 4.2 — face walk
# ---------------------------------------------------------------------------


def test_triangle_zone_unbounded_interior_unbounded():
    """Design §5.1: line through the triangle."""
    arr = build_line_arrangement(_triangle())
    query = _query_y_equals("1/2")
    zone = compute_zone(arr, query)
    verify_zone_invariants(arr, zone)
    assert len(zone.faces) == 3
    assert len(zone.edges) == 2
    assert zone.vertices == ()
    assert zone.faces[0].unbounded
    assert not zone.faces[1].unbounded
    assert zone.faces[2].unbounded
    interior = arr.bounded_faces[0]
    assert zone.faces[1].id == interior.id
    assert [f.id for f in zone.faces] == list(brute_force_face_ids(arr, query))
    start = point_on_line_at_parameter(query, -1)
    assert start.x < 0
    assert face_containing_point(arr, start).id == zone.faces[0].id


def test_face_on_other_side_is_involution():
    arr = build_line_arrangement(_triangle())
    query = _query_y_equals("1/2")
    zone = compute_zone(arr, query)
    for face, edge in zip(zone.faces, zone.edges):
        other = face_on_other_side(arr, face, edge.id)
        assert face_on_other_side(arr, other, edge.id).id == face.id


def test_single_line_zone_two_halfplanes():
    arr = build_line_arrangement([Line2D(a=1, b=0, c=0, id=0)])
    query = _query_y_equals(0)
    zone = compute_zone(arr, query)
    verify_zone_invariants(arr, zone)
    assert len(zone.crossings) == 1
    assert len(zone.faces) == 2
    assert all(f.unbounded for f in zone.faces)
    assert zone.faces[0].id != zone.faces[1].id


def test_no_crossings_stays_in_one_face():
    arr = build_line_arrangement(
        [Line2D(a=0, b=1, c=0, id=0), Line2D(a=0, b=1, c=-2, id=1)]
    )
    query = _query_y_equals(1)
    zone = compute_zone(arr, query)
    verify_zone_invariants(arr, zone)
    assert zone.crossings == ()
    assert len(zone.faces) == 1
    assert zone.faces[0].unbounded


def test_zone_five_lines_matches_brute_force():
    """Design Test 13."""
    rng = random.Random(13)
    arr = build_line_arrangement(random_simple_lines(rng, 5))
    query = random_query_missing_vertices(rng, arr)
    zone = compute_zone(arr, query)
    verify_zone_invariants(arr, zone)
    _assert_crossings_match(zone.crossings, brute_force_crossings(arr, query))
    assert [f.id for f in zone.faces] == list(brute_force_face_ids(arr, query))
    assert len(zone.faces) == len(zone.crossings) + 1


@pytest.mark.parametrize("n", range(2, 8))
@pytest.mark.parametrize("seed", [0, 1, 7])
def test_random_zone_matches_oracle(n: int, seed: int):
    rng = random.Random(seed)
    arr = build_line_arrangement(random_simple_lines(rng, n))
    query = random_query_missing_vertices(rng, arr)
    zone = compute_zone(arr, query)
    verify_zone_invariants(arr, zone)
    assert [f.id for f in zone.faces] == list(brute_force_face_ids(arr, query))


# ---------------------------------------------------------------------------
# Step 4.3 — boundary cases (well-defined under general position)
# ---------------------------------------------------------------------------


def test_zone_far_from_vertices_below_triangle():
    arr = build_line_arrangement(_triangle())
    query = _query_y_equals(-2)
    zone = compute_zone(arr, query)
    verify_zone_invariants(arr, zone)
    assert all(f.unbounded for f in zone.faces)
    assert arr.bounded_faces[0].id not in {f.id for f in zone.faces}
    assert [f.id for f in zone.faces] == list(brute_force_face_ids(arr, query))
    for crossing in zone.crossings:
        for vertex in arr.vertices:
            assert crossing.point != vertex.point


def test_zone_near_miss_still_enters_triangle():
    """Exact near-miss: y = 1/1000, just above the x-axis."""
    arr = build_line_arrangement(_triangle())
    query = _query_y_equals("1/1000")
    zone = compute_zone(arr, query)
    verify_zone_invariants(arr, zone)
    assert [f.unbounded for f in zone.faces] == [True, False, True]
    assert zone.faces[1].id == arr.bounded_faces[0].id
    assert [f.id for f in zone.faces] == list(brute_force_face_ids(arr, query))
    origin = Point2D(0, 0)
    assert not query.contains(origin)
    for vertex in arr.vertices:
        assert not query.contains(vertex.point)


def test_through_vertex_rejected():
    arr = build_line_arrangement(_triangle())
    query = Line2D(a=1, b=-1, c=0)
    assert any(query.contains(v.point) for v in arr.vertices)
    with pytest.raises(QueryThroughVertex):
        compute_zone(arr, query)
    crossings = compute_crossings(arr, query)
    assert any(c.vertex_id is not None for c in crossings)


def test_overlap_rejected_by_compute_zone():
    arr = build_line_arrangement(_triangle())
    with pytest.raises(QueryOverlapsArrangement):
        compute_zone(arr, Line2D(a=1, b=0, c=0))


# ---------------------------------------------------------------------------
# Step 4.4 — full zone of a supporting line
# ---------------------------------------------------------------------------


def test_supporting_two_lines_no_opposite_vertices():
    arr = build_line_arrangement(_two_axes())
    for line_index in (0, 1):
        zone = compute_supporting_line_zone(arr, line_index)
        verify_supporting_line_zone(arr, zone)
        assert zone.opposite_vertices == ()
        assert len(zone.vertices_on_line) == 1
        assert zone.vertices_on_line[0].point == Point2D(0, 0)
        assert len(zone.faces) == 4


def test_supporting_triangle_opposite_is_the_far_vertex():
    """y=0 is a side of the triangle. The other corner is (0,1)."""
    arr = build_line_arrangement(_triangle())
    zone = compute_supporting_line_zone(arr, line_index=1)
    verify_supporting_line_zone(arr, zone)
    on_line = {v.point for v in zone.vertices_on_line}
    opposite = {v.point for v in zone.opposite_vertices}
    assert on_line == {Point2D(0, 0), Point2D(1, 0)}
    assert opposite == {Point2D(0, 1)}


def test_supporting_triangle_each_side_has_one_opposite():
    arr = build_line_arrangement(_triangle())
    expected = {
        0: ( {Point2D(0, 0), Point2D(0, 1)}, {Point2D(1, 0)} ),
        1: ( {Point2D(0, 0), Point2D(1, 0)}, {Point2D(0, 1)} ),
        2: ( {Point2D(1, 0), Point2D(0, 1)}, {Point2D(0, 0)} ),
    }
    for line_index, (on_line, opposite) in expected.items():
        zone = compute_supporting_line_zone(arr, line_index)
        verify_supporting_line_zone(arr, zone)
        assert {v.point for v in zone.vertices_on_line} == on_line
        assert {v.point for v in zone.opposite_vertices} == opposite


def test_supporting_line_index_out_of_range():
    arr = build_line_arrangement(_two_axes())
    with pytest.raises(IndexError):
        compute_supporting_line_zone(arr, 2)
