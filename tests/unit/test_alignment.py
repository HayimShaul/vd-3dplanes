"""Phase 6 — alignment events. Design tests 17–18."""

from __future__ import annotations

import random

import pytest

from tests.oracles.alignment import (
    alignment_keys,
    brute_force_alignment_events,
    planes_alignment_at_z2,
    planes_alignment_blocked,
    planes_never_align,
)
from tests.oracles.events import random_planes_general_position
from vd3d.events import (
    AlignmentGeometry,
    EventType,
    VerticalIntersectionLine,
    alignment_visible,
    alignment_z,
    build_vertical_wall,
    compute_intersection_lines,
    compute_line_zone_on_wall,
    generate_alignment_events,
    generate_triple_events,
    lines_meet,
    point_on_line_at_z,
    represent_line_on_wall,
    slice_planes_by_wall,
    validate_alignment_event,
    verify_alignment_events,
    wall_frame_point,
    x_of_line,
    y_of_line,
)
from vd3d.geometry import Line3D, Plane, Point3D


def _line_l12() -> Line3D:
    """``(t, 0, t)``: ``x = z``, ``y = 0``."""
    return Line3D(point=Point3D(0, 0, 0), direction=(1, 0, 1), plane_a=1, plane_b=2, id=0)


def _line_l34() -> Line3D:
    """``(4 - s, 3, s)``: ``x = 4 - z``, ``y = 3``."""
    return Line3D(point=Point3D(4, 3, 0), direction=(-1, 0, 1), plane_a=3, plane_b=4, id=1)


# ---------------------------------------------------------------------------
# Step 6.1 — x_L(z) and the y-parallel wall
# ---------------------------------------------------------------------------


def test_x_y_of_line_parametrization():
    line = _line_l12()
    assert x_of_line(line, 0) == 0
    assert x_of_line(line, 2) == 2
    assert x_of_line(line, -1) == -1
    assert y_of_line(line, 5) == 0
    assert point_on_line_at_z(line, 2) == Point3D(2, 0, 2)
    assert line.contains(point_on_line_at_z(line, 2))


def test_x_of_line_rejects_vertical():
    vertical = Line3D(point=Point3D(1, 2, 3), direction=(1, 1, 0))
    with pytest.raises(VerticalIntersectionLine):
        x_of_line(vertical, 0)
    with pytest.raises(VerticalIntersectionLine):
        build_vertical_wall(vertical)


def test_vertical_wall_contains_line_and_has_no_y_term():
    line = _line_l12()
    wall = build_vertical_wall(line)
    assert wall.b == 0
    assert wall.a != 0
    for t in (-2, 0, 1, "5/2"):
        assert wall.contains(line.point_at(t))
    # x = z  →  x - z = 0
    assert wall.a == 1
    assert wall.c == -1
    assert wall.d == 0


def test_vertical_wall_of_l34():
    wall = build_vertical_wall(_line_l34())
    assert wall.b == 0
    for t in (-1, 0, 2):
        assert wall.contains(_line_l34().point_at(t))
    # x = 4 - z  →  x + z - 4 = 0
    assert point_on_line_at_z(_line_l34(), 2) == Point3D(2, 3, 2)
    assert wall.contains(Point3D(2, 99, 2))


# ---------------------------------------------------------------------------
# Step 6.2 — brute-force oracle (definition of correctness)
# ---------------------------------------------------------------------------


def test_alignment_z_known_pair():
    assert alignment_z(_line_l12(), _line_l34()) == 2
    assert not lines_meet(_line_l12(), _line_l34())


def test_alignment_z_parallel_walls_never():
    other = Line3D(point=Point3D(1, 3, 0), direction=(1, 0, 1), id=2)
    assert alignment_z(_line_l12(), other) is None


def _named_line(planes: list[Plane], plane_ids: set[int]):
    for line in compute_intersection_lines(planes):
        if {line.plane_a, line.plane_b} == plane_ids:
            return line
    raise AssertionError(f"missing line for planes {plane_ids}")


def test_known_alignment_at_z2():
    """Design Test 17: one visible alignment at ``z = 2``."""
    planes = planes_alignment_at_z2()
    l12 = _named_line(planes, {1, 2})
    l34 = _named_line(planes, {3, 4})
    assert alignment_z(l12, l34) == 2
    assert alignment_visible(l12, l34, 2, planes)
    events = brute_force_alignment_events(planes)
    verify_alignment_events(events)
    assert len(events) == 1
    event = events[0]
    assert event.z == 2
    assert event.type is EventType.VERTICAL_ALIGNMENT
    assert isinstance(event.geometric_data, AlignmentGeometry)
    data = event.geometric_data
    assert data.point_a.x == data.point_b.x == 2
    assert {data.point_a.y, data.point_b.y} == {0, 3}
    assert tuple(sorted(event.line_ids)) == tuple(sorted((l12.id, l34.id)))


def test_never_align_is_empty():
    events = brute_force_alignment_events(planes_never_align())
    assert events == ()


def test_visibility_rejects_vertex_between():
    planes = planes_alignment_blocked()
    lines = compute_intersection_lines(planes)
    l12 = next(line for line in lines if {line.plane_a, line.plane_b} == {1, 2})
    l34 = next(line for line in lines if {line.plane_a, line.plane_b} == {3, 4})
    assert alignment_z(l12, l34) == 2
    assert not alignment_visible(l12, l34, 2, planes, lines)
    keys = alignment_keys(brute_force_alignment_events(planes))
    assert (int(EventType.VERTICAL_ALIGNMENT), 2, tuple(sorted((l12.id, l34.id)))) not in keys


def test_three_planes_have_no_alignments():
    planes = [
        Plane(id=1, a=1, b=0, c=1, d=-4),
        Plane(id=2, a=0, b=1, c=1, d=-5),
        Plane(id=3, a=1, b=1, c=1, d=-6),
    ]
    assert brute_force_alignment_events(planes) == ()
    assert generate_alignment_events(planes) == ()


# ---------------------------------------------------------------------------
# Step 6.3 — slice by wall
# ---------------------------------------------------------------------------


def test_slice_wall_has_l_and_other_planes():
    planes = planes_alignment_at_z2()
    lines = compute_intersection_lines(planes)
    l12 = next(line for line in lines if {line.plane_a, line.plane_b} == {1, 2})
    wall = build_vertical_wall(l12)
    wall_lines = slice_planes_by_wall(planes, wall, l12)
    assert wall_lines[0] == represent_line_on_wall(l12)
    sources = {line.source_plane_id for line in wall_lines[1:]}
    assert sources == {3, 4}
    assert len(wall_lines) == 3


def test_wall_frame_alignment_is_vertical():
    """Same ``z`` is the same wall-x, so alignment is a vertical segment."""
    z = alignment_z(_line_l12(), _line_l34())
    assert z == 2
    on_l = wall_frame_point(_line_l12(), z)
    on_other = wall_frame_point(_line_l34(), z)
    assert on_l.x == on_other.x == z
    assert on_l.y == 0
    assert on_other.y == 3


def test_wall_arrangement_l_is_one_line():
    planes = planes_alignment_at_z2()
    lines = compute_intersection_lines(planes)
    l12 = next(line for line in lines if {line.plane_a, line.plane_b} == {1, 2})
    wall = build_vertical_wall(l12)
    zone = compute_line_zone_on_wall(l12, wall, planes)
    assert zone.query_index == 0
    assert len(zone.arrangement.lines) == 3
    query = zone.arrangement.lines[0]
    for t in (0, 1, -2):
        p = point_on_line_at_z(l12, t)
        assert query.contains(wall_frame_point(l12, t))
        assert wall.contains(p)


# ---------------------------------------------------------------------------
# Step 6.4 — zone extraction vs oracle
# ---------------------------------------------------------------------------


def test_zone_matches_oracle_test17():
    """Design Test 18 on the hand fixture."""
    planes = planes_alignment_at_z2()
    zone_events = generate_alignment_events(planes)
    oracle = brute_force_alignment_events(planes)
    verify_alignment_events(zone_events)
    assert alignment_keys(zone_events) == alignment_keys(oracle)
    assert len(zone_events) == 1
    assert zone_events[0].z == 2
    lines = compute_intersection_lines(planes)
    assert validate_alignment_event(zone_events[0], planes, lines)


def test_zone_matches_oracle_never_and_blocked():
    for planes in (planes_never_align(), planes_alignment_blocked()):
        assert alignment_keys(generate_alignment_events(planes)) == alignment_keys(
            brute_force_alignment_events(planes)
        )


def test_zone_matches_oracle_random():
    rng = random.Random(18)
    for n in (4, 5, 6):
        for _ in range(4):
            planes = random_planes_general_position(rng, n)
            zone_events = generate_alignment_events(planes)
            oracle = brute_force_alignment_events(planes)
            verify_alignment_events(zone_events)
            assert alignment_keys(zone_events) == alignment_keys(oracle)
            for event in zone_events:
                assert validate_alignment_event(
                    event, planes, compute_intersection_lines(planes)
                )


def test_dedup_from_both_walls():
    planes = planes_alignment_at_z2()
    events = generate_alignment_events(planes)
    assert len(events) == 1
    assert events[0].stable_id.startswith("align:")


def test_triples_unchanged_on_alignment_fixture():
    planes = planes_alignment_at_z2()
    triples = generate_triple_events(planes)
    # Four planes in GP: C(4,3)=4 triples.
    assert len(triples) == 4
