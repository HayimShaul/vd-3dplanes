"""Phase 5 — intersection lines and triple events. Design tests 15–16."""

from __future__ import annotations

import random

import pytest

from tests.oracles.events import (
    brute_force_intersecting_pairs,
    brute_force_triple_points,
    random_planes_general_position,
)
from vd3d.events import (
    Event,
    EventType,
    compute_intersection_lines,
    generate_triple_events,
    verify_intersection_lines,
    verify_triple_events,
)
from vd3d.geometry import Plane, Point3D, line_lies_on_both_planes


def _planes_through_123() -> list[Plane]:
    """Three non-vertical planes meeting at ``(1, 2, 3)`` (design Test 15)."""
    return [
        Plane(id=1, a=1, b=0, c=1, d=-4),  # x + z = 4
        Plane(id=2, a=0, b=1, c=1, d=-5),  # y + z = 5
        Plane(id=3, a=1, b=1, c=1, d=-6),  # x + y + z = 6
    ]


def _parallel_family() -> list[Plane]:
    """Three pairwise-parallel planes (design Test 16)."""
    return [
        Plane(id=1, a=1, b=1, c=1, d=0),
        Plane(id=2, a=1, b=1, c=1, d=-1),
        Plane(id=3, a=1, b=1, c=1, d=-2),
    ]


# ---------------------------------------------------------------------------
# Step 5.1 — intersection lines
# ---------------------------------------------------------------------------


def test_intersection_lines_empty_and_singleton():
    assert compute_intersection_lines([]) == ()
    assert compute_intersection_lines([Plane(id=1, a=1, b=0, c=1, d=0)]) == ()


def test_intersection_lines_two_crossing():
    p = Plane(id=1, a=1, b=0, c=1, d=-4)
    q = Plane(id=2, a=0, b=1, c=1, d=-5)
    lines = compute_intersection_lines([p, q])
    assert len(lines) == 1
    assert lines[0].id == 0
    assert {lines[0].plane_a, lines[0].plane_b} == {1, 2}
    assert line_lies_on_both_planes(lines[0], p, q)
    verify_intersection_lines([p, q], lines)


def test_intersection_lines_skip_parallels():
    planes = [
        Plane(id=1, a=1, b=0, c=1, d=0),
        Plane(id=2, a=1, b=0, c=1, d=-2),
        Plane(id=3, a=0, b=1, c=1, d=0),
    ]
    lines = compute_intersection_lines(planes)
    assert len(lines) == 2
    pairs = {(line.plane_a, line.plane_b) for line in lines}
    assert pairs == {(1, 3), (2, 3)}
    verify_intersection_lines(planes, lines)


def test_intersection_lines_three_general_position():
    planes = _planes_through_123()
    lines = compute_intersection_lines(planes)
    assert len(lines) == 3
    assert [line.id for line in lines] == [0, 1, 2]
    verify_intersection_lines(planes, lines)
    assert brute_force_intersecting_pairs(planes) == tuple(
        tuple(sorted((line.plane_a, line.plane_b))) for line in lines
    )
    meeting = Point3D(1, 2, 3)
    for line in lines:
        assert line.contains(meeting)


def test_intersection_lines_count_formula_random():
    rng = random.Random(15)
    for n in range(2, 7):
        planes = random_planes_general_position(rng, n)
        lines = compute_intersection_lines(planes)
        assert len(lines) == n * (n - 1) // 2
        verify_intersection_lines(planes, lines)
        expected_pairs = brute_force_intersecting_pairs(planes)
        got_pairs = tuple(tuple(sorted((line.plane_a, line.plane_b))) for line in lines)
        assert got_pairs == expected_pairs


def test_intersection_lines_reject_duplicate_ids():
    p = Plane(id=1, a=1, b=0, c=1, d=0)
    q = Plane(id=1, a=0, b=1, c=1, d=0)
    with pytest.raises(ValueError, match="unique"):
        compute_intersection_lines([p, q])


# ---------------------------------------------------------------------------
# Step 5.2 — triple events
# ---------------------------------------------------------------------------


def test_triple_event_known_point_z3():
    """Design Test 15: planes through (1, 2, 3) emit one event at z=3."""
    planes = _planes_through_123()
    events = generate_triple_events(planes)
    verify_triple_events(planes, events)
    assert len(events) == 1
    event = events[0]
    assert event.z == 3
    assert event.type is EventType.TRIPLE_INTERSECTION
    assert event.geometric_data == Point3D(1, 2, 3)
    assert event.plane_ids == (1, 2, 3)
    assert event.stable_id == "triple:1:2:3"
    assert event.sort_key == (event.z, int(EventType.TRIPLE_INTERSECTION), "triple:1:2:3")


def test_parallel_family_has_zero_triples():
    """Design Test 16: a parallel family produces no triple events."""
    planes = _parallel_family()
    assert compute_intersection_lines(planes) == ()
    events = generate_triple_events(planes)
    assert events == ()
    verify_triple_events(planes, events)


def test_parallel_pencil_has_zero_triples():
    """Three planes through one line: three pairwise lines, no unique point."""
    planes = [
        Plane(id=1, a=1, b=0, c=0, d=0),
        Plane(id=2, a=0, b=1, c=0, d=0),
        Plane(id=3, a=1, b=1, c=0, d=0),
    ]
    lines = compute_intersection_lines(planes)
    assert len(lines) == 3
    assert generate_triple_events(planes) == ()
    verify_intersection_lines(planes, lines)
    verify_triple_events(planes, generate_triple_events(planes))


def test_two_planes_have_no_triple():
    planes = _planes_through_123()[:2]
    assert generate_triple_events(planes) == ()


def test_every_triple_point_lies_on_its_planes_random():
    rng = random.Random(16)
    for n in range(3, 7):
        planes = random_planes_general_position(rng, n)
        events = generate_triple_events(planes)
        assert len(events) == n * (n - 1) * (n - 2) // 6
        verify_triple_events(planes, events)
        oracle = brute_force_triple_points(planes)
        assert len(events) == len(oracle)
        for event, (point, plane_ids) in zip(events, oracle):
            assert event.geometric_data == point
            assert event.plane_ids == plane_ids
            assert event.z == point.z


def test_triple_events_sorted_by_z():
    planes = [
        Plane(id=1, a=1, b=0, c=0, d=-1),  # x = 1
        Plane(id=2, a=0, b=1, c=0, d=-2),  # y = 2
        Plane(id=3, a=0, b=0, c=1, d=-3),  # z = 3  → (1,2,3)
        Plane(id=4, a=0, b=0, c=1, d=-5),  # z = 5  → (1,2,5) with 1,2
    ]
    # Triples: (1,2,3) at z=3; (1,2,4) at z=5; (1,3,4) singular (two z=const);
    # (2,3,4) singular.
    events = generate_triple_events(planes)
    verify_triple_events(planes, events)
    assert [event.z for event in events] == [3, 5]
    assert events[0].sort_key < events[1].sort_key


def test_event_sort_key_orders_type_then_id():
    low = Event(
        z=1,
        type=EventType.TRIPLE_INTERSECTION,
        geometric_data=Point3D(0, 0, 1),
        plane_ids=(1, 2, 3),
    )
    high = Event(
        z=1,
        type=EventType.VERTICAL_ALIGNMENT,
        geometric_data=Point3D(0, 0, 1),
        plane_ids=(1, 2),
    )
    assert low.sort_key < high.sort_key


def test_event_rejects_float_z():
    with pytest.raises(TypeError):
        Event(
            z=1.0,  # type: ignore[arg-type]
            type=EventType.TRIPLE_INTERSECTION,
            geometric_data=Point3D(0, 0, 1),
            plane_ids=(1, 2, 3),
        )


def test_generate_triple_events_reject_duplicate_ids():
    planes = [
        Plane(id=7, a=1, b=0, c=1, d=0),
        Plane(id=7, a=0, b=1, c=1, d=0),
        Plane(id=8, a=1, b=1, c=1, d=-1),
    ]
    with pytest.raises(ValueError, match="unique"):
        generate_triple_events(planes)
