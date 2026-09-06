"""Event-generation invariants checked by unit tests."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events.types import AlignmentGeometry, Event, EventType
from vd3d.geometry.invariants import line_lies_on_both_planes, point_lies_on_planes
from vd3d.geometry.points import Point3D
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane, normals_parallel


def planes_by_id(planes: Sequence[Plane]) -> dict[int, Plane]:
    return {plane.id: plane for plane in planes}


def intersection_line_count_matches_pairs(planes: Sequence[Plane], lines: Sequence[Line3D]) -> bool:
    """``|lines|`` equals the number of non-parallel pairs."""
    expected = 0
    for i, p in enumerate(planes):
        for q in planes[i + 1 :]:
            if not normals_parallel(p, q):
                expected += 1
    return len(lines) == expected


def each_line_lies_on_source_planes(planes: Sequence[Plane], lines: Sequence[Line3D]) -> bool:
    lookup = planes_by_id(planes)
    seen_ids: set[int] = set()
    for line in lines:
        if line.id is None or line.id in seen_ids:
            return False
        seen_ids.add(line.id)
        if line.plane_a is None or line.plane_b is None:
            return False
        if line.plane_a not in lookup or line.plane_b not in lookup:
            return False
        if not line_lies_on_both_planes(line, lookup[line.plane_a], lookup[line.plane_b]):
            return False
    return True


def triple_event_point_on_planes(event: Event, planes: Sequence[Plane]) -> bool:
    if event.type is not EventType.TRIPLE_INTERSECTION:
        return False
    point = event.geometric_data
    if not isinstance(point, Point3D):
        return False
    if point.z != event.z:
        return False
    lookup = planes_by_id(planes)
    try:
        trio = tuple(lookup[i] for i in event.plane_ids)
    except KeyError:
        return False
    if len(trio) != 3:
        return False
    return point_lies_on_planes(point, *trio)


def events_sorted_by_key(events: Sequence[Event]) -> bool:
    keys = [event.sort_key for event in events]
    return keys == sorted(keys) and len(keys) == len(set(keys))


def verify_intersection_lines(planes: Sequence[Plane], lines: Sequence[Line3D]) -> None:
    checks = (
        (
            "line count matches non-parallel pairs",
            lambda: intersection_line_count_matches_pairs(planes, lines),
        ),
        (
            "each line lies on both source planes",
            lambda: each_line_lies_on_source_planes(planes, lines),
        ),
    )
    for name, pred in checks:
        if not pred():
            raise AssertionError(f"intersection-line invariant failed: {name}")


def alignment_event_same_x(event: Event) -> bool:
    if event.type is not EventType.VERTICAL_ALIGNMENT:
        return False
    data = event.geometric_data
    if not isinstance(data, AlignmentGeometry):
        return False
    return (
        data.point_a.z == event.z
        and data.point_b.z == event.z
        and data.point_a.x == data.point_b.x
        and data.point_a.y != data.point_b.y
    )


def verify_triple_events(planes: Sequence[Plane], events: Sequence[Event]) -> None:
    checks = (
        ("events sorted by (z, type, stable_id)", lambda: events_sorted_by_key(events)),
        (
            "each triple point lies on its three planes",
            lambda: all(triple_event_point_on_planes(event, planes) for event in events),
        ),
    )
    for name, pred in checks:
        if not pred():
            raise AssertionError(f"triple-event invariant failed: {name}")


def verify_alignment_events(events: Sequence[Event]) -> None:
    checks = (
        ("events sorted by (z, type, stable_id)", lambda: events_sorted_by_key(events)),
        (
            "each alignment has the same x and distinct y",
            lambda: all(alignment_event_same_x(event) for event in events),
        ),
    )
    for name, pred in checks:
        if not pred():
            raise AssertionError(f"alignment-event invariant failed: {name}")


def events_have_unique_z(events: Sequence[Event]) -> bool:
    zs = [event.z for event in events]
    return len(zs) == len(set(zs))


def z_strictly_below_all_events(z, events: Sequence[Event]) -> bool:
    from vd3d.geometry.scalar import as_scalar

    height = as_scalar(z)
    return all(height < event.z for event in events)


def verify_all_events(events: Sequence[Event]) -> None:
    triples = [event for event in events if event.type is EventType.TRIPLE_INTERSECTION]
    alignments = [event for event in events if event.type is EventType.VERTICAL_ALIGNMENT]
    checks = (
        ("events sorted by unique (z, type, stable_id)", lambda: events_sorted_by_key(events)),
        (
            "every event is a triple or an alignment",
            lambda: len(triples) + len(alignments) == len(events),
        ),
        (
            "each alignment has the same x and distinct y",
            lambda: all(alignment_event_same_x(event) for event in alignments),
        ),
    )
    for name, pred in checks:
        if not pred():
            raise AssertionError(f"all-event invariant failed: {name}")
