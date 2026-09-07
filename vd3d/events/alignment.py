"""Vertical-alignment events from the wall zone (design §8.2–8.5)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from vd3d.arrangement2d import build_line_arrangement, intersect_lines_2d
from vd3d.arrangement2d.intersect import COINCIDENT
from vd3d.arrangement2d.types import Arrangement2D
from vd3d.events.ids import require_unique_plane_ids
from vd3d.events.lines import compute_intersection_lines
from vd3d.events.types import AlignmentGeometry, Event, EventType
from vd3d.events.wall import (
    VerticalIntersectionLine,
    alignment_z,
    build_vertical_wall,
    lift_wall_point,
    line_xz_param,
    lines_meet,
    point_on_line_at_z,
)
from vd3d.geometry.intersections import slice_plane_at_z
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane, normals_parallel
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar, as_scalar
from vd3d.zone import compute_supporting_line_zone
from vd3d.zone.types import SupportingLineZone


def represent_line_on_wall(line: Line3D) -> Line2D:
    """``L`` in the wall frame ``(z, y)``: ``y - ay z - by = 0``."""
    _ax, _bx, ay, by = line_xz_param(line)
    return Line2D(a=-ay, b=1, c=-by, source_plane_id=None, id=0)


def plane_on_wall(plane: Plane, wall: Plane) -> Line2D | None:
    """Trace of ``plane`` on the y-parallel wall, in the ``(z, y)`` frame.

    Substitutes ``x = -(c_w z + d_w) / a_w`` into the plane. ``None`` if
    the result is degenerate (plane parallel to the wall).
    """
    if wall.a == 0 or wall.b != 0:
        raise ValueError("wall must be y-parallel with a ≠ 0")
    ax = -wall.c / wall.a
    bx = -wall.d / wall.a
    a_z = plane.a * ax + plane.c
    a_y = plane.b
    a_c = plane.a * bx + plane.d
    if a_y == 0 and a_z == 0:
        return None
    return Line2D(a=a_z, b=a_y, c=a_c, source_plane_id=plane.id)


def slice_planes_by_wall(planes: Sequence[Plane], wall: Plane, line: Line3D) -> tuple[Line2D, ...]:
    """Wall-frame lines: ``L`` first, then every other non-parallel plane.

    Source planes of ``L`` are omitted; they are coincident with ``L``.
    Other traces coincident with an already-included wall-line are
    dropped (four planes through one point can produce duplicate traces).
    """
    query = represent_line_on_wall(line)
    source = {line.plane_a, line.plane_b}
    out: list[Line2D] = [query]
    next_id = 1
    for plane in planes:
        if plane.id in source:
            continue
        if normals_parallel(plane, wall):
            continue
        traced = plane_on_wall(plane, wall)
        if traced is None:
            continue
        skip = False
        for existing in out:
            if intersect_lines_2d(existing, traced) is COINCIDENT:
                skip = True
                break
        if skip:
            continue
        out.append(Line2D(
            a=traced.a,
            b=traced.b,
            c=traced.c,
            source_plane_id=traced.source_plane_id,
            id=next_id,
        ))
        next_id += 1
    return tuple(out)


@dataclass(frozen=True, slots=True)
class WallLineZone:
    """Supporting-line zone of ``L`` in the arrangement on its wall."""

    line: Line3D
    wall: Plane
    arrangement: Arrangement2D
    query_index: int
    supporting: SupportingLineZone


def compute_line_zone_on_wall(
    line: Line3D, wall: Plane, planes: Sequence[Plane]
) -> WallLineZone:
    """``COMPUTE_LINE_ZONE_ON_WALL``: arrangement on ``W``, zone of ``L``."""
    wall_lines = slice_planes_by_wall(planes, wall, line)
    arrangement = build_line_arrangement(wall_lines)
    query_index = 0
    supporting = compute_supporting_line_zone(arrangement, query_index)
    return WallLineZone(
        line=line,
        wall=wall,
        arrangement=arrangement,
        query_index=query_index,
        supporting=supporting,
    )


def alignment_visible(
    left: Line3D,
    right: Line3D,
    z: int | Scalar | str,
    planes: Sequence[Plane],
    lines: Sequence[Line3D] | None = None,
) -> bool:
    """Open vertical segment between the two slice-vertices hits no feature."""
    z = as_scalar(z)
    p = point_on_line_at_z(left, z)
    q = point_on_line_at_z(right, z)
    if p.x != q.x or p.y == q.y:
        return False
    y_lo, y_hi = (p.y, q.y) if p.y < q.y else (q.y, p.y)
    x = p.x
    others = lines if lines is not None else compute_intersection_lines(list(planes))
    for other in others:
        if other.id is not None and other.id in {left.id, right.id}:
            continue
        if other is left or other is right:
            continue
        try:
            mid = point_on_line_at_z(other, z)
        except VerticalIntersectionLine:
            continue
        if mid.x == x and y_lo < mid.y < y_hi:
            return False
    for plane in planes:
        sliced = slice_plane_at_z(plane, z)
        if sliced is None:
            continue
        if sliced.b != 0:
            y_hit = -(sliced.a * x + sliced.c) / sliced.b
            if y_lo < y_hit < y_hi:
                return False
            continue
        if sliced.contains(Point2D(x, (y_lo + y_hi) / 2)):
            return False
    return True


def make_alignment_event(left: Line3D, right: Line3D, z: Scalar) -> Event:
    if left.id is None or right.id is None:
        raise ValueError("alignment lines must have ids")
    ids = tuple(sorted((left.id, right.id)))
    plane_ids = tuple(
        sorted({left.plane_a, left.plane_b, right.plane_a, right.plane_b} - {None})
    )
    first, second = (left, right) if left.id <= right.id else (right, left)
    p = point_on_line_at_z(first, z)
    q = point_on_line_at_z(second, z)
    return Event(
        z=z,
        type=EventType.VERTICAL_ALIGNMENT,
        geometric_data=AlignmentGeometry(point_a=p, point_b=q),
        plane_ids=plane_ids,
        line_ids=ids,
        stable_id=f"align:{ids[0]}:{ids[1]}",
    )


def validate_alignment_event(
    event: Event,
    planes: Sequence[Plane],
    lines: Sequence[Line3D],
) -> bool:
    """``VALIDATE_ALIGNMENT_EVENT``: same-x, not a triple, visible."""
    if event.type is not EventType.VERTICAL_ALIGNMENT:
        return False
    if len(event.line_ids) != 2:
        return False
    by_id = {line.id: line for line in lines}
    try:
        left = by_id[event.line_ids[0]]
        right = by_id[event.line_ids[1]]
    except KeyError:
        return False
    try:
        z = alignment_z(left, right)
    except VerticalIntersectionLine:
        return False
    if z is None or z != event.z:
        return False
    if lines_meet(left, right):
        return False
    p = point_on_line_at_z(left, z)
    q = point_on_line_at_z(right, z)
    if p.x != q.x:
        return False
    return alignment_visible(left, right, z, planes, lines)


def alignment_event_for_opposite_vertex(
    zone: WallLineZone,
    vertex,
    planes: Sequence[Plane],
    lines: Sequence[Line3D],
) -> Event | None:
    """Event if this opposite vertex is a visible alignment, else ``None``."""
    lifted = lift_wall_point(vertex.point, zone.wall)
    z = lifted.z
    for other in lines:
        if other.id == zone.line.id:
            continue
        try:
            at_z = point_on_line_at_z(other, z)
        except VerticalIntersectionLine:
            continue
        if at_z != lifted:
            continue
        event = make_alignment_event(zone.line, other, z)
        if validate_alignment_event(event, planes, lines):
            return event
        return None
    return None


def extract_alignment_events_from_zone(
    zone: WallLineZone,
    planes: Sequence[Plane],
    lines: Sequence[Line3D],
) -> tuple[Event, ...]:
    """Opposite wall-vertices that validate as alignments."""
    events: list[Event] = []
    for vertex in zone.supporting.opposite_vertices:
        event = alignment_event_for_opposite_vertex(zone, vertex, planes, lines)
        if event is not None:
            events.append(event)
    return tuple(events)


def canonical_event_key(event: Event) -> tuple:
    if event.type is EventType.TRIPLE_INTERSECTION:
        return (int(event.type), event.z, event.plane_ids)
    return (int(event.type), event.z, tuple(sorted(event.line_ids)))


def deduplicate_events(events: Sequence[Event]) -> tuple[Event, ...]:
    """Keep one event per canonical key (same alignment is found twice)."""
    seen: dict[tuple, Event] = {}
    for event in events:
        seen[canonical_event_key(event)] = event
    return tuple(seen.values())


def enumerate_alignment_pairs(planes: Sequence[Plane]) -> tuple[Event, ...]:
    """All-pairs alignment generator. Definition of correctness (the oracle)."""
    require_unique_plane_ids(planes)
    lines = compute_intersection_lines(planes)
    events: list[Event] = []
    for i, left in enumerate(lines):
        for right in lines[i + 1 :]:
            try:
                z = alignment_z(left, right)
            except VerticalIntersectionLine:
                continue
            if z is None:
                continue
            if lines_meet(left, right):
                continue
            if not alignment_visible(left, right, z, planes, lines):
                continue
            events.append(make_alignment_event(left, right, z))
    return tuple(sorted(deduplicate_events(events), key=lambda event: event.sort_key))


def generate_alignment_events(planes: Sequence[Plane]) -> tuple[Event, ...]:
    """``GENERATE_ALIGNMENT_EVENTS`` via the wall-zone of every intersection line."""
    require_unique_plane_ids(planes)
    lines = compute_intersection_lines(planes)
    events: list[Event] = []
    for line in lines:
        try:
            wall = build_vertical_wall(line)
        except VerticalIntersectionLine:
            continue
        zone = compute_line_zone_on_wall(line, wall, planes)
        events.extend(extract_alignment_events_from_zone(zone, planes, lines))
    return tuple(sorted(deduplicate_events(events), key=lambda event: event.sort_key))
