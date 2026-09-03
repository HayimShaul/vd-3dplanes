"""Arrangement invariants from design §21, as exact predicates."""

from __future__ import annotations

from vd3d.arrangement2d.predicates import cross2
from vd3d.arrangement2d.types import Arrangement2D, Face, HalfEdge
from vd3d.geometry.points import Point2D


def twin_involution(arrangement: Arrangement2D) -> bool:
    """``twin(twin(e)) == e`` for every half-edge."""
    for he in arrangement.half_edges:
        twin = arrangement.half_edges[he.twin_id]
        if twin.twin_id != he.id:
            return False
        if twin.edge_id != he.edge_id:
            return False
        if twin.line_index != he.line_index:
            return False
        if twin.direction != (-he.direction[0], -he.direction[1]):
            return False
    return True


def next_prev_consistent(arrangement: Arrangement2D) -> bool:
    for he in arrangement.half_edges:
        if he.next_id is None or he.prev_id is None:
            return False
        nxt = arrangement.half_edges[he.next_id]
        prv = arrangement.half_edges[he.prev_id]
        if nxt.prev_id != he.id or prv.next_id != he.id:
            return False
    return True


def vertex_lies_on_incident_edges(arrangement: Arrangement2D) -> bool:
    for vertex in arrangement.vertices:
        for line_index in vertex.line_indices:
            if not arrangement.lines[line_index].contains(vertex.point):
                return False
        for he_id in vertex.outgoing:
            he = arrangement.half_edges[he_id]
            if he.origin_id != vertex.id:
                return False
            if not arrangement.lines[he.line_index].contains(vertex.point):
                return False
            edge = arrangement.edges[he.edge_id]
            if edge.line_index != he.line_index:
                return False
    return True


def edge_lies_on_supporting_line(arrangement: Arrangement2D) -> bool:
    for edge in arrangement.edges:
        line = arrangement.lines[edge.line_index]
        for vid in (edge.start_vertex, edge.end_vertex):
            if vid is None:
                continue
            if not line.contains(arrangement.vertices[vid].point):
                return False
        for he in arrangement.half_edges:
            if he.edge_id != edge.id:
                continue
            if he.line_index != edge.line_index:
                return False
            if he.origin_id is None:
                continue
            if not line.contains(arrangement.vertices[he.origin_id].point):
                return False
    return True


def face_cycles_valid(arrangement: Arrangement2D) -> bool:
    seen: set[int] = set()
    for face in arrangement.faces:
        cycle = arrangement.cycle(face)
        if face.half_edge_id is None:
            if cycle:
                return False
            continue
        if not cycle:
            return False
        for he in cycle:
            if he.id in seen:
                return False
            seen.add(he.id)
            if he.face_id != face.id:
                return False
        if cycle[0].id != face.half_edge_id:
            # start may not be the recorded half-edge if we re-walk; allow any
            # as long as the recorded one is in the cycle.
            if all(he.id != face.half_edge_id for he in cycle):
                return False
    return seen == {he.id for he in arrangement.half_edges}


def euler_characteristic_plane(arrangement: Arrangement2D) -> bool:
    """``V - E + F = 1`` for an arrangement of lines in the Euclidean plane."""
    v = len(arrangement.vertices)
    e = len(arrangement.edges)
    f = len(arrangement.faces)
    return v - e + f == 1


def _point_on_half_edge_line(arrangement: Arrangement2D, he: HalfEdge) -> Point2D:
    if he.origin_id is not None:
        return arrangement.vertices[he.origin_id].point
    twin = arrangement.half_edges[he.twin_id]
    if twin.origin_id is not None:
        return arrangement.vertices[twin.origin_id].point
    sample, _ = arrangement.lines[he.line_index].sample_points()
    return sample


def side_of_half_edge(arrangement: Arrangement2D, he: HalfEdge, point: Point2D) -> int:
    """Sign of the 2D cross product: +1 left, 0 on the line, -1 right."""
    origin = _point_on_half_edge_line(arrangement, he)
    cr = cross2(he.direction, (point.x - origin.x, point.y - origin.y))
    if cr > 0:
        return 1
    if cr < 0:
        return -1
    return 0


def point_in_face(
    arrangement: Arrangement2D,
    face: Face,
    point: Point2D,
    *,
    closed: bool = False,
) -> bool:
    """True iff ``point`` is on the face side of every boundary half-edge.

    Line-arrangement faces are convex, so this is an exact half-plane test.
    ``closed=True`` includes the boundary.
    """
    cycle = arrangement.cycle(face)
    if not cycle:
        return True
    for he in cycle:
        side = side_of_half_edge(arrangement, he, point)
        if closed:
            if side < 0:
                return False
        elif side <= 0:
            return False
    return True


def representative_in_unique_face(arrangement: Arrangement2D) -> bool:
    for face in arrangement.faces:
        if not point_in_face(arrangement, face, face.representative, closed=False):
            return False
        for other in arrangement.faces:
            if other.id == face.id:
                continue
            if point_in_face(arrangement, other, face.representative, closed=False):
                return False
    return True


def verify_arrangement_invariants(arrangement: Arrangement2D) -> None:
    """Raise ``AssertionError`` if any checked invariant fails."""
    checks = (
        ("twin involution", twin_involution),
        ("next/prev", next_prev_consistent),
        ("vertex on incident edges", vertex_lies_on_incident_edges),
        ("edge on supporting line", edge_lies_on_supporting_line),
        ("face cycles", face_cycles_valid),
        ("Euler V-E+F=1", euler_characteristic_plane),
        ("representative in unique face", representative_in_unique_face),
    )
    for name, pred in checks:
        if not pred(arrangement):
            raise AssertionError(f"arrangement invariant failed: {name}")
