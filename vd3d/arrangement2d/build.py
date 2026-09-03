"""Build a DCEL for an arrangement of unbounded 2D lines.

Orientation (locked):

- Outgoing half-edges around a vertex are ordered **counterclockwise**,
  starting from the +x axis (see ``predicates.quadrant``).
- Each half-edge has its incident face on the **left**.
- Bounded-face cycles therefore walk counterclockwise.

No bounding box is part of the combinatorics. Infinity is handled by
linking half-edges whose destination is not a finite vertex.
"""

from __future__ import annotations

from functools import cmp_to_key

from vd3d.arrangement2d.pieces import compute_edge_pieces, compute_vertices
from vd3d.arrangement2d.predicates import (
    Direction,
    cmp_direction_ccw,
    left_normal,
    line_direction,
)
from vd3d.arrangement2d.types import Arrangement2D, EdgePiece, Face, HalfEdge, Vertex
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar


def build_line_arrangement(lines: list[Line2D] | tuple[Line2D, ...]) -> Arrangement2D:
    """``BUILD_LINE_ARRANGEMENT``: vertices, pieces, half-edges, faces."""
    packed = tuple(lines)
    vertices = compute_vertices(packed)
    edges = compute_edge_pieces(packed, vertices)
    half_edges, vertices = _build_half_edges(packed, vertices, edges)
    faces = _trace_faces(packed, vertices, edges, half_edges)
    return Arrangement2D(
        lines=packed,
        vertices=vertices,
        edges=edges,
        half_edges=tuple(half_edges),
        faces=faces,
    )


def _build_half_edges(
    lines: tuple[Line2D, ...],
    vertices: tuple[Vertex, ...],
    edges: tuple[EdgePiece, ...],
) -> tuple[list[HalfEdge], tuple[Vertex, ...]]:
    half_edges: list[HalfEdge] = []
    for edge in edges:
        direction = line_direction(lines[edge.line_index])
        opposite = (-direction[0], -direction[1])
        # Forward half-edge travels in ``direction`` (increasing t).
        forward_origin = edge.start_vertex
        backward_origin = edge.end_vertex
        fwd_id = 2 * edge.id
        bwd_id = 2 * edge.id + 1
        half_edges.append(
            HalfEdge(
                id=fwd_id,
                origin_id=forward_origin,
                twin_id=bwd_id,
                edge_id=edge.id,
                line_index=edge.line_index,
                direction=direction,
            )
        )
        half_edges.append(
            HalfEdge(
                id=bwd_id,
                origin_id=backward_origin,
                twin_id=fwd_id,
                edge_id=edge.id,
                line_index=edge.line_index,
                direction=opposite,
            )
        )

    if any(he.id != i for i, he in enumerate(half_edges)):
        raise RuntimeError("half-edge ids must equal their indices")

    mutable = list(vertices)
    _link_at_vertices(mutable, half_edges)
    _link_at_infinity(lines, mutable, half_edges)
    return half_edges, tuple(mutable)


def _link_at_vertices(vertices: list[Vertex], half_edges: list[HalfEdge]) -> None:
    """At each finite vertex: next of the incoming twin is the clockwise outgoing."""
    by_origin: dict[int, list[HalfEdge]] = {}
    for he in half_edges:
        if he.origin_id is None:
            continue
        by_origin.setdefault(he.origin_id, []).append(he)

    for vertex in vertices:
        outgoing = by_origin.get(vertex.id, [])
        outgoing.sort(key=cmp_to_key(lambda a, b: cmp_direction_ccw(a.direction, b.direction)))
        vertex.outgoing = tuple(he.id for he in outgoing)
        n = len(outgoing)
        if n == 0:
            continue
        for i, he in enumerate(outgoing):
            # Face on the left ⇒ arriving via twin(he) turns onto the
            # previous (clockwise) outgoing half-edge.
            pred = outgoing[(i - 1) % n]
            twin = half_edges[he.twin_id]
            twin.next_id = pred.id
            pred.prev_id = twin.id


def _link_at_infinity(
    lines: tuple[Line2D, ...],
    vertices: list[Vertex],
    half_edges: list[HalfEdge],
) -> None:
    """Connect half-edges that depart to infinity, in CCW direction order.

    Parallel rays share a direction; they are ordered by increasing offset
    along the left normal of that direction.
    """
    departing = [he for he in half_edges if he.next_id is None]
    if not departing:
        return

    def cmp(a: HalfEdge, b: HalfEdge) -> int:
        keyed = cmp_direction_ccw(a.direction, b.direction)
        if keyed != 0:
            return keyed
        oa = _left_offset(a, lines, vertices)
        ob = _left_offset(b, lines, vertices)
        if oa < ob:
            return -1
        if oa > ob:
            return 1
        return a.id - b.id

    departing.sort(key=cmp_to_key(cmp))
    n = len(departing)
    for i, he in enumerate(departing):
        nxt = departing[(i + 1) % n]
        twin_next = half_edges[nxt.twin_id]
        he.next_id = twin_next.id
        twin_next.prev_id = he.id


def _left_offset(
    he: HalfEdge,
    lines: tuple[Line2D, ...],
    vertices: list[Vertex],
) -> Scalar:
    lx, ly = left_normal(he.direction)
    point = _point_on_supporting_line(he, lines, vertices)
    return lx * point.x + ly * point.y


def _point_on_supporting_line(
    he: HalfEdge,
    lines: tuple[Line2D, ...],
    vertices: list[Vertex],
) -> Point2D:
    if he.origin_id is not None:
        return vertices[he.origin_id].point
    sample, _ = lines[he.line_index].sample_points()
    return sample


def _trace_faces(
    lines: tuple[Line2D, ...],
    vertices: tuple[Vertex, ...],
    edges: tuple[EdgePiece, ...],
    half_edges: list[HalfEdge],
) -> tuple[Face, ...]:
    if not half_edges:
        return (
            Face(id=0, half_edge_id=None, unbounded=True, representative=Point2D(0, 0)),
        )

    visited: set[int] = set()
    faces: list[Face] = []
    for start in half_edges:
        if start.id in visited:
            continue
        cycle = _walk_cycle(start, half_edges)
        for he in cycle:
            visited.add(he.id)
        unbounded = any(edges[he.edge_id].kind != "segment" for he in cycle)
        face = Face(
            id=len(faces),
            half_edge_id=start.id,
            unbounded=unbounded,
            representative=Point2D(0, 0),
        )
        for he in cycle:
            he.face_id = face.id
        face.representative = _representative_point(
            cycle, lines, vertices, edges, half_edges
        )
        faces.append(face)
    return tuple(faces)


def _walk_cycle(start: HalfEdge, half_edges: list[HalfEdge]) -> list[HalfEdge]:
    if start.next_id is None:
        raise RuntimeError(f"half-edge {start.id} has no next")
    cycle = [start]
    current = half_edges[start.next_id]
    guard = len(half_edges) + 1
    while current.id != start.id:
        cycle.append(current)
        if current.next_id is None:
            raise RuntimeError(f"half-edge {current.id} has no next")
        current = half_edges[current.next_id]
        if len(cycle) > guard:
            raise RuntimeError("half-edge next pointers did not form a cycle")
    return cycle


def _representative_point(
    cycle: list[HalfEdge],
    lines: tuple[Line2D, ...],
    vertices: tuple[Vertex, ...],
    edges: tuple[EdgePiece, ...],
    half_edges: list[HalfEdge],
) -> Point2D:
    """An exact interior sample: shoot from a boundary point along the left normal.

    Faces of a line arrangement are convex, so the first positive hit against
    another cycle line, halved, is inside. No hit means the face is unbounded
    in that direction; step by 1.
    """
    seed, left = _interior_seed(cycle, lines, vertices, edges, half_edges)
    best_t: Scalar | None = None
    for he in cycle:
        line = lines[he.line_index]
        denom = line.a * left[0] + line.b * left[1]
        if denom == 0:
            continue
        t = -line.eval(seed) / denom
        if t <= 0:
            continue
        if best_t is None or t < best_t:
            best_t = t
    step = 1 if best_t is None else best_t / 2
    return Point2D(seed.x + step * left[0], seed.y + step * left[1])


def _interior_seed(
    cycle: list[HalfEdge],
    lines: tuple[Line2D, ...],
    vertices: tuple[Vertex, ...],
    edges: tuple[EdgePiece, ...],
    half_edges: list[HalfEdge],
) -> tuple[Point2D, Direction]:
    """A point on the boundary together with the left-normal of that half-edge."""
    for he in cycle:
        kind = edges[he.edge_id].kind
        dx, dy = he.direction
        left = left_normal(he.direction)
        if kind == "segment" and he.origin_id is not None:
            twin_origin = half_edges[he.twin_id].origin_id
            if twin_origin is None:
                raise RuntimeError("segment half-edge twin missing origin")
            a = vertices[he.origin_id].point
            b = vertices[twin_origin].point
            return Point2D((a.x + b.x) / 2, (a.y + b.y) / 2), left
        if he.origin_id is not None:
            origin = vertices[he.origin_id].point
            return Point2D(origin.x + dx, origin.y + dy), left
        twin = half_edges[he.twin_id]
        if twin.origin_id is not None:
            dest = vertices[twin.origin_id].point
            return Point2D(dest.x - dx, dest.y - dy), left
        sample, _ = lines[he.line_index].sample_points()
        return sample, left
    raise RuntimeError("empty face cycle")
