"""``FACE_CONTAINING_POINT``, ``FACE_ON_OTHER_SIDE``, and ``COMPUTE_ZONE``."""

from __future__ import annotations

from vd3d.arrangement2d.invariants import point_in_face
from vd3d.arrangement2d.predicates import parameter_on_line
from vd3d.arrangement2d.types import Arrangement2D, EdgePiece, Face
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.zone.crossings import (
    coincident_line_indices,
    collapse_crossings,
    compute_crossings,
)
from vd3d.zone.geom import point_on_line_at_parameter
from vd3d.zone.types import PointNotInOpenFace, Zone


def face_containing_point(arrangement: Arrangement2D, point: Point2D) -> Face:
    """The unique open face that contains ``point``.

    Line-arrangement faces are convex open cells. A point on a supporting
    line belongs to no open face.
    """
    matches = [
        face
        for face in arrangement.faces
        if point_in_face(arrangement, face, point, closed=False)
    ]
    if len(matches) != 1:
        raise PointNotInOpenFace(
            f"point {point} lies in {len(matches)} open faces (expected 1)"
        )
    return matches[0]


def face_on_other_side(arrangement: Arrangement2D, face: Face, edge_id: int) -> Face:
    """The face incident to ``edge_id`` that is not ``face``."""
    sides = [
        he for he in arrangement.half_edges if he.edge_id == edge_id
    ]
    if len(sides) != 2:
        raise RuntimeError(f"edge {edge_id} does not have two half-edges")
    face_ids = [he.face_id for he in sides]
    if face_ids[0] is None or face_ids[1] is None:
        raise RuntimeError(f"edge {edge_id} is missing incident faces")
    if face_ids[0] == face_ids[1]:
        raise RuntimeError(f"edge {edge_id} bounds the same face twice")
    if face.id == face_ids[0]:
        return arrangement.faces[face_ids[1]]
    if face.id == face_ids[1]:
        return arrangement.faces[face_ids[0]]
    raise RuntimeError(f"face {face.id} is not incident to edge {edge_id}")


def compute_zone(arrangement: Arrangement2D, query: Line2D) -> Zone:
    """Walk the faces of ``arrangement`` crossed by ``query``.

    Start in the face that contains a point of ``query`` strictly before
    the first crossing (parameter ``t_min - 1``).

    Through a vertex: incident-edge hits at the same ``t`` collapse to
    one feature; the next face is located by sampling ``query`` just after
    that ``t`` (the line does not enter every wedge at the vertex).

    Overlap with an arrangement line: the zone is the incident faces of
    that line, in order of first appearance along ``query``.
    """
    coincident = coincident_line_indices(arrangement, query)
    if coincident:
        return _zone_along_supporting_line(arrangement, query, coincident[0])

    raw = compute_crossings(arrangement, query)
    crossings = collapse_crossings(raw)
    vertices = tuple(
        arrangement.vertices[c.vertex_id]
        for c in crossings
        if c.vertex_id is not None
    )

    if not crossings:
        sample, _ = query.sample_points()
        start = face_containing_point(arrangement, sample)
        return Zone(
            query=query,
            faces=(start,),
            edges=(),
            vertices=(),
            crossings=(),
        )

    before = point_on_line_at_parameter(query, crossings[0].t - 1)
    current = face_containing_point(arrangement, before)
    faces = [current]
    edges = []
    for index, crossing in enumerate(crossings):
        edges.append(arrangement.edges[crossing.edge_id])
        if crossing.vertex_id is not None:
            if index + 1 < len(crossings):
                sample_t = (crossing.t + crossings[index + 1].t) / 2
            else:
                sample_t = crossing.t + 1
            current = face_containing_point(
                arrangement, point_on_line_at_parameter(query, sample_t)
            )
        else:
            current = face_on_other_side(arrangement, current, crossing.edge_id)
        faces.append(current)

    return Zone(
        query=query,
        faces=tuple(faces),
        edges=tuple(edges),
        vertices=vertices,
        crossings=crossings,
    )


def _zone_along_supporting_line(
    arrangement: Arrangement2D, query: Line2D, line_index: int
) -> Zone:
    """Ordered incident faces of a line that coincides with ``query``."""
    crossings = compute_crossings(arrangement, query)
    edges = [edge for edge in arrangement.edges if edge.line_index == line_index]
    edges.sort(key=lambda edge: _edge_sort_key(arrangement, edge, query))
    faces: list[Face] = []
    seen: set[int] = set()
    for edge in edges:
        for he in arrangement.half_edges:
            if he.edge_id != edge.id or he.face_id is None:
                continue
            if he.face_id in seen:
                continue
            seen.add(he.face_id)
            faces.append(arrangement.faces[he.face_id])
    vertices = tuple(
        arrangement.vertices[c.vertex_id]
        for c in crossings
        if c.vertex_id is not None
    )
    return Zone(
        query=query,
        faces=tuple(faces),
        edges=tuple(edges),
        vertices=vertices,
        crossings=crossings,
    )


def _edge_sort_key(
    arrangement: Arrangement2D, edge: EdgePiece, query: Line2D
) -> tuple:
    """Sort key of an overlapped edge along ``query`` (a sample point's ``t``)."""
    line = arrangement.lines[edge.line_index]
    if edge.kind == "segment" and edge.t_min is not None and edge.t_max is not None:
        t = (edge.t_min + edge.t_max) / 2
    elif edge.kind == "ray" and edge.t_min is not None:
        t = edge.t_min + 1
    elif edge.kind == "ray" and edge.t_max is not None:
        t = edge.t_max - 1
    else:
        sample, _ = line.sample_points()
        return (parameter_on_line(query, sample), edge.id)
    point = point_on_line_at_parameter(line, t)
    return (parameter_on_line(query, point), edge.id)
