"""``FACE_CONTAINING_POINT``, ``FACE_ON_OTHER_SIDE``, and ``COMPUTE_ZONE``."""

from __future__ import annotations

from vd3d.arrangement2d.invariants import point_in_face
from vd3d.arrangement2d.types import Arrangement2D, Face
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.zone.crossings import compute_crossings
from vd3d.zone.geom import point_on_line_at_parameter
from vd3d.zone.types import PointNotInOpenFace, QueryThroughVertex, Zone


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

    General-position query: not coincident with an arrangement line, and
    not through a vertex. Start in the face that contains a point of
    ``query`` strictly before the first crossing (parameter ``t_min - 1``).
    """
    crossings = compute_crossings(arrangement, query)
    if any(c.vertex_id is not None for c in crossings):
        raise QueryThroughVertex("query line passes through an arrangement vertex")

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
    for crossing in crossings:
        edges.append(arrangement.edges[crossing.edge_id])
        current = face_on_other_side(arrangement, current, crossing.edge_id)
        faces.append(current)

    return Zone(
        query=query,
        faces=tuple(faces),
        edges=tuple(edges),
        vertices=(),
        crossings=crossings,
    )
