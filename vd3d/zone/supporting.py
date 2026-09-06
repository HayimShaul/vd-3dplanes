"""Full zone of a supporting line already in the arrangement.

``COMPUTE_ZONE`` cannot be used here: the query overlaps its own edges.
Instead collect every finite vertex of every face incident to the line,
and split those vertices into *on the line* versus *opposite*.
"""

from __future__ import annotations

from vd3d.arrangement2d.types import Arrangement2D, Face, Vertex
from vd3d.zone.types import SupportingLineZone


def compute_supporting_line_zone(
    arrangement: Arrangement2D, line_index: int
) -> SupportingLineZone:
    """Vertices of faces incident to arrangement line ``line_index``.

    ``vertices_on_line`` lie on the supporting line. ``opposite_vertices``
    are the other corners of those faces — the alignment candidates.
    """
    if line_index < 0 or line_index >= len(arrangement.lines):
        raise IndexError(f"line_index {line_index} out of range")

    line = arrangement.lines[line_index]
    incident: dict[int, Face] = {}
    for edge in arrangement.edges:
        if edge.line_index != line_index:
            continue
        for he in arrangement.half_edges:
            if he.edge_id != edge.id:
                continue
            if he.face_id is None:
                continue
            face = arrangement.faces[he.face_id]
            incident[face.id] = face

    faces = tuple(incident[i] for i in sorted(incident))
    seen: dict[int, Vertex] = {}
    for face in faces:
        for vid in _finite_vertex_ids(arrangement, face):
            seen[vid] = arrangement.vertices[vid]

    on_line: list[Vertex] = []
    opposite: list[Vertex] = []
    for vid in sorted(seen):
        vertex = seen[vid]
        if line.contains(vertex.point):
            on_line.append(vertex)
        else:
            opposite.append(vertex)

    return SupportingLineZone(
        line_index=line_index,
        faces=faces,
        vertices_on_line=tuple(on_line),
        opposite_vertices=tuple(opposite),
    )


def _finite_vertex_ids(arrangement: Arrangement2D, face: Face) -> set[int]:
    ids: set[int] = set()
    for he in arrangement.cycle(face):
        if he.origin_id is not None:
            ids.add(he.origin_id)
        twin = arrangement.half_edges[he.twin_id]
        if twin.origin_id is not None:
            ids.add(twin.origin_id)
    return ids
