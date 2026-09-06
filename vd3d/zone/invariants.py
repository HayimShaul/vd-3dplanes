"""Zone invariants checked by unit tests."""

from __future__ import annotations

from vd3d.arrangement2d.predicates import parameter_on_line
from vd3d.arrangement2d.types import Arrangement2D
from vd3d.zone.geom import edge_contains_point
from vd3d.zone.types import SupportingLineZone, Zone
from vd3d.zone.walk import face_on_other_side


def crossings_sorted(zone: Zone) -> bool:
    crossings = zone.crossings
    for left, right in zip(crossings, crossings[1:]):
        if (left.t, left.edge_id) > (right.t, right.edge_id):
            return False
        if left.t == right.t and left.edge_id == right.edge_id:
            return False
    return True


def crossings_lie_on_query_and_edge(arrangement: Arrangement2D, zone: Zone) -> bool:
    for crossing, edge in zip(zone.crossings, zone.edges):
        if crossing.edge_id != edge.id:
            return False
        if not zone.query.contains(crossing.point):
            return False
        if not edge_contains_point(arrangement, edge, crossing.point):
            return False
        if parameter_on_line(zone.query, crossing.point) != crossing.t:
            return False
    return True


def zone_faces_match_crossings(arrangement: Arrangement2D, zone: Zone) -> bool:
    if len(zone.faces) != len(zone.crossings) + 1:
        return False
    if len(zone.edges) != len(zone.crossings):
        return False
    for i, crossing in enumerate(zone.crossings):
        nxt = face_on_other_side(arrangement, zone.faces[i], crossing.edge_id)
        if nxt.id != zone.faces[i + 1].id:
            return False
        if zone.faces[i].id == zone.faces[i + 1].id:
            return False
    return True


def supporting_split_matches_line(
    arrangement: Arrangement2D, zone: SupportingLineZone
) -> bool:
    line = arrangement.lines[zone.line_index]
    on_ids = {v.id for v in zone.vertices_on_line}
    opp_ids = {v.id for v in zone.opposite_vertices}
    if on_ids & opp_ids:
        return False
    for vertex in zone.vertices_on_line:
        if not line.contains(vertex.point):
            return False
    for vertex in zone.opposite_vertices:
        if line.contains(vertex.point):
            return False
    return True


def verify_zone_invariants(arrangement: Arrangement2D, zone: Zone) -> None:
    checks = (
        ("crossings sorted along query", lambda: crossings_sorted(zone)),
        (
            "crossings on query and edge",
            lambda: crossings_lie_on_query_and_edge(arrangement, zone),
        ),
        (
            "faces match edge sides",
            lambda: zone_faces_match_crossings(arrangement, zone),
        ),
    )
    for name, pred in checks:
        if not pred():
            raise AssertionError(f"zone invariant failed: {name}")


def verify_supporting_line_zone(
    arrangement: Arrangement2D, zone: SupportingLineZone
) -> None:
    if not supporting_split_matches_line(arrangement, zone):
        raise AssertionError("supporting-line zone: on-line / opposite split is wrong")
    if not zone.faces:
        if arrangement.lines:
            raise AssertionError("supporting line has no incident faces")
