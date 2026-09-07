"""Intersect a query line with every arrangement edge and sort along the query."""

from __future__ import annotations

from vd3d.arrangement2d.intersect import COINCIDENT, PARALLEL, intersect_lines_2d
from vd3d.arrangement2d.predicates import parameter_on_line
from vd3d.arrangement2d.types import Arrangement2D, EdgePiece
from vd3d.geometry.line2d import Line2D
from vd3d.zone.geom import edge_contains_point, vertex_id_at
from vd3d.zone.types import Crossing


def coincident_line_indices(arrangement: Arrangement2D, query: Line2D) -> tuple[int, ...]:
    """Arrangement lines whose equation is the same as ``query`` (any orientation)."""
    found: list[int] = []
    for index, line in enumerate(arrangement.lines):
        if intersect_lines_2d(query, line) is COINCIDENT:
            found.append(index)
    return tuple(found)


def compute_crossings(arrangement: Arrangement2D, query: Line2D) -> tuple[Crossing, ...]:
    """Intersections of ``query`` with arrangement edges, sorted by parameter.

    Parallel edges contribute nothing. A coincident supporting line is an
    overlap: crossings are the vertices of that line (one per vertex), not
    a continuum of edge points. A hit at an arrangement vertex is recorded
    (``vertex_id`` set) on every incident edge that contains the point;
    ``compute_zone`` collapses those to one feature per ``t``.
    """
    coincident = coincident_line_indices(arrangement, query)
    if coincident:
        return _overlap_crossings(arrangement, query, coincident)

    crossings: list[Crossing] = []
    for edge in arrangement.edges:
        supporting = arrangement.lines[edge.line_index]
        result = intersect_lines_2d(query, supporting)
        if result is COINCIDENT:
            continue
        if result is PARALLEL:
            continue
        if not edge_contains_point(arrangement, edge, result):
            continue
        crossings.append(
            Crossing(
                point=result,
                t=parameter_on_line(query, result),
                edge_id=edge.id,
                vertex_id=vertex_id_at(arrangement, result),
            )
        )
    crossings.sort(key=lambda c: (c.t, c.edge_id))
    return tuple(crossings)


def collapse_crossings(crossings: tuple[Crossing, ...]) -> tuple[Crossing, ...]:
    """One crossing per distinct parameter ``t``.

    Several edges meet at a vertex, so a through-vertex query records many
    raw crossings at the same ``t``. The zone walk uses one feature there.
    """
    if not crossings:
        return ()
    groups: list[list[Crossing]] = []
    for crossing in crossings:
        if groups and groups[-1][0].t == crossing.t:
            groups[-1].append(crossing)
        else:
            groups.append([crossing])
    collapsed: list[Crossing] = []
    for group in groups:
        vertex_hits = [c for c in group if c.vertex_id is not None]
        pool = vertex_hits if vertex_hits else group
        chosen = min(pool, key=lambda c: (c.vertex_id is None, c.vertex_id or 0, c.edge_id))
        collapsed.append(chosen)
    return tuple(collapsed)


def _overlap_crossings(
    arrangement: Arrangement2D,
    query: Line2D,
    line_indices: tuple[int, ...],
) -> tuple[Crossing, ...]:
    """Vertices on each coincident supporting line, as crossings along ``query``."""
    wanted = set(line_indices)
    crossings: list[Crossing] = []
    seen: set[int] = set()
    for vertex in arrangement.vertices:
        if vertex.id in seen:
            continue
        if not any(index in vertex.line_indices for index in wanted):
            continue
        edge = _incident_edge_on_lines(arrangement, vertex, wanted)
        if edge is None:
            continue
        seen.add(vertex.id)
        crossings.append(
            Crossing(
                point=vertex.point,
                t=parameter_on_line(query, vertex.point),
                edge_id=edge.id,
                vertex_id=vertex.id,
            )
        )
    crossings.sort(key=lambda c: (c.t, c.edge_id))
    return tuple(crossings)


def _incident_edge_on_lines(
    arrangement: Arrangement2D,
    vertex,
    line_indices: set[int],
) -> EdgePiece | None:
    for edge in arrangement.edges:
        if edge.line_index not in line_indices:
            continue
        if edge.start_vertex == vertex.id or edge.end_vertex == vertex.id:
            return edge
    return None
