"""Intersect a query line with every arrangement edge and sort along the query."""

from __future__ import annotations

from vd3d.arrangement2d.intersect import COINCIDENT, PARALLEL, intersect_lines_2d
from vd3d.arrangement2d.predicates import parameter_on_line
from vd3d.arrangement2d.types import Arrangement2D
from vd3d.geometry.line2d import Line2D
from vd3d.zone.geom import edge_contains_point, vertex_id_at
from vd3d.zone.types import Crossing, QueryOverlapsArrangement


def compute_crossings(arrangement: Arrangement2D, query: Line2D) -> tuple[Crossing, ...]:
    """Intersections of ``query`` with arrangement edges, sorted by parameter.

    Parallel edges contribute nothing. A coincident supporting line raises
    ``QueryOverlapsArrangement``. A hit at an arrangement vertex is recorded
    (``vertex_id`` set) on every incident edge that contains the point;
    ``compute_zone`` rejects that case.
    """
    crossings: list[Crossing] = []
    for edge in arrangement.edges:
        supporting = arrangement.lines[edge.line_index]
        result = intersect_lines_2d(query, supporting)
        if result is COINCIDENT:
            raise QueryOverlapsArrangement(
                f"query coincides with arrangement line {edge.line_index}"
            )
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
