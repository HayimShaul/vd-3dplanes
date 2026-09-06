"""Brute-force zone checkers. Allowed to be slow and stupid."""

from __future__ import annotations

from vd3d.arrangement2d.intersect import COINCIDENT, PARALLEL, intersect_lines_2d
from vd3d.arrangement2d.predicates import parameter_on_line
from vd3d.arrangement2d.sampling import random_query_missing_vertices
from vd3d.arrangement2d.types import Arrangement2D
from vd3d.geometry.line2d import Line2D
from vd3d.zone.geom import edge_contains_point, point_on_line_at_parameter, vertex_id_at
from vd3d.zone.types import Crossing, QueryOverlapsArrangement
from vd3d.zone.walk import face_containing_point

__all__ = [
    "brute_force_crossings",
    "brute_force_face_ids",
    "random_query_missing_vertices",
]


def brute_force_crossings(arrangement: Arrangement2D, query: Line2D) -> tuple[Crossing, ...]:
    """Intersect ``query`` with each supporting line, then attach edge pieces.

    Line-first, unlike ``compute_crossings`` which loops over edges. Same
    geometric answer under general position.
    """
    found: list[Crossing] = []
    for line_index, line in enumerate(arrangement.lines):
        result = intersect_lines_2d(query, line)
        if result is COINCIDENT:
            raise QueryOverlapsArrangement(
                f"query coincides with arrangement line {line_index}"
            )
        if result is PARALLEL:
            continue
        for edge in arrangement.edges:
            if edge.line_index != line_index:
                continue
            if not edge_contains_point(arrangement, edge, result):
                continue
            found.append(
                Crossing(
                    point=result,
                    t=parameter_on_line(query, result),
                    edge_id=edge.id,
                    vertex_id=vertex_id_at(arrangement, result),
                )
            )
    found.sort(key=lambda c: (c.t, c.edge_id))
    return tuple(found)


def brute_force_face_ids(arrangement: Arrangement2D, query: Line2D) -> tuple[int, ...]:
    """Locate the face at a sample of ``query`` in each open interval between crossings.

    Independent of ``FACE_ON_OTHER_SIDE``: midpoint-in-face on the query.
    """
    crossings = brute_force_crossings(arrangement, query)
    if any(c.vertex_id is not None for c in crossings):
        raise ValueError("brute-force face sequence is undefined through a vertex")
    parameters = [c.t for c in crossings]
    if not parameters:
        sample, _ = query.sample_points()
        return (face_containing_point(arrangement, sample).id,)

    samples = [parameters[0] - 1]
    for left, right in zip(parameters, parameters[1:]):
        samples.append((left + right) / 2)
    samples.append(parameters[-1] + 1)
    return tuple(
        face_containing_point(arrangement, point_on_line_at_parameter(query, t)).id
        for t in samples
    )
