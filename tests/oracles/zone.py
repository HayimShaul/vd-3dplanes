"""Brute-force zone checkers. Allowed to be slow and stupid."""

from __future__ import annotations

from vd3d.arrangement2d.intersect import COINCIDENT, PARALLEL, intersect_lines_2d
from vd3d.arrangement2d.predicates import parameter_on_line
from vd3d.arrangement2d.sampling import random_query_missing_vertices
from vd3d.arrangement2d.types import Arrangement2D
from vd3d.geometry.line2d import Line2D
from vd3d.zone.crossings import coincident_line_indices, collapse_crossings
from vd3d.zone.geom import edge_contains_point, point_on_line_at_parameter, vertex_id_at
from vd3d.zone.types import Crossing
from vd3d.zone.walk import face_containing_point

__all__ = [
    "brute_force_crossings",
    "brute_force_face_ids",
    "random_query_missing_vertices",
]


def brute_force_crossings(arrangement: Arrangement2D, query: Line2D) -> tuple[Crossing, ...]:
    """Intersect ``query`` with each supporting line, then attach edge pieces.

    Line-first, unlike ``compute_crossings`` which loops over edges. Same
    geometric answer. Overlap: vertices of the coincident line.
    """
    coincident = coincident_line_indices(arrangement, query)
    if coincident:
        found: list[Crossing] = []
        wanted = set(coincident)
        seen: set[int] = set()
        for vertex in arrangement.vertices:
            if vertex.id in seen:
                continue
            if not any(index in vertex.line_indices for index in wanted):
                continue
            edge_id = None
            for edge in arrangement.edges:
                if edge.line_index not in wanted:
                    continue
                if edge.start_vertex == vertex.id or edge.end_vertex == vertex.id:
                    edge_id = edge.id
                    break
            if edge_id is None:
                continue
            seen.add(vertex.id)
            found.append(
                Crossing(
                    point=vertex.point,
                    t=parameter_on_line(query, vertex.point),
                    edge_id=edge_id,
                    vertex_id=vertex.id,
                )
            )
        found.sort(key=lambda c: (c.t, c.edge_id))
        return tuple(found)

    found = []
    for line_index, line in enumerate(arrangement.lines):
        result = intersect_lines_2d(query, line)
        if result is COINCIDENT:
            continue
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
    Through a vertex, uses unique ``t`` values. Overlap: supporting-line
    incident faces in first-appearance order along the line.
    """
    coincident = coincident_line_indices(arrangement, query)
    if coincident:
        return _brute_force_overlap_face_ids(arrangement, query, coincident[0])

    crossings = brute_force_crossings(arrangement, query)
    collapsed = collapse_crossings(crossings)
    parameters = [c.t for c in collapsed]
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


def _brute_force_overlap_face_ids(
    arrangement: Arrangement2D, query: Line2D, line_index: int
) -> tuple[int, ...]:
    """Incident faces of the coincident line, first appearance along ``query``."""
    edges = [edge for edge in arrangement.edges if edge.line_index == line_index]
    edges.sort(key=lambda edge: _overlap_edge_key(arrangement, edge, query))
    ids: list[int] = []
    seen: set[int] = set()
    for edge in edges:
        for he in arrangement.half_edges:
            if he.edge_id != edge.id or he.face_id is None:
                continue
            if he.face_id in seen:
                continue
            seen.add(he.face_id)
            ids.append(he.face_id)
    return tuple(ids)


def _overlap_edge_key(arrangement: Arrangement2D, edge, query: Line2D) -> tuple:
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
