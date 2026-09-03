"""``VERTICAL_RAYS_FROM`` and ``FIRST_HIT`` against arrangement edges."""

from __future__ import annotations

from vd3d.arrangement2d.predicates import parameter_on_line
from vd3d.arrangement2d.types import Arrangement2D, Vertex
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar
from vd3d.vertical_decomposition.geom import line_is_vertical, y_at_x
from vd3d.vertical_decomposition.types import UNBOUNDED, Hit, VerticalRay


def vertical_rays_from(vertex: Vertex) -> tuple[VerticalRay, VerticalRay]:
    """From ``vertex``, the two y-parallel rays ``+y`` and ``-y`` (``x`` fixed)."""
    plus = VerticalRay(
        origin_vertex_id=vertex.id, origin=vertex.point, direction=1
    )
    minus = VerticalRay(
        origin_vertex_id=vertex.id, origin=vertex.point, direction=-1
    )
    return plus, minus


def first_hit(arrangement: Arrangement2D, ray: VerticalRay) -> Hit:
    """First intersection of ``ray`` with an arrangement edge, or ``UNBOUNDED``.

    Lines through the origin are not hits: they are already incident to the
    vertex. A vertical line at a different ``x`` never meets the ray.
    """
    origin = ray.origin
    best_dist: Scalar | None = None
    best_y: Scalar | None = None
    best_lines: list[int] = []

    for i, line in enumerate(arrangement.lines):
        if line.contains(origin):
            continue
        if line_is_vertical(line):
            continue
        y_hit = y_at_x(line, origin.x)
        dy = y_hit - origin.y
        if ray.direction > 0:
            if dy <= 0:
                continue
            dist = dy
        else:
            if dy >= 0:
                continue
            dist = -dy
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best_y = y_hit
            best_lines = [i]
        elif dist == best_dist:
            best_lines.append(i)

    if best_y is None:
        return UNBOUNDED

    point = Point2D(origin.x, best_y)
    line_index = min(best_lines)
    vertex_id = _vertex_at(arrangement, point)
    edge_id = _edge_containing(arrangement, line_index, point)
    return Hit(
        unbounded=False,
        point=point,
        line_index=line_index,
        edge_id=edge_id,
        vertex_id=vertex_id,
    )


def _vertex_at(arrangement: Arrangement2D, point: Point2D) -> int | None:
    for vertex in arrangement.vertices:
        if vertex.point == point:
            return vertex.id
    return None


def _edge_containing(
    arrangement: Arrangement2D, line_index: int, point: Point2D
) -> int | None:
    line = arrangement.lines[line_index]
    t = parameter_on_line(line, point)
    for edge in arrangement.edges:
        if edge.line_index != line_index:
            continue
        if edge.kind == "line":
            return edge.id
        if edge.kind == "ray":
            if edge.t_min is None and edge.t_max is not None and t <= edge.t_max:
                return edge.id
            if edge.t_max is None and edge.t_min is not None and t >= edge.t_min:
                return edge.id
        if edge.kind == "segment":
            if edge.t_min is None or edge.t_max is None:
                continue
            if edge.t_min <= t <= edge.t_max:
                return edge.id
    return None
