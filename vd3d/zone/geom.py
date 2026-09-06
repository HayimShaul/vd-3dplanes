"""Exact helpers for walking a query line. No floats."""

from __future__ import annotations

from vd3d.arrangement2d.predicates import line_direction, parameter_on_line
from vd3d.arrangement2d.types import Arrangement2D, EdgePiece
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar, as_scalar


def point_on_line_at_parameter(line: Line2D, t: int | Scalar | str) -> Point2D:
    """The unique point of ``line`` whose 1D parameter equals ``t``.

    Parameter is ``(-b, a) · (x, y)``, same as ``parameter_on_line``.
    """
    target = as_scalar(t)
    dx, dy = line_direction(line)
    sample, _ = line.sample_points()
    t0 = parameter_on_line(line, sample)
    denom = dx * dx + dy * dy
    step = (target - t0) / denom
    return Point2D(sample.x + step * dx, sample.y + step * dy)


def edge_contains_point(arrangement: Arrangement2D, edge: EdgePiece, point: Point2D) -> bool:
    """True iff ``point`` lies on the (closed) edge piece."""
    line = arrangement.lines[edge.line_index]
    if not line.contains(point):
        return False
    t = parameter_on_line(line, point)
    if edge.kind == "line":
        return True
    if edge.kind == "ray":
        if edge.t_min is None and edge.t_max is not None:
            return t <= edge.t_max
        if edge.t_max is None and edge.t_min is not None:
            return t >= edge.t_min
        return False
    if edge.kind == "segment":
        if edge.t_min is None or edge.t_max is None:
            return False
        return edge.t_min <= t <= edge.t_max
    raise ValueError(f"unknown edge kind {edge.kind!r}")


def vertex_id_at(arrangement: Arrangement2D, point: Point2D) -> int | None:
    for vertex in arrangement.vertices:
        if vertex.point == point:
            return vertex.id
    return None
