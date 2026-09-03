"""Vertices and undirected edge pieces of a 2D line arrangement.

No half-edges yet: each supporting line is split at its incident vertices
into two rays plus interior segments (or one whole-line piece if it has
no vertices).
"""

from __future__ import annotations

from vd3d.arrangement2d.intersect import COINCIDENT, PARALLEL, intersect_lines_2d
from vd3d.arrangement2d.predicates import parameter_on_line
from vd3d.arrangement2d.types import EdgePiece, Vertex
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D


def compute_vertices(lines: tuple[Line2D, ...]) -> tuple[Vertex, ...]:
    """Pairwise intersections, merged by exact position.

    Concurrent triples (and higher) become one vertex. Coincident input
    lines are rejected; parallels simply contribute no vertex.
    """
    points: dict[Point2D, None] = {}
    for i, left in enumerate(lines):
        for j in range(i + 1, len(lines)):
            result = intersect_lines_2d(left, lines[j])
            if result is COINCIDENT:
                raise ValueError(
                    f"coincident lines are not supported (indices {i}, {j})"
                )
            if result is PARALLEL:
                continue
            points[result] = None

    vertices: list[Vertex] = []
    for point in sorted(points, key=lambda p: (p.x, p.y)):
        incident = tuple(i for i, line in enumerate(lines) if line.contains(point))
        if not incident:
            raise RuntimeError("intersection point lies on no input line")
        vertices.append(
            Vertex(id=len(vertices), point=point, line_indices=incident)
        )
    return tuple(vertices)


def compute_edge_pieces(
    lines: tuple[Line2D, ...],
    vertices: tuple[Vertex, ...],
) -> tuple[EdgePiece, ...]:
    """Split each line at the vertices that lie on it, in increasing ``t``."""
    pieces: list[EdgePiece] = []
    for line_index, line in enumerate(lines):
        on_line = [v for v in vertices if line_index in v.line_indices]
        on_line.sort(key=lambda v: parameter_on_line(line, v.point))
        pieces.extend(_split_one_line(line_index, line, on_line, start_id=len(pieces)))
    return tuple(pieces)


def _split_one_line(
    line_index: int,
    line: Line2D,
    on_line: list[Vertex],
    *,
    start_id: int,
) -> list[EdgePiece]:
    if not on_line:
        return [
            EdgePiece(
                id=start_id,
                kind="line",
                line_index=line_index,
                start_vertex=None,
                end_vertex=None,
                t_min=None,
                t_max=None,
            )
        ]

    out: list[EdgePiece] = []
    first = on_line[0]
    out.append(
        EdgePiece(
            id=start_id,
            kind="ray",
            line_index=line_index,
            start_vertex=None,
            end_vertex=first.id,
            t_min=None,
            t_max=parameter_on_line(line, first.point),
        )
    )
    for left, right in zip(on_line, on_line[1:]):
        out.append(
            EdgePiece(
                id=start_id + len(out),
                kind="segment",
                line_index=line_index,
                start_vertex=left.id,
                end_vertex=right.id,
                t_min=parameter_on_line(line, left.point),
                t_max=parameter_on_line(line, right.point),
            )
        )
    last = on_line[-1]
    out.append(
        EdgePiece(
            id=start_id + len(out),
            kind="ray",
            line_index=line_index,
            start_vertex=last.id,
            end_vertex=None,
            t_min=parameter_on_line(line, last.point),
            t_max=None,
        )
    )
    return out
