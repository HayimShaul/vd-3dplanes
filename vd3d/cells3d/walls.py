"""3D vertical walls of a 2D trapezoid (design §14 ``EXTRACT_VERTICAL_WALLS``)."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events.wall import build_vertical_wall
from vd3d.geometry.intersections import PARALLEL, intersect_planes
from vd3d.geometry.plane import Plane
from vd3d.geometry.scalar import Scalar
from vd3d.vertical_decomposition.cells import y_range_at
from vd3d.vertical_decomposition.geom import line_is_vertical, x_of_vertical_line
from vd3d.vertical_decomposition.types import VDCell2D, VerticalDecomposition, VerticalWall


def supporting_plane_ids(
    vd: VerticalDecomposition, cell: VDCell2D
) -> tuple[int | None, int | None]:
    """Source plane ids of the lower and upper supporting lines."""
    lines = vd.arrangement.lines
    lower = None if cell.lower_line is None else lines[cell.lower_line].source_plane_id
    upper = None if cell.upper_line is None else lines[cell.upper_line].source_plane_id
    return lower, upper


def supporting_planes(
    vd: VerticalDecomposition,
    cell: VDCell2D,
    planes: Sequence[Plane],
) -> tuple[Plane | None, Plane | None]:
    lookup = {plane.id: plane for plane in planes}
    lower_id, upper_id = supporting_plane_ids(vd, cell)
    return (
        None if lower_id is None else lookup[lower_id],
        None if upper_id is None else lookup[upper_id],
    )


def side_planes(
    vd: VerticalDecomposition, cell: VDCell2D, x: Scalar | None
) -> tuple[int, ...]:
    """Source plane ids of vertices/lines on this cell's vertical side at ``x``.

    Vertices that share ``x`` but lie outside the cell's y-span (a blocked
    alignment) are ignored, so the signature stays stable between events.
    """
    if x is None:
        return ()
    lo, hi = y_range_at(cell, x, vd.arrangement.lines)
    ids: set[int] = set()
    for vertex in vd.arrangement.vertices:
        if vertex.point.x != x:
            continue
        y = vertex.point.y
        if lo is not None and y < lo:
            continue
        if hi is not None and y > hi:
            continue
        for line_index in vertex.line_indices:
            pid = vd.arrangement.lines[line_index].source_plane_id
            if pid is not None:
                ids.add(pid)
    for line in vd.arrangement.lines:
        if line_is_vertical(line) and x_of_vertical_line(line) == x:
            if line.source_plane_id is not None:
                ids.add(line.source_plane_id)
    return tuple(sorted(ids))


def planes_at_x(vd: VerticalDecomposition, x: Scalar | None) -> tuple[int, ...]:
    """Source plane ids of arrangement features at a vertical line ``x = const``."""
    if x is None:
        return ()
    ids: set[int] = set()
    for vertex in vd.arrangement.vertices:
        if vertex.point.x != x:
            continue
        for line_index in vertex.line_indices:
            pid = vd.arrangement.lines[line_index].source_plane_id
            if pid is not None:
                ids.add(pid)
    for line in vd.arrangement.lines:
        if line_is_vertical(line) and x_of_vertical_line(line) == x:
            if line.source_plane_id is not None:
                ids.add(line.source_plane_id)
    return tuple(sorted(ids))


def extract_vertical_walls(
    vd: VerticalDecomposition,
    cell: VDCell2D,
    planes: Sequence[Plane],
) -> tuple[Plane, ...]:
    """Steiner (or input) planes for the 2D cell's positive-length vertical sides."""
    lookup = {plane.id: plane for plane in planes}
    seen: dict[tuple[Scalar, Scalar, Scalar, Scalar], Plane] = {}
    for wall_id in cell.vertical_walls:
        wall = vd.walls[wall_id]
        plane = _wall_to_3d_plane(vd, wall, lookup)
        if plane is None:
            continue
        key = (plane.a, plane.b, plane.c, plane.d)
        seen.setdefault(key, plane)
    return tuple(seen.values())


def merge_planes(*groups: Sequence[Plane]) -> tuple[Plane, ...]:
    """Union of planes, identified by exact ``(a, b, c, d)``."""
    seen: dict[tuple[Scalar, Scalar, Scalar, Scalar], Plane] = {}
    for group in groups:
        for plane in group:
            key = (plane.a, plane.b, plane.c, plane.d)
            seen.setdefault(key, plane)
    return tuple(seen.values())


def _wall_to_3d_plane(
    vd: VerticalDecomposition,
    wall: VerticalWall,
    lookup: dict[int, Plane],
) -> Plane | None:
    ids: set[int] = set()
    for vertex_id in (wall.bottom_vertex_id, wall.top_vertex_id):
        if vertex_id is None:
            continue
        vertex = vd.arrangement.vertices[vertex_id]
        for line_index in vertex.line_indices:
            pid = vd.arrangement.lines[line_index].source_plane_id
            if pid is not None:
                ids.add(pid)
    if len(ids) >= 2:
        ordered = sorted(ids)
        line = intersect_planes(lookup[ordered[0]], lookup[ordered[1]])
        if line is PARALLEL:
            return None
        return build_vertical_wall(line)
    for line_index in (wall.bottom_line_index, wall.top_line_index):
        if line_index is None:
            continue
        pid = vd.arrangement.lines[line_index].source_plane_id
        if pid is not None:
            return lookup[pid]
    for line in vd.arrangement.lines:
        if line_is_vertical(line) and x_of_vertical_line(line) == wall.x:
            if line.source_plane_id is not None:
                return lookup[line.source_plane_id]
    return None
