"""Geometry invariants from design §21, as exact predicates."""

from __future__ import annotations

from vd3d.geometry.line2d import Line2D
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane
from vd3d.geometry.points import Point3D
from vd3d.geometry.scalar import Scalar, as_scalar


def line_lies_on_plane(line: Line3D, plane: Plane) -> bool:
    """The supporting line lies on ``plane``: a point and the direction do."""
    if not plane.contains(line.point):
        return False
    dx, dy, dz = line.direction
    moved = Point3D(line.point.x + dx, line.point.y + dy, line.point.z + dz)
    return plane.contains(moved)


def line_lies_on_both_planes(line: Line3D, p: Plane, q: Plane) -> bool:
    return line_lies_on_plane(line, p) and line_lies_on_plane(line, q)


def point_lies_on_planes(point: Point3D, *planes: Plane) -> bool:
    return all(plane.contains(point) for plane in planes)


def slice_lifts_to_plane(line: Line2D, z: int | Scalar | str, plane: Plane) -> bool:
    """Every sample 2D point on ``line`` lifts to a 3D point on ``plane``."""
    z = as_scalar(z)
    for point2d in line.sample_points():
        lifted = Point3D(point2d.x, point2d.y, z)
        if not plane.contains(lifted):
            return False
        if not line.contains(point2d):
            return False
    return True
