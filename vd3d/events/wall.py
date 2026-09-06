"""``x_L(z)`` and the y-parallel wall of an intersection line (design §8.1)."""

from __future__ import annotations

from vd3d.geometry.linalg import cross
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane
from vd3d.geometry.points import Point2D, Point3D
from vd3d.geometry.scalar import Scalar, as_scalar


class VerticalIntersectionLine(ValueError):
    """``L`` is parallel to the xy-plane (``dz = 0``), so ``x(z)`` is undefined."""


def line_xz_param(line: Line3D) -> tuple[Scalar, Scalar, Scalar, Scalar]:
    """Coefficients of ``x = ax z + bx`` and ``y = ay z + by``.

    Requires ``dz ≠ 0``.
    """
    dx, dy, dz = line.direction
    if dz == 0:
        raise VerticalIntersectionLine("intersection line is parallel to the xy-plane")
    ax = dx / dz
    ay = dy / dz
    bx = line.point.x - line.point.z * ax
    by = line.point.y - line.point.z * ay
    return ax, bx, ay, by


def x_of_line(line: Line3D, z: int | Scalar | str) -> Scalar:
    ax, bx, _ay, _by = line_xz_param(line)
    return ax * as_scalar(z) + bx


def y_of_line(line: Line3D, z: int | Scalar | str) -> Scalar:
    _ax, _bx, ay, by = line_xz_param(line)
    return ay * as_scalar(z) + by


def point_on_line_at_z(line: Line3D, z: int | Scalar | str) -> Point3D:
    z = as_scalar(z)
    return Point3D(x_of_line(line, z), y_of_line(line, z), z)


def build_vertical_wall(line: Line3D) -> Plane:
    """The unique y-parallel plane ``x = a z + b`` containing ``L``.

    Equation ``x - a z - b = 0`` (no ``y`` term). The wall id is
    ``-(line.id + 1)`` so it does not collide with input plane ids.
    """
    ax, bx, _ay, _by = line_xz_param(line)
    wall_id = -1 if line.id is None else -(line.id + 1)
    return Plane(id=wall_id, a=1, b=0, c=-ax, d=-bx)


def lift_wall_point(point: Point2D, wall: Plane) -> Point3D:
    """Map wall-frame ``(z, y)`` to 3D ``(a z + b, y, z)``.

    ``Point2D.x`` is the sweep coordinate ``z``. ``Point2D.y`` is 3D ``y``,
    so a 2D-vertical alignment (same ``x``, same ``z``) is a vertical
    segment in the wall plot.
    """
    if wall.a == 0 or wall.b != 0:
        raise ValueError("wall must be y-parallel with a ≠ 0")
    z = point.x
    y = point.y
    x = -(wall.c * z + wall.d) / wall.a
    return Point3D(x, y, z)


def wall_frame_point(line: Line3D, z: int | Scalar | str) -> Point2D:
    """``L`` at height ``z`` in the wall frame ``(z, y)``."""
    z = as_scalar(z)
    return Point2D(z, y_of_line(line, z))


def alignment_z(left: Line3D, right: Line3D) -> Scalar | None:
    """Solve ``x_left(z) = x_right(z)``. ``None`` if the walls are parallel."""
    ax1, bx1, _ay1, _by1 = line_xz_param(left)
    ax2, bx2, _ay2, _by2 = line_xz_param(right)
    if ax1 == ax2:
        return None
    return (bx2 - bx1) / (ax1 - ax2)


def lines_meet(left: Line3D, right: Line3D) -> bool:
    """True iff the supporting lines intersect (or coincide)."""
    cr = cross(left.direction, right.direction)
    if cr == (0, 0, 0):
        return left.contains(right.point)
    offset = (
        right.point.x - left.point.x,
        right.point.y - left.point.y,
        right.point.z - left.point.z,
    )
    return offset[0] * cr[0] + offset[1] * cr[1] + offset[2] * cr[2] == 0
