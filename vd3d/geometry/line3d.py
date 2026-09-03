"""A 3D line given by a point and a direction."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.geometry.linalg import cross
from vd3d.geometry.points import Point3D
from vd3d.geometry.scalar import Scalar, as_scalar

Vec3 = tuple[Scalar, Scalar, Scalar]


@dataclass(frozen=True, slots=True)
class Line3D:
    point: Point3D
    direction: Vec3
    plane_a: int | None = None
    plane_b: int | None = None
    id: int | None = None

    def __post_init__(self) -> None:
        dx, dy, dz = (as_scalar(c) for c in self.direction)
        object.__setattr__(self, "direction", (dx, dy, dz))
        if dx == 0 and dy == 0 and dz == 0:
            raise ValueError("zero direction")

    def point_at(self, t: int | Scalar | str) -> Point3D:
        t = as_scalar(t)
        dx, dy, dz = self.direction
        return Point3D(
            self.point.x + t * dx,
            self.point.y + t * dy,
            self.point.z + t * dz,
        )

    def contains(self, point: Point3D) -> bool:
        offset = (
            point.x - self.point.x,
            point.y - self.point.y,
            point.z - self.point.z,
        )
        return cross(offset, self.direction) == (0, 0, 0)


def directions_parallel(u: Vec3, v: Vec3) -> bool:
    return cross(u, v) == (0, 0, 0)
