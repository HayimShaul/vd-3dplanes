"""Planes given by ``a*x + b*y + c*z + d = 0``."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.geometry.linalg import cross
from vd3d.geometry.points import Point3D
from vd3d.geometry.scalar import Scalar, as_scalar

Vec3 = tuple[Scalar, Scalar, Scalar]


@dataclass(frozen=True, slots=True)
class Plane:
    """Oriented plane ``a*x + b*y + c*z + d = 0``.

    ``eval(p)`` is the residual at ``p``. Its sign is the orientation:
    positive, zero (on the plane), or negative. Coefficients are stored
    exactly as given; they are not normalized.
    """

    id: int
    a: Scalar
    b: Scalar
    c: Scalar
    d: Scalar

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", as_scalar(self.a))
        object.__setattr__(self, "b", as_scalar(self.b))
        object.__setattr__(self, "c", as_scalar(self.c))
        object.__setattr__(self, "d", as_scalar(self.d))
        if self.a == 0 and self.b == 0 and self.c == 0:
            raise ValueError("zero plane: a = b = c = 0")

    @property
    def normal(self) -> Vec3:
        return (self.a, self.b, self.c)

    def eval(self, point: Point3D) -> Scalar:
        return self.a * point.x + self.b * point.y + self.c * point.z + self.d

    def contains(self, point: Point3D) -> bool:
        return self.eval(point) == 0


def normals_parallel(p: Plane, q: Plane) -> bool:
    """True if the normals are parallel or anti-parallel (cross product is 0)."""
    return cross(p.normal, q.normal) == (0, 0, 0)
