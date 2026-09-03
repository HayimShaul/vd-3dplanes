"""Exact 2D and 3D points."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.geometry.scalar import Scalar, as_scalar


@dataclass(frozen=True, slots=True)
class Point2D:
    x: Scalar
    y: Scalar

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", as_scalar(self.x))
        object.__setattr__(self, "y", as_scalar(self.y))


@dataclass(frozen=True, slots=True)
class Point3D:
    x: Scalar
    y: Scalar
    z: Scalar

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", as_scalar(self.x))
        object.__setattr__(self, "y", as_scalar(self.y))
        object.__setattr__(self, "z", as_scalar(self.z))

    def translated(self, dx: int | Scalar | str, dy: int | Scalar | str, dz: int | Scalar | str) -> Point3D:
        return Point3D(
            self.x + as_scalar(dx),
            self.y + as_scalar(dy),
            self.z + as_scalar(dz),
        )
