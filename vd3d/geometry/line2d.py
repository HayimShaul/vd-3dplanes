"""A 2D line given by ``a*x + b*y + c = 0``."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar, as_scalar


@dataclass(frozen=True, slots=True)
class Line2D:
    a: Scalar
    b: Scalar
    c: Scalar
    source_plane_id: int | None = None
    id: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "a", as_scalar(self.a))
        object.__setattr__(self, "b", as_scalar(self.b))
        object.__setattr__(self, "c", as_scalar(self.c))
        if self.a == 0 and self.b == 0:
            raise ValueError("degenerate 2D line: a = b = 0")

    def eval(self, point: Point2D) -> Scalar:
        return self.a * point.x + self.b * point.y + self.c

    def contains(self, point: Point2D) -> bool:
        return self.eval(point) == 0

    def sample_points(self) -> tuple[Point2D, Point2D]:
        """Two distinct points on the line, for lifting and plotting."""
        if self.b != 0:
            return (
                Point2D(0, -self.c / self.b),
                Point2D(1, -(self.a + self.c) / self.b),
            )
        return (
            Point2D(-self.c / self.a, 0),
            Point2D(-(self.b + self.c) / self.a, 1),
        )
