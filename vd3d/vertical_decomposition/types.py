"""Types for the 2D vertical decomposition (y-parallel walls)."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.arrangement2d.types import Arrangement2D
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar


@dataclass(frozen=True, slots=True)
class VerticalRay:
    """A ray ``x = origin.x``, ``y`` increasing (``direction = +1``) or decreasing."""

    origin_vertex_id: int
    origin: Point2D
    direction: int  # +1 = +y, -1 = -y

    def __post_init__(self) -> None:
        if self.direction not in (1, -1):
            raise ValueError("direction must be +1 (+y) or -1 (-y)")


@dataclass(frozen=True, slots=True)
class Hit:
    """First intersection of a vertical ray with an arrangement edge, or unbounded."""

    unbounded: bool
    point: Point2D | None = None
    line_index: int | None = None
    edge_id: int | None = None
    vertex_id: int | None = None


UNBOUNDED = Hit(unbounded=True)


@dataclass(slots=True)
class VerticalWall:
    """A y-parallel wall ``x = const``, possibly unbounded in ``y``."""

    id: int
    x: Scalar
    y_min: Scalar | None  # None = -∞
    y_max: Scalar | None  # None = +∞
    bottom_vertex_id: int | None = None
    top_vertex_id: int | None = None
    bottom_line_index: int | None = None
    top_line_index: int | None = None


@dataclass(slots=True)
class VDCell2D:
    """A trapezoid: at most two vertical sides and two supporting-line sides."""

    id: int
    left_x: Scalar | None
    right_x: Scalar | None
    lower_line: int | None
    upper_line: int | None
    representative: Point2D
    source_face: int
    vertical_walls: tuple[int, ...] = ()
    neighbors: tuple[int, ...] = ()

    @property
    def unbounded(self) -> bool:
        return (
            self.left_x is None
            or self.right_x is None
            or self.lower_line is None
            or self.upper_line is None
        )


@dataclass(slots=True)
class VerticalDecomposition:
    """2D vertical decomposition of an ``Arrangement2D``."""

    arrangement: Arrangement2D
    cells: tuple[VDCell2D, ...]
    walls: tuple[VerticalWall, ...]
