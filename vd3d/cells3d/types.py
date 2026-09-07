"""3D cells grown during the z-sweep (design §1.5–1.6)."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.geometry.plane import Plane
from vd3d.geometry.scalar import Scalar


@dataclass(slots=True)
class Cell3D:
    """A prism: at most one floor, one ceiling, and four vertical walls.

    ``floor`` / ``ceiling`` are the lower / upper supporting input planes of
    the tracked 2D trapezoid (or ``None`` if that side is unbounded).
    ``vertical_walls`` are Steiner planes parallel to the y-axis (the 3D
    walls of the 2D ``x = const`` sides). ``lower_z`` / ``upper_z`` are
    ``None`` for ±∞.
    """

    id: int
    floor: Plane | None
    ceiling: Plane | None
    vertical_walls: tuple[Plane, ...]
    lower_z: Scalar | None
    upper_z: Scalar | None


@dataclass(slots=True)
class ActiveCell:
    """The currently growing 3D cell attached to one 2D cell signature."""

    cell3d_id: int
    current_2d_cell_id: int
    start_z: Scalar | None
    current_floor: Plane | None
    current_ceiling: Plane | None
    current_vertical_walls: tuple[Plane, ...]
