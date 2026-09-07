"""3D-cell invariants from design §21."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.cells3d.types import Cell3D
from vd3d.geometry.plane import Plane


def floor_at_most_one(cell: Cell3D) -> bool:
    return cell.floor is None or isinstance(cell.floor, Plane)


def ceiling_at_most_one(cell: Cell3D) -> bool:
    return cell.ceiling is None or isinstance(cell.ceiling, Plane)


def vertical_walls_at_most_four(cell: Cell3D) -> bool:
    return len(cell.vertical_walls) <= 4


def z_extent_ordered(cell: Cell3D) -> bool:
    """Finite z-interval is strictly ordered. Unbounded sides are allowed."""
    if cell.lower_z is None or cell.upper_z is None:
        return True
    return cell.lower_z < cell.upper_z


def verify_cell3d(cell: Cell3D) -> None:
    checks = (
        ("floor is none or one plane", lambda: floor_at_most_one(cell)),
        ("ceiling is none or one plane", lambda: ceiling_at_most_one(cell)),
        ("vertical walls <= 4", lambda: vertical_walls_at_most_four(cell)),
        ("lower_z < upper_z when both finite", lambda: z_extent_ordered(cell)),
    )
    for name, pred in checks:
        if not pred():
            raise AssertionError(f"3D-cell invariant failed: {name} (cell {cell.id})")


def verify_cells3d(cells: Sequence[Cell3D]) -> None:
    for cell in cells:
        verify_cell3d(cell)
