"""Start / continue / end a 3D cell (design §14–16)."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.cells3d.invariants import verify_cell3d
from vd3d.cells3d.types import ActiveCell, Cell3D
from vd3d.cells3d.walls import extract_vertical_walls, merge_planes, supporting_planes
from vd3d.geometry.plane import Plane
from vd3d.geometry.scalar import Scalar, as_scalar
from vd3d.vertical_decomposition.types import VDCell2D, VerticalDecomposition


def start_3d_cell(
    cells: list[Cell3D],
    cell2d: VDCell2D,
    vd: VerticalDecomposition,
    planes: Sequence[Plane],
    event_z: int | Scalar | str | None,
) -> ActiveCell:
    """``START_3D_CELL``: a new prism whose floor is the 2D lower supporting plane."""
    floor, ceiling = supporting_planes(vd, cell2d, planes)
    walls = extract_vertical_walls(vd, cell2d, planes)
    lower_z = None if event_z is None else as_scalar(event_z)
    cell = Cell3D(
        id=len(cells),
        floor=floor,
        ceiling=ceiling,
        vertical_walls=walls,
        lower_z=lower_z,
        upper_z=None,
    )
    verify_cell3d(cell)
    cells.append(cell)
    return ActiveCell(
        cell3d_id=cell.id,
        current_2d_cell_id=cell2d.id,
        start_z=lower_z,
        current_floor=floor,
        current_ceiling=ceiling,
        current_vertical_walls=walls,
    )


def continue_3d_cell(
    cells: list[Cell3D],
    active: ActiveCell,
    cell2d_after: VDCell2D,
    vd_after: VerticalDecomposition,
    planes: Sequence[Plane],
) -> None:
    """``CONTINUE_3D_CELL``: same 3D cell, possibly updated walls."""
    cell = cells[active.cell3d_id]
    floor, ceiling = supporting_planes(vd_after, cell2d_after, planes)
    if (floor is None) != (cell.floor is None) or (
        floor is not None and cell.floor is not None and floor.id != cell.floor.id
    ):
        raise AssertionError(
            f"continued cell {cell.id} changed floor {cell.floor} -> {floor}"
        )
    if (ceiling is None) != (cell.ceiling is None) or (
        ceiling is not None and cell.ceiling is not None and ceiling.id != cell.ceiling.id
    ):
        raise AssertionError(
            f"continued cell {cell.id} changed ceiling {cell.ceiling} -> {ceiling}"
        )
    walls = extract_vertical_walls(vd_after, cell2d_after, planes)
    merged = merge_planes(cell.vertical_walls, walls)
    cell.vertical_walls = merged
    active.current_2d_cell_id = cell2d_after.id
    active.current_floor = floor
    active.current_ceiling = ceiling
    active.current_vertical_walls = walls
    verify_cell3d(cell)


def end_3d_cell(
    cells: list[Cell3D],
    active: ActiveCell,
    event_z: int | Scalar | str | None,
) -> Cell3D:
    """``END_3D_CELL``: freeze ``upper_z`` (``None`` = +∞)."""
    cell = cells[active.cell3d_id]
    cell.upper_z = None if event_z is None else as_scalar(event_z)
    verify_cell3d(cell)
    return cell
