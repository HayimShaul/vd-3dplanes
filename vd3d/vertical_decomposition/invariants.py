"""Vertical-decomposition invariants from design §21, as exact predicates."""

from __future__ import annotations

from vd3d.arrangement2d.invariants import point_in_face
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar
from vd3d.vertical_decomposition.cells import (
    interval_positive,
    y_range_at,
)
from vd3d.vertical_decomposition.geom import y_at_x
from vd3d.vertical_decomposition.types import VDCell2D, VerticalDecomposition


def point_in_cell(
    vd: VerticalDecomposition,
    cell: VDCell2D,
    point: Point2D,
    *,
    closed: bool = False,
) -> bool:
    """True iff ``point`` lies in the trapezoid (open or closed)."""
    lines = vd.arrangement.lines
    if cell.left_x is not None:
        if closed:
            if point.x < cell.left_x:
                return False
        elif point.x <= cell.left_x:
            return False
    if cell.right_x is not None:
        if closed:
            if point.x > cell.right_x:
                return False
        elif point.x >= cell.right_x:
            return False
    if cell.lower_line is not None:
        y_line = y_at_x(lines[cell.lower_line], point.x)
        if closed:
            if point.y < y_line:
                return False
        elif point.y <= y_line:
            return False
    if cell.upper_line is not None:
        y_line = y_at_x(lines[cell.upper_line], point.x)
        if closed:
            if point.y > y_line:
                return False
        elif point.y >= y_line:
            return False
    return True


def cell_has_valid_boundary(vd: VerticalDecomposition, cell: VDCell2D) -> bool:
    if cell.left_x is not None and cell.right_x is not None:
        if not (cell.left_x < cell.right_x):
            return False
    lines = vd.arrangement.lines
    sample_x = _sample_x_of(cell)
    if cell.lower_line is not None and cell.upper_line is not None:
        if not (y_at_x(lines[cell.lower_line], sample_x) < y_at_x(lines[cell.upper_line], sample_x)):
            return False
    if cell.lower_line is not None:
        if lines[cell.lower_line].b == 0:
            return False
    if cell.upper_line is not None:
        if lines[cell.upper_line].b == 0:
            return False
    return True


def vertical_walls_at_most_four(cell: VDCell2D) -> bool:
    return len(cell.vertical_walls) <= 4


def wall_belongs_to_decomposition(vd: VerticalDecomposition) -> bool:
    xs = {vertex.point.x for vertex in vd.arrangement.vertices}
    for line in vd.arrangement.lines:
        if line.b == 0:
            xs.add(-line.c / line.a)
    used = {wid for cell in vd.cells for wid in cell.vertical_walls}
    for wall in vd.walls:
        if wall.x not in xs:
            return False
        if wall.y_min is not None and wall.y_max is not None and wall.y_min >= wall.y_max:
            return False
        if wall.id not in used and not _wall_is_degenerate_for_all_cells(vd, wall):
            return False
    return True


def _wall_is_degenerate_for_all_cells(
    vd: VerticalDecomposition, wall
) -> bool:
    """A wall that only meets cells at a point is unused; that is allowed."""
    lines = vd.arrangement.lines
    for cell in vd.cells:
        if cell.left_x != wall.x and cell.right_x != wall.x:
            continue
        lo, hi = y_range_at(cell, wall.x, lines)
        if interval_positive(lo, hi):
            return False
    return True


def interiors_disjoint_at_representatives(vd: VerticalDecomposition) -> bool:
    for cell in vd.cells:
        if not point_in_cell(vd, cell, cell.representative, closed=False):
            return False
        for other in vd.cells:
            if other.id == cell.id:
                continue
            if point_in_cell(vd, other, cell.representative, closed=False):
                return False
    return True


def source_face_matches_representative(vd: VerticalDecomposition) -> bool:
    arr = vd.arrangement
    for cell in vd.cells:
        face = arr.faces[cell.source_face]
        if not point_in_face(arr, face, cell.representative, closed=False):
            return False
        if not point_in_cell(vd, cell, cell.representative, closed=False):
            return False
    return True


def adjacency_symmetric(vd: VerticalDecomposition) -> bool:
    ids = {cell.id for cell in vd.cells}
    for cell in vd.cells:
        for nid in cell.neighbors:
            if nid not in ids:
                return False
            other = vd.cells[nid]
            if cell.id not in other.neighbors:
                return False
    return True


def face_representatives_covered(vd: VerticalDecomposition) -> bool:
    """Every arrangement-face sample is in the VD of that face.

    A face representative may land on a Steiner vertical wall (open count 0).
    Then it must lie in at least one closed cell, all with that source face.
    """
    for face in vd.arrangement.faces:
        open_ids = [
            cell.id
            for cell in vd.cells
            if point_in_cell(vd, cell, face.representative, closed=False)
        ]
        if len(open_ids) == 1:
            if vd.cells[open_ids[0]].source_face != face.id:
                return False
            continue
        if len(open_ids) != 0:
            return False
        closed_ids = [
            cell.id
            for cell in vd.cells
            if point_in_cell(vd, cell, face.representative, closed=True)
        ]
        if not closed_ids:
            return False
        if any(vd.cells[i].source_face != face.id for i in closed_ids):
            return False
    return True


def verify_vd_invariants(vd: VerticalDecomposition) -> None:
    """Raise ``AssertionError`` if any checked invariant fails."""
    checks = (
        ("cell valid boundary", lambda: all(cell_has_valid_boundary(vd, c) for c in vd.cells)),
        ("vertical walls <= 4", lambda: all(vertical_walls_at_most_four(c) for c in vd.cells)),
        ("walls belong to decomposition", lambda: wall_belongs_to_decomposition(vd)),
        ("interiors disjoint at representatives", lambda: interiors_disjoint_at_representatives(vd)),
        ("source face matches representative", lambda: source_face_matches_representative(vd)),
        ("adjacency symmetric", lambda: adjacency_symmetric(vd)),
        ("face representatives covered", lambda: face_representatives_covered(vd)),
    )
    for name, pred in checks:
        if not pred():
            raise AssertionError(f"VD invariant failed: {name}")


def _sample_x_of(cell: VDCell2D) -> Scalar:
    if cell.left_x is None and cell.right_x is None:
        return Scalar(0)
    if cell.left_x is None:
        return cell.right_x - 1  # type: ignore[operator]
    if cell.right_x is None:
        return cell.left_x + 1
    return (cell.left_x + cell.right_x) / 2
