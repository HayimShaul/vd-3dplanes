"""Brute-force partition checks for a 2D vertical decomposition."""

from __future__ import annotations

from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar
from vd3d.vertical_decomposition.invariants import point_in_cell
from vd3d.vertical_decomposition.types import VerticalDecomposition


def cells_containing(
    vd: VerticalDecomposition, point: Point2D, *, closed: bool = False
) -> list[int]:
    return [
        cell.id
        for cell in vd.cells
        if point_in_cell(vd, cell, point, closed=closed)
    ]


def point_on_vd_boundary(vd: VerticalDecomposition, point: Point2D) -> bool:
    """True iff ``point`` lies on a supporting line or a vertical wall line."""
    for line in vd.arrangement.lines:
        if line.contains(point):
            return True
    for wall in vd.walls:
        if point.x != wall.x:
            continue
        if wall.y_min is not None and point.y < wall.y_min:
            continue
        if wall.y_max is not None and point.y > wall.y_max:
            continue
        return True
    return False


def assert_point_in_exactly_one_open_cell(
    vd: VerticalDecomposition, point: Point2D
) -> None:
    """Interior points sit in one cell; boundary points sit in none (open)."""
    open_ids = cells_containing(vd, point, closed=False)
    if point_on_vd_boundary(vd, point):
        if open_ids:
            raise AssertionError(
                f"boundary point {point} is in open cells {open_ids}"
            )
        closed_ids = cells_containing(vd, point, closed=True)
        if not closed_ids:
            raise AssertionError(f"boundary point {point} is in no closed cell")
        return
    if len(open_ids) != 1:
        raise AssertionError(
            f"point {point} is in {len(open_ids)} open cells {open_ids}"
        )


def sample_grid(half: int, *, step: Scalar | int = 1) -> list[Point2D]:
    """Axis-aligned grid of exact points in ``[-half, half]^2``."""
    step_s = step if isinstance(step, Scalar) else Scalar(step)
    points: list[Point2D] = []
    n = int(Scalar(2 * half) / step_s)
    x = Scalar(-half)
    for _ in range(n + 1):
        y = Scalar(-half)
        for _ in range(n + 1):
            points.append(Point2D(x, y))
            y += step_s
        x += step_s
    return points


def assert_grid_partition(
    vd: VerticalDecomposition, *, half: int = 6, step: Scalar | int = 1
) -> None:
    for point in sample_grid(half, step=step):
        assert_point_in_exactly_one_open_cell(vd, point)
    for cell in vd.cells:
        assert_point_in_exactly_one_open_cell(vd, cell.representative)
