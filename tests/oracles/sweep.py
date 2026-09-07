"""Brute-force checkers for the reference 3D sweep."""

from __future__ import annotations

import random
from collections.abc import Sequence

from vd3d.events.slice import compute_vd_at_z
from vd3d.geometry.points import Point2D, Point3D
from vd3d.geometry.scalar import Scalar, as_scalar
from vd3d.sweep.algorithm import SweepResult
from vd3d.sweep.locate import interval_containing, locate_cell3d
from vd3d.sweep.matching import cell_signature
from vd3d.vertical_decomposition.invariants import point_in_cell


def mid_interval_z(lower: Scalar | None, upper: Scalar | None) -> Scalar:
    """A sample height strictly inside an open z-interval."""
    if lower is None and upper is None:
        return as_scalar(0)
    if lower is None:
        return upper - 1  # type: ignore[operator]
    if upper is None:
        return lower + 1
    return (lower + upper) / 2


def assert_mid_interval_matches_recompute(result: SweepResult) -> None:
    """At a sample ``z`` in each interval, signatures equal ``compute_vd_at_z``."""
    for interval in result.intervals:
        z = mid_interval_z(interval.lower_z, interval.upper_z)
        vd = compute_vd_at_z(result.planes, z)
        keys = {cell_signature(vd, cell) for cell in vd.cells}
        bound = {sig for sig, _ in interval.binding}
        if keys != bound:
            raise AssertionError(
                f"interval ({interval.lower_z}, {interval.upper_z}) signatures "
                f"differ from compute_vd_at_z({z})"
            )
        if len(interval.binding) != len(vd.cells):
            raise AssertionError("interval binding is not 1-1 with 2D cells")


def sample_box_points(
    rng: random.Random,
    *,
    n: int,
    half: int = 3,
) -> list[Point3D]:
    """Exact random points in ``[-half, half]^3``."""
    points: list[Point3D] = []
    for _ in range(n):
        points.append(
            Point3D(
                rng.randint(-half, half),
                rng.randint(-half, half),
                rng.randint(-half, half),
            )
        )
    return points


def assert_point_location_partition(
    result: SweepResult, points: Sequence[Point3D]
) -> None:
    """Each interior sample is in exactly one 3D cell; boundary samples in none."""
    for point in points:
        cell_id = locate_cell3d(result, point)
        interval = interval_containing(result, point.z)
        try:
            vd = compute_vd_at_z(result.planes, point.z)
        except (ValueError, RuntimeError):
            if cell_id is not None:
                raise AssertionError(
                    f"degenerate-z point {point} located in 3D cell {cell_id}"
                )
            continue
        xy = Point2D(point.x, point.y)
        open_2d = [
            cell
            for cell in vd.cells
            if point_in_cell(vd, cell, xy, closed=False)
        ]
        on_boundary = interval is None or len(open_2d) != 1
        if on_boundary:
            if cell_id is not None:
                raise AssertionError(
                    f"boundary point {point} located in 3D cell {cell_id}"
                )
            continue
        if cell_id is None:
            sig = cell_signature(vd, open_2d[0])
            bound = {s for s, _ in interval.binding}
            raise AssertionError(
                f"interior point {point} is in no 3D cell; "
                f"interval=({interval.lower_z},{interval.upper_z}) "
                f"sig={sig!r} in_binding={sig in bound} "
                f"nbind={len(bound)} n2d={len(vd.cells)}"
            )
        expected = dict(interval.binding)[cell_signature(vd, open_2d[0])]
        if cell_id != expected:
            raise AssertionError(
                f"point {point} located in {cell_id}, expected {expected}"
            )


def assert_grid_partition_3d(
    result: SweepResult, *, half: int = 2, step: int = 1
) -> None:
    """Axis-aligned exact grid in ``[-half, half]^3`` is a partition."""
    points: list[Point3D] = []
    z = as_scalar(-half)
    while z <= half:
        y = as_scalar(-half)
        while y <= half:
            x = as_scalar(-half)
            while x <= half:
                points.append(Point3D(x, y, z))
                x += step
            y += step
        z += step
    assert_point_location_partition(result, points)
