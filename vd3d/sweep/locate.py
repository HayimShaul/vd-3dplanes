"""Point location in the reference 3D decomposition (design §30)."""

from __future__ import annotations

from vd3d.events.slice import compute_vd_at_z
from vd3d.geometry.points import Point2D, Point3D
from vd3d.geometry.scalar import Scalar
from vd3d.sweep.algorithm import SweepResult
from vd3d.sweep.matching import cell_signature, locate_cell
from vd3d.sweep.types import ZInterval


def interval_containing(result: SweepResult, z: Scalar) -> ZInterval | None:
    """The open interval with ``lower < z < upper``, or ``None`` on an event."""
    for interval in result.intervals:
        if interval.lower_z is not None and z <= interval.lower_z:
            continue
        if interval.upper_z is not None and z >= interval.upper_z:
            continue
        return interval
    return None


def locate_cell3d(result: SweepResult, point: Point3D) -> int | None:
    """3D cell id containing ``point``, or ``None`` if it is on a boundary.

    Independent path: 2D VD at ``point.z``, then the interval's signature map.
    """
    interval = interval_containing(result, point.z)
    if interval is None:
        return None
    try:
        vd = compute_vd_at_z(result.planes, point.z)
    except (ValueError, RuntimeError):
        return None
    cell = locate_cell(vd, Point2D(point.x, point.y), closed=False)
    if cell is None:
        return None
    sig = cell_signature(vd, cell)
    lookup = interval.lookup()
    return lookup.get(sig)
