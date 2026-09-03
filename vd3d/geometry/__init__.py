"""Exact geometric kernel.

2D and 3D types live here. This package must not import ``vd3d.viz``
or any floating-point plotting library.
"""

from vd3d.geometry.intersections import (
    PARALLEL,
    ParallelPlanes,
    build_slice_lines,
    intersect_planes,
    intersect_three_planes,
    slice_plane_at_z,
)
from vd3d.geometry.invariants import (
    line_lies_on_both_planes,
    line_lies_on_plane,
    point_lies_on_planes,
    slice_lifts_to_plane,
)
from vd3d.geometry.linalg import cross, det2, det3, solve2, solve3
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.line3d import Line3D, directions_parallel
from vd3d.geometry.plane import Plane, normals_parallel
from vd3d.geometry.points import Point2D, Point3D
from vd3d.geometry.scalar import Scalar, as_scalar

__all__ = [
    "Line2D",
    "Line3D",
    "PARALLEL",
    "ParallelPlanes",
    "Plane",
    "Point2D",
    "Point3D",
    "Scalar",
    "as_scalar",
    "build_slice_lines",
    "cross",
    "det2",
    "det3",
    "directions_parallel",
    "intersect_planes",
    "intersect_three_planes",
    "line_lies_on_both_planes",
    "line_lies_on_plane",
    "normals_parallel",
    "point_lies_on_planes",
    "slice_lifts_to_plane",
    "slice_plane_at_z",
    "solve2",
    "solve3",
]
