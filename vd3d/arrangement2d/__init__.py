"""2D line arrangement (DCEL). Independent of the 3D sweep.

Must not import ``vd3d.sweep`` or ``vd3d.cells3d``.
"""

from vd3d.arrangement2d.build import build_line_arrangement
from vd3d.arrangement2d.intersect import (
    COINCIDENT,
    PARALLEL,
    CoincidentLines,
    ParallelLines,
    intersect_lines_2d,
    lines_parallel,
)
from vd3d.arrangement2d.invariants import (
    euler_characteristic_plane,
    point_in_face,
    verify_arrangement_invariants,
)
from vd3d.arrangement2d.pieces import compute_edge_pieces, compute_vertices
from vd3d.arrangement2d.predicates import (
    cmp_direction_ccw,
    cross2,
    line_direction,
    parameter_on_line,
    sort_directions_ccw,
)
from vd3d.arrangement2d.types import Arrangement2D, EdgePiece, Face, HalfEdge, Vertex

__all__ = [
    "COINCIDENT",
    "PARALLEL",
    "Arrangement2D",
    "CoincidentLines",
    "EdgePiece",
    "Face",
    "HalfEdge",
    "ParallelLines",
    "Vertex",
    "build_line_arrangement",
    "cmp_direction_ccw",
    "compute_edge_pieces",
    "compute_vertices",
    "cross2",
    "euler_characteristic_plane",
    "intersect_lines_2d",
    "line_direction",
    "lines_parallel",
    "parameter_on_line",
    "point_in_face",
    "sort_directions_ccw",
    "verify_arrangement_invariants",
]
