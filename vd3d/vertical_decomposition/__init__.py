"""2D vertical decomposition of a line arrangement.

Must not import ``vd3d.sweep`` or ``vd3d.cells3d``.
"""

from vd3d.vertical_decomposition.cells import (
    collect_decomposition_walls,
    compute_vertical_decomposition,
    finalize_decomposition_walls,
    insert_decomposition_segment,
)
from vd3d.vertical_decomposition.geom import (
    line_is_vertical,
    x_of_vertical_line,
    y_at_x,
)
from vd3d.vertical_decomposition.invariants import (
    point_in_cell,
    verify_vd_invariants,
)
from vd3d.vertical_decomposition.rays import first_hit, vertical_rays_from
from vd3d.vertical_decomposition.types import (
    UNBOUNDED,
    Hit,
    VDCell2D,
    VerticalDecomposition,
    VerticalRay,
    VerticalWall,
)

__all__ = [
    "UNBOUNDED",
    "Hit",
    "VDCell2D",
    "VerticalDecomposition",
    "VerticalRay",
    "VerticalWall",
    "collect_decomposition_walls",
    "compute_vertical_decomposition",
    "finalize_decomposition_walls",
    "first_hit",
    "insert_decomposition_segment",
    "line_is_vertical",
    "point_in_cell",
    "verify_vd_invariants",
    "vertical_rays_from",
    "x_of_vertical_line",
    "y_at_x",
]
