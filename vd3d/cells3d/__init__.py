"""3D cells grown during the sweep.

Must not import ``vd3d.viz``. Must not import ``vd3d.sweep`` (the sweep
calls into this package).
"""

from vd3d.cells3d.invariants import (
    ceiling_at_most_one,
    floor_at_most_one,
    verify_cell3d,
    verify_cells3d,
    vertical_walls_at_most_four,
    z_extent_ordered,
)
from vd3d.cells3d.lifecycle import continue_3d_cell, end_3d_cell, start_3d_cell
from vd3d.cells3d.types import ActiveCell, Cell3D
from vd3d.cells3d.walls import (
    extract_vertical_walls,
    merge_planes,
    planes_at_x,
    side_planes,
    supporting_plane_ids,
    supporting_planes,
)

__all__ = [
    "ActiveCell",
    "Cell3D",
    "ceiling_at_most_one",
    "continue_3d_cell",
    "end_3d_cell",
    "extract_vertical_walls",
    "floor_at_most_one",
    "merge_planes",
    "planes_at_x",
    "side_planes",
    "start_3d_cell",
    "supporting_plane_ids",
    "supporting_planes",
    "verify_cell3d",
    "verify_cells3d",
    "vertical_walls_at_most_four",
    "z_extent_ordered",
]
