"""Zone of a query line in a 2D arrangement.

Must not import ``vd3d.sweep`` or ``vd3d.cells3d``.
"""

from vd3d.zone.crossings import compute_crossings
from vd3d.zone.geom import (
    edge_contains_point,
    point_on_line_at_parameter,
    vertex_id_at,
)
from vd3d.zone.invariants import verify_supporting_line_zone, verify_zone_invariants
from vd3d.zone.supporting import compute_supporting_line_zone
from vd3d.zone.types import (
    Crossing,
    PointNotInOpenFace,
    QueryOverlapsArrangement,
    QueryThroughVertex,
    SupportingLineZone,
    Zone,
    ZoneError,
)
from vd3d.zone.walk import compute_zone, face_containing_point, face_on_other_side

__all__ = [
    "Crossing",
    "PointNotInOpenFace",
    "QueryOverlapsArrangement",
    "QueryThroughVertex",
    "SupportingLineZone",
    "Zone",
    "ZoneError",
    "compute_crossings",
    "compute_supporting_line_zone",
    "compute_zone",
    "edge_contains_point",
    "face_containing_point",
    "face_on_other_side",
    "point_on_line_at_parameter",
    "verify_supporting_line_zone",
    "verify_zone_invariants",
    "vertex_id_at",
]
