"""Pairwise intersection lines and sweep events.

Phase 5: intersection lines and triple events.
Phase 6: y-parallel walls and vertical-alignment events.
Phase 7: combined event list and the 2D VD at a given ``z``.
Phase 8 lives in ``vd3d.sweep`` / ``vd3d.cells3d``.
Phase 9 incremental 2D updates live in ``vd3d.sweep.update``.
Must not import ``vd3d.viz``.
"""

from vd3d.events.all import choose_z_below_all_events, generate_all_events
from vd3d.events.alignment import (
    WallLineZone,
    alignment_event_for_opposite_vertex,
    alignment_visible,
    canonical_event_key,
    compute_line_zone_on_wall,
    deduplicate_events,
    enumerate_alignment_pairs,
    extract_alignment_events_from_zone,
    generate_alignment_events,
    make_alignment_event,
    plane_on_wall,
    represent_line_on_wall,
    slice_planes_by_wall,
    validate_alignment_event,
)
from vd3d.events.invariants import (
    events_have_unique_z,
    verify_alignment_events,
    verify_all_events,
    verify_intersection_lines,
    verify_triple_events,
    z_strictly_below_all_events,
)
from vd3d.events.lines import compute_intersection_lines
from vd3d.events.slice import build_2d_arrangement_at_z, compute_vd_at_z, initial_slice
from vd3d.events.triples import generate_triple_events
from vd3d.events.types import AlignmentGeometry, Event, EventType, default_event_stable_id
from vd3d.events.wall import (
    VerticalIntersectionLine,
    alignment_z,
    build_vertical_wall,
    lift_wall_point,
    line_xz_param,
    lines_meet,
    point_on_line_at_z,
    wall_frame_point,
    x_of_line,
    y_of_line,
)

__all__ = [
    "AlignmentGeometry",
    "Event",
    "EventType",
    "VerticalIntersectionLine",
    "WallLineZone",
    "alignment_event_for_opposite_vertex",
    "alignment_visible",
    "alignment_z",
    "build_2d_arrangement_at_z",
    "build_vertical_wall",
    "canonical_event_key",
    "choose_z_below_all_events",
    "compute_intersection_lines",
    "compute_vd_at_z",
    "compute_line_zone_on_wall",
    "deduplicate_events",
    "default_event_stable_id",
    "enumerate_alignment_pairs",
    "events_have_unique_z",
    "extract_alignment_events_from_zone",
    "generate_alignment_events",
    "generate_all_events",
    "generate_triple_events",
    "initial_slice",
    "lift_wall_point",
    "line_xz_param",
    "lines_meet",
    "make_alignment_event",
    "plane_on_wall",
    "point_on_line_at_z",
    "represent_line_on_wall",
    "slice_planes_by_wall",
    "validate_alignment_event",
    "verify_alignment_events",
    "verify_all_events",
    "verify_intersection_lines",
    "verify_triple_events",
    "z_strictly_below_all_events",
    "wall_frame_point",
    "x_of_line",
    "y_of_line",
]
