"""z-sweep of the 2D vertical decomposition.

Reference implementation: incremental 2D updates at every event, checked
against a recomputed ``z+`` slice (design §12, §19).

Must not import ``vd3d.viz``. 2D packages must not import this package.
"""

from vd3d.sweep.algorithm import (
    SweepResult,
    group_events_by_z,
    vertical_decomposition_3d,
)
from vd3d.sweep.around import (
    DEFAULT_EVENT_EPS,
    compute_vd_around_event,
    z_before_after,
)
from vd3d.sweep.invariants import (
    active_count_matches_2d,
    equivalent_vd,
    verify_active_against_vd,
    verify_incremental_matches_recompute,
    verify_interval_matches_vd,
)
from vd3d.sweep.locate import interval_containing, locate_cell3d
from vd3d.sweep.matching import (
    cell_near_event,
    cell_signature,
    event_points_2d,
    local_cell_ids,
    local_cell_ids_for_events,
    locate_cell,
    match_cells,
    same_combinatorics,
    signatures_unique,
    unmatched_after,
    unmatched_before,
)
from vd3d.sweep.process import process_event, process_event_group, write_snapshot_json
from vd3d.sweep.types import (
    CellMatch,
    CellSignature,
    EventGroup,
    SimultaneousEvents,
    SweepSnapshot,
    ZInterval,
)
from vd3d.sweep.update import (
    event_vertex_plane_keys,
    local_vertex_plane_keys,
    local_vertex_plane_keys_for_events,
    update_2d_decomposition,
    update_2d_for_event_group,
    update_for_triple_intersection,
    update_for_vertical_alignment,
    vertex_plane_key,
)

__all__ = [
    "DEFAULT_EVENT_EPS",
    "CellMatch",
    "CellSignature",
    "EventGroup",
    "SimultaneousEvents",
    "SweepResult",
    "SweepSnapshot",
    "ZInterval",
    "active_count_matches_2d",
    "cell_near_event",
    "cell_signature",
    "compute_vd_around_event",
    "equivalent_vd",
    "event_points_2d",
    "event_vertex_plane_keys",
    "group_events_by_z",
    "interval_containing",
    "local_cell_ids",
    "local_cell_ids_for_events",
    "local_vertex_plane_keys",
    "local_vertex_plane_keys_for_events",
    "locate_cell",
    "locate_cell3d",
    "match_cells",
    "process_event",
    "process_event_group",
    "same_combinatorics",
    "signatures_unique",
    "unmatched_after",
    "unmatched_before",
    "update_2d_decomposition",
    "update_2d_for_event_group",
    "update_for_triple_intersection",
    "update_for_vertical_alignment",
    "verify_active_against_vd",
    "verify_incremental_matches_recompute",
    "verify_interval_matches_vd",
    "vertex_plane_key",
    "vertical_decomposition_3d",
    "write_snapshot_json",
    "z_before_after",
]
