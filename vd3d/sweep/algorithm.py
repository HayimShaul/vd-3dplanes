"""Reference 3D vertical decomposition: recompute 2D VD at every event (design §18–19)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from vd3d.cells3d.invariants import verify_cells3d
from vd3d.cells3d.lifecycle import end_3d_cell, start_3d_cell
from vd3d.cells3d.types import ActiveCell, Cell3D
from vd3d.events.all import choose_z_below_all_events, generate_all_events
from vd3d.events.ids import require_unique_plane_ids
from vd3d.events.slice import compute_vd_at_z
from vd3d.events.types import Event
from vd3d.geometry.plane import Plane
from vd3d.sweep.matching import cell_signature, signatures_unique
from vd3d.sweep.process import process_event
from vd3d.sweep.types import CellSignature, SimultaneousEvents, SweepSnapshot, ZInterval
from vd3d.vertical_decomposition.types import VerticalDecomposition


@dataclass(slots=True)
class SweepResult:
    """Completed 3D cells plus the debugging trail of the reference sweep."""

    planes: tuple[Plane, ...]
    cells: tuple[Cell3D, ...]
    events: tuple[Event, ...]
    snapshots: tuple[SweepSnapshot, ...]
    intervals: tuple[ZInterval, ...]
    initial_vd: VerticalDecomposition


def group_events_by_z(events: Sequence[Event]) -> tuple[tuple[Event, ...], ...]:
    """Consecutive events that share a ``z``. General position: each group has one."""
    groups: list[tuple[Event, ...]] = []
    for event in events:
        if groups and groups[-1][0].z == event.z:
            groups[-1] = (*groups[-1], event)
        else:
            groups.append((event,))
    return tuple(groups)


def vertical_decomposition_3d(
    planes: Sequence[Plane],
    *,
    snapshot_dir: Path | None = None,
    require_general_position: bool = True,
) -> SweepResult:
    """``VERTICAL_DECOMPOSITION_3D`` with recomputed 2D slices (no incremental update)."""
    require_unique_plane_ids(planes)
    plane_tuple = tuple(planes)
    events = generate_all_events(plane_tuple)
    groups = group_events_by_z(events)
    if require_general_position:
        for group in groups:
            if len(group) != 1:
                raise SimultaneousEvents(
                    f"{len(group)} events share z={group[0].z}; "
                    "simultaneous groups are Phase 10"
                )

    z0 = choose_z_below_all_events(events)
    initial_vd = compute_vd_at_z(plane_tuple, z0)
    if not signatures_unique(initial_vd):
        raise RuntimeError("initial 2D cell signatures are not unique")

    cells: list[Cell3D] = []
    active: dict[CellSignature, ActiveCell] = {}
    for cell2d in initial_vd.cells:
        record = start_3d_cell(cells, cell2d, initial_vd, plane_tuple, None)
        active[cell_signature(initial_vd, cell2d)] = record
    if len(active) != len(initial_vd.cells):
        raise AssertionError("|active 3D cells| != |initial 2D cells|")

    snapshots: list[SweepSnapshot] = []
    intervals: list[ZInterval] = []
    current_vd = initial_vd

    def _bind(vd: VerticalDecomposition) -> tuple[tuple[CellSignature, int], ...]:
        return tuple(
            (cell_signature(vd, cell), active[cell_signature(vd, cell)].cell3d_id)
            for cell in vd.cells
        )

    intervals.append(
        ZInterval(lower_z=None, upper_z=_first_z(events), binding=_bind(current_vd))
    )

    for group in groups:
        event = group[0]
        current_vd, active, snapshot = process_event(
            plane_tuple,
            event,
            events,
            cells,
            active,
            current_vd,
            snapshot_dir=snapshot_dir,
        )
        snapshots.append(snapshot)
        intervals.append(
            ZInterval(
                lower_z=event.z,
                upper_z=_next_z(events, event),
                binding=_bind(current_vd),
            )
        )

    for record in list(active.values()):
        end_3d_cell(cells, record, None)

    verify_cells3d(cells)
    return SweepResult(
        planes=plane_tuple,
        cells=tuple(cells),
        events=events,
        snapshots=tuple(snapshots),
        intervals=tuple(intervals),
        initial_vd=initial_vd,
    )


def _first_z(events: Sequence[Event]):
    return None if not events else events[0].z


def _next_z(events: Sequence[Event], event: Event):
    remaining = [other.z for other in events if other.z > event.z]
    return None if not remaining else min(remaining)
