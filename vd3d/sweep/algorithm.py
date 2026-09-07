"""3D vertical decomposition: incremental 2D updates, recompute as oracle (design §18–19)."""

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
from vd3d.sweep.process import process_event_group
from vd3d.sweep.types import CellSignature, EventGroup, SimultaneousEvents, SweepSnapshot, ZInterval
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


def group_events_by_z(events: Sequence[Event]) -> tuple[EventGroup, ...]:
    """Consecutive events that share a ``z``. General position: each group has one."""
    groups: list[EventGroup] = []
    for event in events:
        if groups and groups[-1].z == event.z:
            groups[-1] = EventGroup(z=event.z, events=(*groups[-1].events, event))
        else:
            groups.append(EventGroup(z=event.z, events=(event,)))
    return tuple(groups)


def vertical_decomposition_3d(
    planes: Sequence[Plane],
    *,
    snapshot_dir: Path | None = None,
    require_general_position: bool = False,
    incremental: bool = True,
) -> SweepResult:
    """``VERTICAL_DECOMPOSITION_3D``. After each event group the 2D VD is
    updated incrementally and checked against ``compute_vd_at_z(z+)``.

    Events that share a ``z`` are one transaction. Pass
    ``require_general_position=True`` to reject those groups.
    """
    require_unique_plane_ids(planes)
    plane_tuple = tuple(planes)
    events = generate_all_events(plane_tuple)
    groups = group_events_by_z(events)
    if require_general_position:
        for group in groups:
            if len(group.events) != 1:
                raise SimultaneousEvents(
                    f"{len(group.events)} events share z={group.z}; "
                    "pass require_general_position=False to process the group"
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
        current_vd, active, snapshot = process_event_group(
            plane_tuple,
            group.events,
            events,
            cells,
            active,
            current_vd,
            snapshot_dir=snapshot_dir,
            incremental=incremental,
        )
        snapshots.append(snapshot)
        intervals.append(
            ZInterval(
                lower_z=group.z,
                upper_z=_next_z(events, group.z),
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


def _next_z(events: Sequence[Event], z):
    remaining = [other.z for other in events if other.z > z]
    return None if not remaining else min(remaining)
