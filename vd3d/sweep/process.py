"""Reference ``PROCESS_EVENT``: recompute 2D VD, match, update 3D cells (design §17, §19)."""

from __future__ import annotations

import json
from collections.abc import Sequence
from pathlib import Path

from vd3d.cells3d.lifecycle import continue_3d_cell, end_3d_cell, start_3d_cell
from vd3d.cells3d.types import ActiveCell, Cell3D
from vd3d.events.types import Event
from vd3d.geometry.plane import Plane
from vd3d.sweep.around import compute_vd_around_event
from vd3d.sweep.matching import (
    cell_signature,
    match_cells,
    same_combinatorics,
    signatures_unique,
    unmatched_after,
    unmatched_before,
)
from vd3d.sweep.types import CellSignature, SweepSnapshot
from vd3d.vertical_decomposition.types import VerticalDecomposition


def process_event(
    planes: Sequence[Plane],
    event: Event,
    events: Sequence[Event],
    cells: list[Cell3D],
    active: dict[CellSignature, ActiveCell],
    current_vd: VerticalDecomposition | None = None,
    *,
    snapshot_dir: Path | None = None,
) -> tuple[VerticalDecomposition, dict[CellSignature, ActiveCell], SweepSnapshot]:
    """Recompute ``vd_before`` / ``vd_after``, match cells, update 3D lifecycle."""
    vd_before, vd_after, _z_minus, _z_plus = compute_vd_around_event(
        planes, event, events
    )
    if not signatures_unique(vd_before) or not signatures_unique(vd_after):
        raise RuntimeError("2D cell signatures are not unique at this event")
    if current_vd is not None and not same_combinatorics(current_vd, vd_before):
        raise AssertionError(
            "current 2D VD is not combinatorially equal to the recomputed z− slice"
        )

    matches = match_cells(vd_before, vd_after)
    dying = unmatched_before(vd_before, matches)
    born = unmatched_after(vd_after, matches)

    for cell in dying:
        sig = cell_signature(vd_before, cell)
        end_3d_cell(cells, active[sig], event.z)

    new_active: dict[CellSignature, ActiveCell] = {}
    continued: list[tuple[int, int]] = []
    for pair in matches:
        sig_before = cell_signature(vd_before, pair.before)
        record = active[sig_before]
        continue_3d_cell(cells, record, pair.after, vd_after, planes)
        new_active[cell_signature(vd_after, pair.after)] = record
        continued.append((pair.before.id, pair.after.id))

    started_ids: list[int] = []
    for cell in born:
        record = start_3d_cell(cells, cell, vd_after, planes, event.z)
        new_active[cell_signature(vd_after, cell)] = record
        started_ids.append(cell.id)

    if len(new_active) != len(vd_after.cells):
        raise AssertionError(
            f"|active 3D cells|={len(new_active)} != |2D cells|={len(vd_after.cells)}"
        )

    snapshot = SweepSnapshot(
        event=event,
        n_before=len(vd_before.cells),
        n_after=len(vd_after.cells),
        matches=tuple(continued),
        started=tuple(started_ids),
        ended=tuple(cell.id for cell in dying),
        n_active=len(new_active),
    )
    if snapshot_dir is not None:
        write_snapshot_json(snapshot, snapshot_dir)
    return vd_after, new_active, snapshot


def write_snapshot_json(snapshot: SweepSnapshot, directory: Path) -> Path:
    """Write ``<directory>/<event_id>.json``. PNG companions live in viz."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{_safe_name(snapshot.event.stable_id)}.json"
    path.write_text(json.dumps(snapshot.as_dict(), indent=2) + "\n", encoding="utf-8")
    return path


def _safe_name(stable_id: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-._" else "_" for ch in stable_id)
