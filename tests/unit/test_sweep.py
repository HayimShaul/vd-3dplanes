"""Phase 8 — reference sweep and 3D cells. Design §13–19, tests 19–22."""

from __future__ import annotations

import random
from pathlib import Path

import pytest

from tests.oracles.events import random_planes_general_position
from tests.oracles.sweep import (
    assert_grid_partition_3d,
    assert_mid_interval_matches_recompute,
    assert_point_location_partition,
    sample_box_points,
)
from vd3d.cells3d import (
    continue_3d_cell,
    end_3d_cell,
    start_3d_cell,
    verify_cell3d,
    verify_cells3d,
    vertical_walls_at_most_four,
)
from vd3d.events import (
    EventType,
    compute_vd_at_z,
    events_have_unique_z,
    generate_all_events,
)
from vd3d.events.samples import (
    planes_alignment_at_z2,
    planes_one,
    planes_through_123,
    planes_two,
)
from vd3d.geometry import Plane
from vd3d.sweep import (
    SimultaneousEvents,
    cell_signature,
    compute_vd_around_event,
    local_cell_ids,
    match_cells,
    process_event,
    same_combinatorics,
    signatures_unique,
    unmatched_after,
    unmatched_before,
    vertical_decomposition_3d,
    z_before_after,
)
from vd3d.sweep.matching import cell_near_event
from vd3d.vertical_decomposition import verify_vd_invariants


def _gp_planes(rng: random.Random, n: int) -> list[Plane]:
    for _ in range(40):
        planes = random_planes_general_position(rng, n)
        events = generate_all_events(planes)
        if events_have_unique_z(events):
            return planes
    raise AssertionError(f"failed to sample {n} general-position planes with unique event z")


# ---------------------------------------------------------------------------
# Step 8.1 — cell matching
# ---------------------------------------------------------------------------


def test_match_far_from_event_is_identity_on_same_slice():
    planes = planes_through_123()
    vd = compute_vd_at_z(planes, 0)
    verify_vd_invariants(vd)
    matches = match_cells(vd, vd)
    assert len(matches) == len(vd.cells)
    assert all(pair.before.id == pair.after.id for pair in matches)
    assert signatures_unique(vd)


def test_triple_fixture_only_local_neighbourhood_unmatched():
    """Test 21: far cells match 1-1; only the triple neighbourhood fails."""
    planes = planes_through_123()
    events = generate_all_events(planes)
    assert len(events) == 1
    event = events[0]
    vd_before, vd_after, z_minus, z_plus = compute_vd_around_event(planes, event, events)
    assert z_minus < event.z < z_plus
    verify_vd_invariants(vd_before)
    verify_vd_invariants(vd_after)
    matches = match_cells(vd_before, vd_after)
    dying = unmatched_before(vd_before, matches)
    born = unmatched_after(vd_after, matches)
    assert len(matches) >= 1
    assert len(dying) >= 1
    assert len(born) >= 1
    local_before = local_cell_ids(vd_before, event)
    local_after = local_cell_ids(vd_after, event)
    assert {cell.id for cell in dying} <= local_before
    assert {cell.id for cell in born} <= local_after
    far_before = [cell for cell in vd_before.cells if cell.id not in local_before]
    matched_before = {pair.before.id for pair in matches}
    assert all(cell.id in matched_before for cell in far_before)


def test_z_before_after_rejects_simultaneous():
    planes = planes_alignment_at_z2()
    events = generate_all_events(planes)
    simultaneous = [event for event in events if event.z == events[0].z]
    assert len(simultaneous) >= 2
    with pytest.raises(SimultaneousEvents):
        z_before_after(simultaneous[0], events)


# ---------------------------------------------------------------------------
# Step 8.2 — 3D cell lifecycle
# ---------------------------------------------------------------------------


def test_start_continue_end_on_two_cell_transition():
    planes = planes_one()
    vd = compute_vd_at_z(planes, 0)
    assert len(vd.cells) == 2
    cells = []
    active = []
    for cell in vd.cells:
        record = start_3d_cell(cells, cell, vd, planes, None)
        active.append(record)
    assert len(cells) == 2
    assert len(active) == len(vd.cells)
    for record, cell in zip(active, vd.cells):
        continue_3d_cell(cells, record, cell, vd, planes)
    assert len(active) == 2
    for record in active:
        ended = end_3d_cell(cells, record, None)
        verify_cell3d(ended)
        assert ended.lower_z is None
        assert ended.upper_z is None
        assert ended.vertical_walls == ()
    floors = {cell.floor.id if cell.floor else None for cell in cells}
    ceilings = {cell.ceiling.id if cell.ceiling else None for cell in cells}
    assert floors == {1, None}
    assert ceilings == {1, None}
    verify_cells3d(cells)


def test_lifecycle_triple_keeps_active_count():
    planes = planes_through_123()
    events = generate_all_events(planes)
    from vd3d.events import choose_z_below_all_events
    from vd3d.sweep.matching import cell_signature as sig
    from vd3d.cells3d.types import ActiveCell

    z0 = choose_z_below_all_events(events)
    vd0 = compute_vd_at_z(planes, z0)
    cells = []
    active = {sig(vd0, cell): start_3d_cell(cells, cell, vd0, planes, None) for cell in vd0.cells}
    assert len(active) == len(vd0.cells)
    _vd, new_active, snapshot = process_event(planes, events[0], events, cells, active, vd0)
    assert snapshot.n_active == snapshot.n_after
    assert len(new_active) == snapshot.n_after
    assert isinstance(next(iter(new_active.values())), ActiveCell)


# ---------------------------------------------------------------------------
# Step 8.3 — PROCESS_EVENT snapshots
# ---------------------------------------------------------------------------


def test_process_event_writes_snapshot_json(tmp_path: Path):
    planes = planes_through_123()
    events = generate_all_events(planes)
    from vd3d.events import choose_z_below_all_events
    from vd3d.sweep.matching import cell_signature as sig

    z0 = choose_z_below_all_events(events)
    vd0 = compute_vd_at_z(planes, z0)
    cells = []
    active = {sig(vd0, cell): start_3d_cell(cells, cell, vd0, planes, None) for cell in vd0.cells}
    _vd, new_active, snapshot = process_event(
        planes, events[0], events, cells, active, vd0, snapshot_dir=tmp_path
    )
    files = list(tmp_path.glob("*.json"))
    assert len(files) == 1
    text = files[0].read_text(encoding="utf-8")
    assert events[0].stable_id.split(":")[0] in text or "TRIPLE" in text
    assert '"n_active"' in text
    assert len(new_active) == snapshot.n_after


def test_full_sweep_snapshots_and_active_count():
    planes = planes_through_123()
    result = vertical_decomposition_3d(planes)
    assert len(result.snapshots) == 1
    snap = result.snapshots[0]
    assert snap.n_active == snap.n_after
    assert snap.ended
    assert snap.started
    verify_cells3d(result.cells)


# ---------------------------------------------------------------------------
# Step 8.4 — main algorithm, tests 19–22
# ---------------------------------------------------------------------------


def test_19_one_plane_two_cells_no_walls():
    planes = planes_one()
    result = vertical_decomposition_3d(planes)
    assert result.events == ()
    assert len(result.cells) == 2
    for cell in result.cells:
        assert cell.vertical_walls == ()
        assert vertical_walls_at_most_four(cell)
        assert cell.lower_z is None
        assert cell.upper_z is None
    sides = {
        (cell.floor.id if cell.floor else None, cell.ceiling.id if cell.ceiling else None)
        for cell in result.cells
    }
    assert sides == {(1, None), (None, 1)}


def test_20_two_planes():
    planes = planes_two()
    result = vertical_decomposition_3d(planes)
    assert result.events == ()
    assert len(result.cells) == 4
    for cell in result.cells:
        verify_cell3d(cell)
        assert len(cell.vertical_walls) <= 1
        assert cell.lower_z is None
        assert cell.upper_z is None
    floors = {cell.floor.id for cell in result.cells if cell.floor is not None}
    ceilings = {cell.ceiling.id for cell in result.cells if cell.ceiling is not None}
    assert floors == {1, 2} or floors <= {1, 2}
    assert ceilings == {1, 2} or ceilings <= {1, 2}
    with_floor = sum(cell.floor is not None for cell in result.cells)
    with_ceiling = sum(cell.ceiling is not None for cell in result.cells)
    assert with_floor == 2
    assert with_ceiling == 2


def test_21_three_planes_local_change_only():
    planes = planes_through_123()
    result = vertical_decomposition_3d(planes)
    assert len(result.events) == 1
    assert result.events[0].type is EventType.TRIPLE_INTERSECTION
    snap = result.snapshots[0]
    assert snap.n_before == snap.n_after
    assert snap.ended and snap.started
    assert snap.n_active == snap.n_after
    verify_cells3d(result.cells)
    ended_ids = set(snap.ended)
    started_ids = set(snap.started)
    vd_before, vd_after, *_ = compute_vd_around_event(
        planes, result.events[0], result.events
    )
    local_before = local_cell_ids(vd_before, result.events[0])
    local_after = local_cell_ids(vd_after, result.events[0])
    assert ended_ids <= local_before
    assert started_ids <= local_after


def test_22_alignment_changes_locally():
    """Alignment at z=2: combinatorics change locally; far cells keep signatures."""
    planes = planes_alignment_at_z2()
    events = generate_all_events(planes)
    alignments = [event for event in events if event.type is EventType.VERTICAL_ALIGNMENT]
    assert len(alignments) == 1
    event = alignments[0]
    assert event.z == 2
    # Isolate this event so z±ε does not collide with the simultaneous triples.
    vd_before, vd_after, *_ = compute_vd_around_event(planes, event, (event,))
    verify_vd_invariants(vd_before)
    verify_vd_invariants(vd_after)
    assert not same_combinatorics(vd_before, vd_after)
    matches = match_cells(vd_before, vd_after)
    changed = [
        pair
        for pair in matches
        if cell_signature(vd_before, pair.before) != cell_signature(vd_after, pair.after)
    ]
    unmatched_b = unmatched_before(vd_before, matches)
    unmatched_a = unmatched_after(vd_after, matches)
    local_before = local_cell_ids(vd_before, event)
    local_after = local_cell_ids(vd_after, event)
    assert changed or unmatched_b or unmatched_a
    for pair in changed:
        assert pair.before.id in local_before or cell_near_event(
            vd_before, pair.before, event
        )
    assert {cell.id for cell in unmatched_b} <= local_before | {
        cell.id for cell in vd_before.cells if cell_near_event(vd_before, cell, event)
    }
    assert {cell.id for cell in unmatched_a} <= local_after | {
        cell.id for cell in vd_after.cells if cell_near_event(vd_after, cell, event)
    }
    far_before = [cell for cell in vd_before.cells if cell.id not in local_before]
    matched_before = {pair.before.id for pair in matches}
    for cell in far_before:
        assert cell.id in matched_before


def test_simultaneous_events_rejected_on_full_sweep():
    with pytest.raises(SimultaneousEvents):
        vertical_decomposition_3d(planes_alignment_at_z2())


# ---------------------------------------------------------------------------
# Step 8.5 — global oracles
# ---------------------------------------------------------------------------


def test_oracles_one_two_three_planes():
    for planes in (planes_one(), planes_two(), planes_through_123()):
        result = vertical_decomposition_3d(planes)
        verify_cells3d(result.cells)
        assert_mid_interval_matches_recompute(result)
        assert_grid_partition_3d(result, half=2, step=1)
        rng = random.Random(81)
        assert_point_location_partition(result, sample_box_points(rng, n=30, half=3))


def test_oracles_random_gp_n_le_5():
    rng = random.Random(82)
    seen = 0
    for n in (3, 4, 5):
        for _ in range(3):
            planes = _gp_planes(rng, n)
            result = vertical_decomposition_3d(planes)
            verify_cells3d(result.cells)
            assert_mid_interval_matches_recompute(result)
            assert_grid_partition_3d(result, half=2, step=1)
            assert_point_location_partition(
                result, sample_box_points(rng, n=20, half=3)
            )
            seen += 1
    assert seen == 9
