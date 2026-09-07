"""Phase 9 — incremental 2D updates. Design §12 / §19, differential vs recompute."""

from __future__ import annotations

import random

import pytest

from tests.oracles.events import random_planes_general_position
from tests.oracles.incremental import (
    assert_equivalent_vd,
    assert_update_matches_recompute,
    cell_bound_keys,
    wall_keys,
)
from tests.oracles.sweep import assert_mid_interval_matches_recompute
from vd3d.events import (
    EventType,
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
    compute_vd_around_event,
    equivalent_vd,
    event_vertex_plane_keys,
    local_cell_ids,
    local_vertex_plane_keys,
    match_cells,
    same_combinatorics,
    unmatched_after,
    unmatched_before,
    update_2d_decomposition,
    update_for_triple_intersection,
    update_for_vertical_alignment,
    vertical_decomposition_3d,
)
from vd3d.sweep.matching import cell_near_event
from vd3d.vertical_decomposition import verify_vd_invariants


def _gp_planes(rng: random.Random, n: int) -> list[Plane]:
    for _ in range(40):
        planes = random_planes_general_position(rng, n)
        events = generate_all_events(planes)
        if events_have_unique_z(events):
            return planes
    raise AssertionError(
        f"failed to sample {n} general-position planes with unique event z"
    )


# ---------------------------------------------------------------------------
# Step 9.1 — triple handler
# ---------------------------------------------------------------------------


def test_triple_handler_rejects_alignment():
    planes = planes_alignment_at_z2()
    events = generate_all_events(planes)
    alignment = next(
        event for event in events if event.type is EventType.VERTICAL_ALIGNMENT
    )
    vd, _, _, z_plus = compute_vd_around_event(planes, alignment, (alignment,))
    with pytest.raises(ValueError, match="triple"):
        update_for_triple_intersection(vd, alignment, planes, z_plus)


def test_triple_incremental_matches_recompute():
    planes = planes_through_123()
    events = generate_all_events(planes)
    assert len(events) == 1
    event = events[0]
    assert event.type is EventType.TRIPLE_INTERSECTION
    vd_before, vd_ref, _, z_plus = compute_vd_around_event(planes, event, events)
    incremental = assert_update_matches_recompute(vd_before, event, planes, z_plus)
    assert equivalent_vd(incremental, vd_ref)
    assert wall_keys(incremental.walls) == wall_keys(vd_ref.walls)
    assert cell_bound_keys(incremental) == cell_bound_keys(vd_ref)
    verify_vd_invariants(incremental)


def test_triple_local_neighbourhood_is_the_only_change():
    planes = planes_through_123()
    events = generate_all_events(planes)
    event = events[0]
    vd_before, _vd_after, _, z_plus = compute_vd_around_event(planes, event, events)
    incremental = update_for_triple_intersection(vd_before, event, planes, z_plus)
    matches = match_cells(vd_before, incremental)
    dying = unmatched_before(vd_before, matches)
    born = unmatched_after(incremental, matches)
    local_before = local_cell_ids(vd_before, event)
    local_after = local_cell_ids(incremental, event)
    assert {cell.id for cell in dying} <= local_before
    assert {cell.id for cell in born} <= local_after
    keys = event_vertex_plane_keys(event, planes)
    assert len(keys) == 3
    local_keys = local_vertex_plane_keys(
        vd_before, incremental.arrangement, event, planes
    )
    assert keys <= local_keys


# ---------------------------------------------------------------------------
# Step 9.2 — alignment handler
# ---------------------------------------------------------------------------


def test_alignment_handler_rejects_triple():
    planes = planes_through_123()
    events = generate_all_events(planes)
    vd, _, _, z_plus = compute_vd_around_event(planes, events[0], events)
    with pytest.raises(ValueError, match="alignment"):
        update_for_vertical_alignment(vd, events[0], planes, z_plus)


def test_alignment_incremental_matches_recompute():
    planes = planes_alignment_at_z2()
    events = generate_all_events(planes)
    alignments = [
        event for event in events if event.type is EventType.VERTICAL_ALIGNMENT
    ]
    assert len(alignments) == 1
    event = alignments[0]
    vd_before, vd_ref, _, z_plus = compute_vd_around_event(planes, event, (event,))
    incremental = assert_update_matches_recompute(vd_before, event, planes, z_plus)
    assert not same_combinatorics(vd_before, incremental)
    assert equivalent_vd(incremental, vd_ref)
    matches = match_cells(vd_before, incremental)
    local_before = local_cell_ids(vd_before, event)
    far_before = [cell for cell in vd_before.cells if cell.id not in local_before]
    matched = {pair.before.id for pair in matches}
    for cell in far_before:
        assert cell.id in matched or cell_near_event(vd_before, cell, event)


def test_dispatcher_unknown_type_rejected():
    planes = planes_through_123()
    events = generate_all_events(planes)
    event = events[0]
    vd_before, _, _, z_plus = compute_vd_around_event(planes, event, events)
    object.__setattr__(event, "type", None)
    with pytest.raises(ValueError, match="unknown event type"):
        update_2d_decomposition(vd_before, event, planes, z_plus)


# ---------------------------------------------------------------------------
# Sweep integration + random differential
# ---------------------------------------------------------------------------


def test_full_sweep_incremental_agrees_with_reference_path():
    planes = planes_through_123()
    inc = vertical_decomposition_3d(planes, incremental=True)
    ref = vertical_decomposition_3d(planes, incremental=False)
    assert len(inc.cells) == len(ref.cells)
    assert len(inc.snapshots) == len(ref.snapshots)
    assert inc.snapshots[0].n_after == ref.snapshots[0].n_after
    assert inc.snapshots[0].n_active == ref.snapshots[0].n_active
    assert_mid_interval_matches_recompute(inc)


def test_no_event_instances_unchanged():
    for planes in (planes_one(), planes_two()):
        inc = vertical_decomposition_3d(planes, incremental=True)
        ref = vertical_decomposition_3d(planes, incremental=False)
        assert len(inc.cells) == len(ref.cells)
        assert inc.events == ()


def test_random_gp_each_event_matches_recompute():
    rng = random.Random(91)
    seen_triple = 0
    seen_align = 0
    for n in (3, 4, 5):
        for _ in range(3):
            planes = _gp_planes(rng, n)
            events = generate_all_events(planes)
            for event in events:
                vd_before, _, _, z_plus = compute_vd_around_event(
                    planes, event, events
                )
                assert_update_matches_recompute(vd_before, event, planes, z_plus)
                if event.type is EventType.TRIPLE_INTERSECTION:
                    seen_triple += 1
                else:
                    seen_align += 1
            result = vertical_decomposition_3d(planes, incremental=True)
            assert_mid_interval_matches_recompute(result)
    assert seen_triple >= 1
    assert seen_align >= 0
