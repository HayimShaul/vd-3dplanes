"""Phase 10 — simultaneous groups, zone degeneracies, non-general position."""

from __future__ import annotations

import random

import pytest

from tests.oracles.incremental import assert_equivalent_vd, cell_bound_keys, wall_keys
from tests.oracles.simple_arrangement import random_simple_lines
from tests.oracles.sweep import (
    assert_grid_partition_3d,
    assert_mid_interval_matches_recompute,
    assert_point_location_partition,
    sample_box_points,
)
from tests.oracles.zone import brute_force_crossings, brute_force_face_ids
from vd3d.arrangement2d import build_line_arrangement
from vd3d.cells3d import verify_cells3d
from vd3d.events import EventType, generate_all_events, generate_triple_events
from vd3d.events.samples import (
    planes_alignment_at_z2,
    planes_four_through_123,
    planes_through_123,
    planes_vertical_and_slanted,
)
from vd3d.geometry import Line2D, Point2D, slice_plane_at_z
from vd3d.sweep import (
    SimultaneousEvents,
    compute_vd_around_event,
    group_events_by_z,
    update_2d_for_event_group,
    vertical_decomposition_3d,
    z_before_after,
)
from vd3d.vertical_decomposition import verify_vd_invariants
from vd3d.zone import (
    coincident_line_indices,
    collapse_crossings,
    compute_crossings,
    compute_supporting_line_zone,
    compute_zone,
    verify_zone_invariants,
)


def _assert_crossings_match(left, right) -> None:
    assert len(left) == len(right)
    for a, b in zip(left, right):
        assert a.point == b.point
        assert a.t == b.t
        assert a.edge_id == b.edge_id
        assert a.vertex_id == b.vertex_id


# ---------------------------------------------------------------------------
# Step 10.1 — zone through a vertex
# ---------------------------------------------------------------------------


def test_through_vertex_random_matches_oracle():
    rng = random.Random(101)
    for n in range(3, 7):
        arr = build_line_arrangement(random_simple_lines(rng, n))
        vertex = rng.choice(arr.vertices)
        point = vertex.point
        query = Line2D(a=1, b=1, c=-(point.x + point.y))
        if coincident_line_indices(arr, query):
            query = Line2D(a=1, b=2, c=-(point.x + 2 * point.y))
        assert query.contains(point)
        assert not coincident_line_indices(arr, query)
        zone = compute_zone(arr, query)
        verify_zone_invariants(arr, zone)
        assert any(v.id == vertex.id for v in zone.vertices)
        assert [f.id for f in zone.faces] == list(brute_force_face_ids(arr, query))
        raw = compute_crossings(arr, query)
        assert len(collapse_crossings(raw)) == len(zone.crossings)
        _assert_crossings_match(raw, brute_force_crossings(arr, query))


# ---------------------------------------------------------------------------
# Step 10.2 — overlapping query
# ---------------------------------------------------------------------------


def test_overlap_each_triangle_side():
    arr = build_line_arrangement(
        [
            Line2D(a=1, b=0, c=0, id=0),
            Line2D(a=0, b=1, c=0, id=1),
            Line2D(a=1, b=1, c=-1, id=2),
        ]
    )
    for line_index in (0, 1, 2):
        query = arr.lines[line_index]
        zone = compute_zone(arr, query)
        verify_zone_invariants(arr, zone)
        supporting = compute_supporting_line_zone(arr, line_index)
        assert {f.id for f in zone.faces} == {f.id for f in supporting.faces}
        assert {v.id for v in zone.vertices} == {v.id for v in supporting.vertices_on_line}
        assert [f.id for f in zone.faces] == list(brute_force_face_ids(arr, query))


def test_overlap_opposite_orientation():
    arr = build_line_arrangement(
        [Line2D(a=1, b=0, c=0, id=0), Line2D(a=0, b=1, c=0, id=1)]
    )
    query = Line2D(a=-1, b=0, c=0)
    assert coincident_line_indices(arr, query) == (0,)
    zone = compute_zone(arr, query)
    verify_zone_invariants(arr, zone)
    supporting = compute_supporting_line_zone(arr, 0)
    assert {f.id for f in zone.faces} == {f.id for f in supporting.faces}


# ---------------------------------------------------------------------------
# Step 10.3 — simultaneous event groups
# ---------------------------------------------------------------------------


def test_group_events_by_z_clusters_same_height():
    planes = planes_alignment_at_z2()
    events = generate_all_events(planes)
    groups = group_events_by_z(events)
    assert any(len(group.events) > 1 for group in groups)
    for group in groups:
        assert all(event.z == group.z for event in group.events)
    assert sum(len(group.events) for group in groups) == len(events)


def test_z_before_after_uses_gap_to_other_heights():
    planes = planes_alignment_at_z2()
    events = generate_all_events(planes)
    z_minus, z_plus = z_before_after(events[0], events)
    assert z_minus < events[0].z < z_plus
    for event in events:
        if event.z == events[0].z:
            continue
        assert event.z < z_minus or event.z > z_plus


def test_general_position_flag_still_rejects_groups():
    with pytest.raises(SimultaneousEvents):
        vertical_decomposition_3d(
            planes_alignment_at_z2(), require_general_position=True
        )


def test_alignment_fixture_full_sweep():
    planes = planes_alignment_at_z2()
    result = vertical_decomposition_3d(planes)
    verify_cells3d(result.cells)
    groups = group_events_by_z(result.events)
    assert len(result.snapshots) == len(groups)
    assert any(len(snap.group_ids) > 1 for snap in result.snapshots)
    assert_mid_interval_matches_recompute(result)
    assert_grid_partition_3d(result, half=2, step=1)
    rng = random.Random(110)
    assert_point_location_partition(result, sample_box_points(rng, n=24, half=3))


def test_group_incremental_matches_recompute():
    planes = planes_alignment_at_z2()
    events = generate_all_events(planes)
    groups = group_events_by_z(events)
    multi = next(group for group in groups if len(group.events) > 1)
    vd_before, vd_ref, _, z_plus = compute_vd_around_event(
        planes, multi.representative, events
    )
    incremental = update_2d_for_event_group(vd_before, multi.events, planes, z_plus)
    assert_equivalent_vd(incremental, vd_ref, label=f"group z={multi.z}")
    assert wall_keys(incremental.walls) == wall_keys(vd_ref.walls)
    assert cell_bound_keys(incremental) == cell_bound_keys(vd_ref)
    verify_vd_invariants(incremental)


# ---------------------------------------------------------------------------
# Step 10.4 — four planes at a point, vertical input planes
# ---------------------------------------------------------------------------


def test_four_planes_four_triples_same_z():
    planes = planes_four_through_123()
    triples = generate_triple_events(planes)
    at_3 = [event for event in triples if event.z == 3]
    assert len(at_3) == 4
    assert all(event.type is EventType.TRIPLE_INTERSECTION for event in at_3)
    groups = group_events_by_z(generate_all_events(planes))
    at_event = next(group for group in groups if group.z == 3)
    assert len(at_event.events) >= 4


def test_four_planes_sweep_oracles():
    planes = planes_four_through_123()
    result = vertical_decomposition_3d(planes)
    verify_cells3d(result.cells)
    assert_mid_interval_matches_recompute(result)
    assert_grid_partition_3d(result, half=2, step=1)
    rng = random.Random(111)
    assert_point_location_partition(result, sample_box_points(rng, n=24, half=3))


def test_vertical_plane_slice_is_z_invariant():
    planes = planes_vertical_and_slanted()
    vertical = planes[0]
    assert vertical.c == 0
    line0 = slice_plane_at_z(vertical, 0)
    line5 = slice_plane_at_z(vertical, 5)
    assert line0 is not None and line5 is not None
    assert line0.a == line5.a and line0.b == line5.b and line0.c == line5.c
    sample = Point2D(0, 1)
    assert line0.contains(sample)


def test_vertical_planes_sweep_oracles():
    planes = planes_vertical_and_slanted()
    triples = generate_triple_events(planes)
    assert len(triples) == 1
    assert triples[0].z == 0
    result = vertical_decomposition_3d(planes)
    verify_cells3d(result.cells)
    assert_mid_interval_matches_recompute(result)
    assert_grid_partition_3d(result, half=2, step=1)
    rng = random.Random(112)
    assert_point_location_partition(result, sample_box_points(rng, n=20, half=3))


def test_three_plane_gp_still_one_snapshot():
    result = vertical_decomposition_3d(planes_through_123())
    assert len(result.snapshots) == 1
    assert result.snapshots[0].group_ids == (result.events[0].stable_id,)
