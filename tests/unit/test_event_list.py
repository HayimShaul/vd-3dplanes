"""Phase 7 — combined event list and VD at a given z. Design §9–10."""

from __future__ import annotations

import random

import pytest

from tests.oracles.alignment import planes_alignment_at_z2
from tests.oracles.event_list import brute_force_all_events, event_keys
from tests.oracles.events import random_planes_general_position
from vd3d.events import (
    EventType,
    build_2d_arrangement_at_z,
    choose_z_below_all_events,
    compute_vd_at_z,
    events_have_unique_z,
    generate_alignment_events,
    generate_all_events,
    generate_triple_events,
    initial_slice,
    verify_all_events,
    z_strictly_below_all_events,
)
from vd3d.events.samples import planes_parallel_family, planes_through_123
from vd3d.geometry import Plane
from vd3d.geometry.scalar import as_scalar
from vd3d.vertical_decomposition import verify_vd_invariants


# ---------------------------------------------------------------------------
# Step 7.1 — generate_all_events and choose_z_below_all_events
# ---------------------------------------------------------------------------


def test_all_events_known_triple():
    planes = planes_through_123()
    events = generate_all_events(planes)
    verify_all_events(events)
    assert events_have_unique_z(events)
    assert len(events) == 1
    assert events[0].type is EventType.TRIPLE_INTERSECTION
    assert events[0].z == 3
    assert events[0].plane_ids == (1, 2, 3)
    assert event_keys(events) == event_keys(brute_force_all_events(planes))


def test_all_events_parallel_family_empty():
    assert generate_all_events([]) == ()
    events = generate_all_events(planes_parallel_family())
    assert events == ()
    verify_all_events(events)
    assert choose_z_below_all_events(events) == 0


def test_all_events_is_triples_union_alignments():
    planes = planes_alignment_at_z2()
    events = generate_all_events(planes)
    verify_all_events(events)
    assert event_keys(events) == event_keys(generate_triple_events(planes)) | event_keys(
        generate_alignment_events(planes)
    )
    assert event_keys(events) == event_keys(brute_force_all_events(planes))
    types = {event.type for event in events}
    assert EventType.TRIPLE_INTERSECTION in types
    assert EventType.VERTICAL_ALIGNMENT in types
    keys = [event.sort_key for event in events]
    assert keys == sorted(keys)
    assert len(keys) == len(set(keys))


def test_all_events_sorted_unique_keys_random():
    rng = random.Random(71)
    for n in (3, 4, 5):
        for _ in range(4):
            planes = random_planes_general_position(rng, n)
            events = generate_all_events(planes)
            verify_all_events(events)
            assert event_keys(events) == event_keys(brute_force_all_events(planes))
            keys = [event.sort_key for event in events]
            assert keys == sorted(keys)
            assert len(keys) == len(set(keys))


def test_general_position_fixtures_have_unique_z():
    planes = planes_through_123()
    events = generate_all_events(planes)
    assert events_have_unique_z(events)

    rng = random.Random(72)
    seen = 0
    for n in (3, 4, 5):
        for _ in range(8):
            planes = random_planes_general_position(rng, n)
            events = generate_all_events(planes)
            if not events_have_unique_z(events):
                continue
            zs = [event.z for event in events]
            assert zs == sorted(zs)
            assert len(zs) == len(set(zs))
            seen += 1
    assert seen >= 6


def test_choose_z_below_all_events():
    planes = planes_through_123()
    events = generate_all_events(planes)
    z0 = choose_z_below_all_events(events)
    assert z0 == 2
    assert z_strictly_below_all_events(z0, events)
    assert choose_z_below_all_events(events, margin="3/2") == events[0].z - as_scalar("3/2")
    with pytest.raises(ValueError, match="positive"):
        choose_z_below_all_events(events, margin=0)


def test_generate_all_events_reject_duplicate_ids():
    planes = [
        Plane(id=1, a=1, b=0, c=1, d=0),
        Plane(id=1, a=0, b=1, c=1, d=0),
        Plane(id=2, a=1, b=1, c=1, d=-1),
    ]
    with pytest.raises(ValueError, match="unique"):
        generate_all_events(planes)


# ---------------------------------------------------------------------------
# Step 7.2 — slice + VD at a given z
# ---------------------------------------------------------------------------


def test_vd_at_triple_concurrent_vs_triangle():
    """Just below / at / just above Test 15: 3 vertices vs 1 concurrent."""
    planes = planes_through_123()
    events = generate_all_events(planes)
    z_event = events[0].z
    assert z_event == 3

    vd_event = compute_vd_at_z(planes, z_event)
    vd_below = compute_vd_at_z(planes, z_event - 1)
    vd_above = compute_vd_at_z(planes, z_event + 1)
    for vd in (vd_event, vd_below, vd_above):
        verify_vd_invariants(vd)
        assert len(vd.arrangement.lines) == 3

    assert len(vd_event.arrangement.vertices) == 1
    assert len(vd_below.arrangement.vertices) == 3
    assert len(vd_above.arrangement.vertices) == 3
    assert len(vd_below.cells) == len(vd_above.cells)
    assert len(vd_event.cells) != len(vd_below.cells)


def test_slice_lines_remember_source_planes():
    planes = planes_through_123()
    arrangement = build_2d_arrangement_at_z(planes, 2)
    assert {line.source_plane_id for line in arrangement.lines} == {1, 2, 3}
    assert [line.id for line in arrangement.lines] == [0, 1, 2]


def test_initial_slice_is_vd_below_events():
    planes = planes_through_123()
    events = generate_all_events(planes)
    z0 = choose_z_below_all_events(events)
    vd = initial_slice(planes, events)
    verify_vd_invariants(vd)
    below = compute_vd_at_z(planes, z0)
    assert len(vd.cells) == len(below.cells)
    assert len(vd.arrangement.vertices) == len(below.arrangement.vertices) == 3


def test_compute_vd_at_z_empty_and_singleton():
    empty = compute_vd_at_z([], 0)
    verify_vd_invariants(empty)
    assert len(empty.cells) == 1
    assert empty.arrangement.vertices == ()

    one = compute_vd_at_z([Plane(id=1, a=1, b=0, c=1, d=0)], 0)
    verify_vd_invariants(one)
    assert len(one.arrangement.lines) == 1
    assert one.arrangement.vertices == ()


def test_compute_vd_at_z_reject_duplicate_ids():
    planes = [
        Plane(id=7, a=1, b=0, c=1, d=0),
        Plane(id=7, a=0, b=1, c=1, d=0),
    ]
    with pytest.raises(ValueError, match="unique"):
        compute_vd_at_z(planes, 0)
