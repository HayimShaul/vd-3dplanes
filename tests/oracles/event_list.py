"""Brute-force checkers for the combined event list and a slice VD."""

from __future__ import annotations

from collections.abc import Sequence

from tests.oracles.alignment import brute_force_alignment_events
from vd3d.events import (
    Event,
    canonical_event_key,
    deduplicate_events,
    generate_triple_events,
)
from vd3d.events.samples import planes_parallel_family, planes_through_123
from vd3d.geometry.plane import Plane

__all__ = [
    "brute_force_all_events",
    "event_keys",
    "planes_parallel_family",
    "planes_through_123",
]


def brute_force_all_events(planes: Sequence[Plane]) -> tuple[Event, ...]:
    """Union of the triple enumerator and the alignment pair-oracle."""
    combined = (
        *generate_triple_events(planes),
        *brute_force_alignment_events(planes),
    )
    return tuple(sorted(deduplicate_events(combined), key=lambda event: event.sort_key))


def event_keys(events: Sequence[Event]) -> set[tuple]:
    return {canonical_event_key(event) for event in events}
