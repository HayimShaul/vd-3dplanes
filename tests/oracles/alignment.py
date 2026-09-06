"""Brute-force alignment oracle. This is the definition of correctness."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events import Event, canonical_event_key, enumerate_alignment_pairs
from vd3d.events.samples import (
    planes_alignment_at_z2,
    planes_alignment_blocked,
    planes_never_align,
)
from vd3d.geometry.plane import Plane

__all__ = [
    "alignment_keys",
    "brute_force_alignment_events",
    "planes_alignment_at_z2",
    "planes_alignment_blocked",
    "planes_never_align",
]


def brute_force_alignment_events(planes: Sequence[Plane]) -> tuple[Event, ...]:
    """All pairs of intersection lines; drop triples; require visibility."""
    return enumerate_alignment_pairs(planes)


def alignment_keys(events: Sequence[Event]) -> set[tuple]:
    return {canonical_event_key(event) for event in events}
