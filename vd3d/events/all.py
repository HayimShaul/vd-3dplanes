"""Complete event list (design §9) and a z below every event (design §10)."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events.alignment import deduplicate_events, generate_alignment_events
from vd3d.events.ids import require_unique_plane_ids
from vd3d.events.triples import generate_triple_events
from vd3d.events.types import Event
from vd3d.geometry.plane import Plane
from vd3d.geometry.scalar import Scalar, as_scalar


def generate_all_events(planes: Sequence[Plane]) -> tuple[Event, ...]:
    """Triples plus alignments, deduplicated, sorted by ``(z, type, id)``."""
    require_unique_plane_ids(planes)
    combined = (
        *generate_triple_events(planes),
        *generate_alignment_events(planes),
    )
    return tuple(sorted(deduplicate_events(combined), key=lambda event: event.sort_key))


def choose_z_below_all_events(
    events: Sequence[Event],
    *,
    margin: int | Scalar | str = 1,
) -> Scalar:
    """A height strictly below every event ``z``.

    ``margin`` is an exact positive offset subtracted from the minimum
    event ``z``. With no events the slice combinatorics never change, so
    the conventional height is ``0``.
    """
    offset = as_scalar(margin)
    if offset <= 0:
        raise ValueError("margin must be positive")
    if not events:
        return as_scalar(0)
    return min(event.z for event in events) - offset
