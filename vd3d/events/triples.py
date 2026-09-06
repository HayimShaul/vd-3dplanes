"""Triple-intersection events (design §7)."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events.ids import require_unique_plane_ids
from vd3d.events.types import Event, EventType, default_event_stable_id
from vd3d.geometry.intersections import intersect_three_planes
from vd3d.geometry.plane import Plane


def generate_triple_events(planes: Sequence[Plane]) -> tuple[Event, ...]:
    """One event per triple of planes that meet at a unique point.

    Events are sorted by ``(z, type, stable_id)``. Parallel pencils and
    any singular 3×3 produce no event.
    """
    require_unique_plane_ids(planes)
    events: list[Event] = []
    for i, p in enumerate(planes):
        for j, q in enumerate(planes[i + 1 :], start=i + 1):
            for r in planes[j + 1 :]:
                point = intersect_three_planes(p, q, r)
                if point is None:
                    continue
                plane_ids = tuple(sorted((p.id, q.id, r.id)))
                events.append(
                    Event(
                        z=point.z,
                        type=EventType.TRIPLE_INTERSECTION,
                        geometric_data=point,
                        plane_ids=plane_ids,
                        stable_id=default_event_stable_id(
                            EventType.TRIPLE_INTERSECTION, plane_ids
                        ),
                    )
                )
    return tuple(sorted(events, key=lambda event: event.sort_key))
