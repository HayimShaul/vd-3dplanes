"""Sweep-event types. Alignment events are generated in Phase 6."""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

from vd3d.geometry.points import Point3D
from vd3d.geometry.scalar import Scalar, as_scalar


class EventType(IntEnum):
    """Kinds of ``z``-sweep events.

    Integer values are the middle component of ``sort_key``.
    """

    TRIPLE_INTERSECTION = 0
    VERTICAL_ALIGNMENT = 1


@dataclass(frozen=True, slots=True)
class AlignmentGeometry:
    """The two slice-vertices that share an ``x`` at the alignment ``z``."""

    point_a: Point3D
    point_b: Point3D


def default_event_stable_id(event_type: EventType, plane_ids: tuple[int, ...]) -> str:
    """Canonical id from the event kind and sorted plane ids."""
    ids = ":".join(str(i) for i in sorted(plane_ids))
    if event_type is EventType.TRIPLE_INTERSECTION:
        return f"triple:{ids}"
    if event_type is EventType.VERTICAL_ALIGNMENT:
        return f"align:{ids}"
    raise ValueError(f"unknown event type: {event_type!r}")


@dataclass(frozen=True, slots=True)
class Event:
    """A sweep event at a single ``z``.

    ``type`` is the design name. ``sort_key`` is ``(z, type, stable_id)``.
    For a triple, ``geometric_data`` is the unique three-plane point.
    For an alignment, it is ``AlignmentGeometry`` (the two slice vertices).
    """

    z: Scalar
    type: EventType
    geometric_data: Point3D | AlignmentGeometry
    plane_ids: tuple[int, ...]
    stable_id: str = ""
    line_ids: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "z", as_scalar(self.z))
        object.__setattr__(self, "plane_ids", tuple(self.plane_ids))
        object.__setattr__(self, "line_ids", tuple(self.line_ids))
        if self.stable_id == "":
            object.__setattr__(
                self,
                "stable_id",
                default_event_stable_id(self.type, self.plane_ids),
            )

    @property
    def sort_key(self) -> tuple[Scalar, int, str]:
        return (self.z, int(self.type), self.stable_id)
