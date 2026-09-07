"""Sweep data types (design §13, §20)."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.events.types import Event
from vd3d.geometry.scalar import Scalar
from vd3d.vertical_decomposition.types import VDCell2D


class SimultaneousEvents(ValueError):
    """Two events share a ``z``. Raised only when general position is required."""


@dataclass(frozen=True, slots=True)
class EventGroup:
    """Events that share a ``z`` (design §1.4). General position: length 1."""

    z: Scalar
    events: tuple[Event, ...]

    @property
    def representative(self) -> Event:
        return self.events[0]


@dataclass(frozen=True, slots=True)
class CellSignature:
    """Combinatorial identity of a 2D trapezoid, stable between events."""

    lower_plane_id: int | None
    upper_plane_id: int | None
    left_unbounded: bool
    right_unbounded: bool
    left_planes: tuple[int, ...]
    right_planes: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class CellMatch:
    """A 2D cell before an event that continues as a 2D cell after it."""

    before: VDCell2D
    after: VDCell2D


@dataclass(frozen=True, slots=True)
class SweepSnapshot:
    """Debugging record of one event (design §20)."""

    event: Event
    n_before: int
    n_after: int
    matches: tuple[tuple[int, int], ...]
    started: tuple[int, ...]
    ended: tuple[int, ...]
    n_active: int
    group_ids: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, object]:
        return {
            "event_id": self.event.stable_id,
            "event_z": str(self.event.z),
            "event_type": self.event.type.name,
            "n_before": self.n_before,
            "n_after": self.n_after,
            "matches": [{"before": a, "after": b} for a, b in self.matches],
            "started": list(self.started),
            "ended": list(self.ended),
            "continued": len(self.matches),
            "n_active": self.n_active,
            "group_ids": list(self.group_ids),
        }


@dataclass(frozen=True, slots=True)
class ZInterval:
    """Open z-interval with the active 3D cell for each 2D signature."""

    lower_z: Scalar | None
    upper_z: Scalar | None
    binding: tuple[tuple[CellSignature, int], ...]

    def lookup(self) -> dict[CellSignature, int]:
        return dict(self.binding)
