"""z±ε around an event and the two recomputed 2D VDs (design §11, §19)."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events.slice import compute_vd_at_z
from vd3d.events.types import Event
from vd3d.geometry.plane import Plane
from vd3d.geometry.scalar import Scalar, as_scalar
from vd3d.sweep.types import SimultaneousEvents
from vd3d.vertical_decomposition.types import VerticalDecomposition

DEFAULT_EVENT_EPS = as_scalar("1/100")


def z_before_after(
    event: Event,
    events: Sequence[Event],
    *,
    default_eps: int | Scalar | str = DEFAULT_EVENT_EPS,
) -> tuple[Scalar, Scalar]:
    """Heights strictly below / above ``event.z``, inside the neighbouring gaps."""
    z = event.z
    eps = as_scalar(default_eps)
    if eps <= 0:
        raise ValueError("default_eps must be positive")
    others = [other.z for other in events if other.stable_id != event.stable_id]
    if others:
        gap = min(abs(other - z) for other in others)
        if gap == 0:
            raise SimultaneousEvents(
                f"event {event.stable_id} shares z={z} with another event"
            )
        half = gap / 2
        if half < eps:
            eps = half
    return z - eps, z + eps


def compute_vd_around_event(
    planes: Sequence[Plane],
    event: Event,
    events: Sequence[Event] | None = None,
    *,
    default_eps: int | Scalar | str = DEFAULT_EVENT_EPS,
) -> tuple[VerticalDecomposition, VerticalDecomposition, Scalar, Scalar]:
    """Recompute the 2D VD at ``z-ε`` and ``z+ε`` (reference, not incremental)."""
    if events is None:
        events = (event,)
    z_minus, z_plus = z_before_after(event, events, default_eps=default_eps)
    return (
        compute_vd_at_z(planes, z_minus),
        compute_vd_at_z(planes, z_plus),
        z_minus,
        z_plus,
    )
