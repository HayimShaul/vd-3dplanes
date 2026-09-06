"""Shared checks for event-generation inputs."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.geometry.plane import Plane


def require_unique_plane_ids(planes: Sequence[Plane]) -> None:
    ids = [plane.id for plane in planes]
    if len(set(ids)) != len(ids):
        raise ValueError("plane ids must be unique")
