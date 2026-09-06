"""Pairwise intersection lines of a set of planes (design §6)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace

from vd3d.events.ids import require_unique_plane_ids
from vd3d.geometry.intersections import PARALLEL, intersect_planes
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane


def compute_intersection_lines(planes: Sequence[Plane]) -> tuple[Line3D, ...]:
    """Every pairwise intersection, skipping parallel (and coincident) pairs.

    Line ids are ``0 .. m-1`` in the order of plane pairs ``(i, j)`` with
    ``i < j``. ``plane_a`` / ``plane_b`` are the source plane ids.
    """
    require_unique_plane_ids(planes)
    lines: list[Line3D] = []
    next_id = 0
    for i, p in enumerate(planes):
        for q in planes[i + 1 :]:
            result = intersect_planes(p, q)
            if result is PARALLEL:
                continue
            lines.append(replace(result, id=next_id))
            next_id += 1
    return tuple(lines)
