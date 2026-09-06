"""Brute-force checkers for pairwise lines and triple events."""

from __future__ import annotations

import random
from collections.abc import Sequence

from vd3d.geometry import (
    PARALLEL,
    Plane,
    intersect_planes,
    intersect_three_planes,
    normals_parallel,
)
from vd3d.geometry.points import Point3D


def brute_force_intersecting_pairs(planes: Sequence[Plane]) -> tuple[tuple[int, int], ...]:
    """Sorted ``(plane_a.id, plane_b.id)`` pairs that are not parallel."""
    pairs: list[tuple[int, int]] = []
    for i, p in enumerate(planes):
        for q in planes[i + 1 :]:
            if intersect_planes(p, q) is PARALLEL:
                continue
            a, b = sorted((p.id, q.id))
            pairs.append((a, b))
    return tuple(pairs)


def brute_force_triple_points(
    planes: Sequence[Plane],
) -> tuple[tuple[Point3D, tuple[int, int, int]], ...]:
    """Every unique three-plane point with sorted plane ids."""
    found: list[tuple[Point3D, tuple[int, int, int]]] = []
    for i, p in enumerate(planes):
        for j, q in enumerate(planes[i + 1 :], start=i + 1):
            for r in planes[j + 1 :]:
                point = intersect_three_planes(p, q, r)
                if point is None:
                    continue
                found.append((point, tuple(sorted((p.id, q.id, r.id)))))
    found.sort(key=lambda item: (item[0].z, item[1]))
    return tuple(found)


def random_planes_general_position(
    rng: random.Random,
    n: int,
    *,
    require_all_triples: bool = True,
) -> list[Plane]:
    """``n`` planes with unique ids, no parallel pair, optionally every triple meets."""
    planes: list[Plane] = []
    for plane_id in range(n):
        for _ in range(400):
            a = rng.randint(-4, 4)
            b = rng.randint(-4, 4)
            c = rng.randint(-4, 4)
            d = rng.randint(-4, 4)
            if a == 0 and b == 0 and c == 0:
                continue
            candidate = Plane(id=plane_id, a=a, b=b, c=c, d=d)
            if any(normals_parallel(candidate, other) for other in planes):
                continue
            if require_all_triples and any(
                intersect_three_planes(p, q, candidate) is None
                for i, p in enumerate(planes)
                for q in planes[i + 1 :]
            ):
                continue
            planes.append(candidate)
            break
        else:
            raise RuntimeError(f"failed to sample plane {plane_id} of {n}")
    return planes
