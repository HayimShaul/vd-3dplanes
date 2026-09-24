"""Exact random plane sampling for the kernel and CLI.

Must not import ``vd3d.viz`` or floating-point plotting libraries.
"""

from __future__ import annotations

import random
from collections.abc import Sequence

from vd3d.geometry.intersections import PARALLEL, intersect_planes, intersect_three_planes
from vd3d.geometry.plane import Plane, normals_parallel
from vd3d.geometry.scalar import as_scalar


def random_int(rng: random.Random, lo: int = -4, hi: int = 4, *, nonzero: bool = False) -> int:
    if not nonzero:
        return rng.randint(lo, hi)
    choices = [n for n in range(lo, hi + 1) if n != 0]
    return rng.choice(choices)


def random_plane(
    rng: random.Random,
    plane_id: int,
    *,
    nonzero_xy: bool = False,
    nonzero_z: bool = False,
) -> Plane:
    for _ in range(200):
        a = random_int(rng)
        b = random_int(rng)
        c = random_int(rng)
        d = random_int(rng)
        if a == 0 and b == 0 and c == 0:
            continue
        if nonzero_xy and a == 0 and b == 0:
            continue
        if nonzero_z and c == 0:
            continue
        return Plane(id=plane_id, a=a, b=b, c=c, d=d)
    raise RuntimeError("failed to sample a plane")


def random_general_position_planes(rng: random.Random, n: int) -> list[Plane]:
    """``n`` planes: no parallels, every triple meets, no vertical intersection line."""
    if n < 1:
        raise ValueError("n must be a positive integer")
    planes: list[Plane] = []
    for plane_id in range(1, n + 1):
        for _ in range(400):
            candidate = random_plane(rng, plane_id, nonzero_xy=True, nonzero_z=True)
            if any(normals_parallel(candidate, other) for other in planes):
                continue
            if any(
                intersect_three_planes(p, q, candidate) is None
                for i, p in enumerate(planes)
                for q in planes[i + 1 :]
            ):
                continue
            if _has_vertical_intersection(candidate, planes):
                continue
            planes.append(candidate)
            break
        else:
            raise RuntimeError(f"failed to sample general-position plane {plane_id}")
    return planes


def _has_vertical_intersection(candidate: Plane, planes: Sequence[Plane]) -> bool:
    for other in planes:
        line = intersect_planes(candidate, other)
        if line is PARALLEL:
            continue
        dx, dy, _dz = line.direction
        if dx == 0 and dy == 0:
            return True
    return False


def random_parallel_family(rng: random.Random, n: int) -> list[Plane]:
    """``n`` pairwise-parallel planes (no intersection lines, no triples)."""
    if n < 1:
        raise ValueError("n must be a positive integer")
    first = random_plane(rng, 1, nonzero_z=True)
    planes = [first]
    used_offsets = {as_scalar(0)}
    for plane_id in range(2, n + 1):
        for _ in range(40):
            offset = as_scalar(random_int(rng, -5, 5, nonzero=True))
            if offset in used_offsets:
                continue
            used_offsets.add(offset)
            planes.append(
                Plane(id=plane_id, a=first.a, b=first.b, c=first.c, d=first.d + offset)
            )
            break
        else:
            raise RuntimeError("failed to sample a parallel family")
    return planes
