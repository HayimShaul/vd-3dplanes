"""Exact random planes and points for the review viewer."""

from __future__ import annotations

import random
import time
from collections.abc import Sequence

from vd3d.geometry.intersections import PARALLEL, intersect_planes, intersect_three_planes
from vd3d.geometry.linalg import det3
from vd3d.geometry.plane import Plane, normals_parallel
from vd3d.geometry.points import Point3D
from vd3d.geometry.scalar import Scalar, as_scalar
from vd3d.viz.convert import to_float


def choose_seed(seed: int | None) -> int:
    """Return ``seed``, or the current time in nanoseconds if none was given."""
    if seed is None:
        return time.time_ns()
    return int(seed)


def resolve_n(n: int | None, rng: random.Random, choices: tuple[int, ...]) -> int:
    """Line/plane count: ``n`` if given, otherwise a value from ``choices``."""
    if n is None:
        return rng.choice(choices)
    if n < 1:
        raise ValueError("n must be a positive integer")
    return int(n)


def format_plane(plane: Plane) -> str:
    return f"{plane.a}x+{plane.b}y+{plane.c}z+{plane.d}=0"


def format_point(point: Point3D) -> str:
    return f"({point.x}, {point.y}, {point.z})"


def view_lim(*points: Point3D, minimum: float = 3.0) -> float:
    lim = minimum
    for point in points:
        lim = max(lim, abs(to_float(point.x)), abs(to_float(point.y)), abs(to_float(point.z)))
    return lim + 1.0


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


def random_point(rng: random.Random, lo: int = -3, hi: int = 3) -> Point3D:
    return Point3D(rng.randint(lo, hi), rng.randint(lo, hi), rng.randint(lo, hi))


def point_on_plane(plane: Plane, rng: random.Random) -> Point3D:
    """A random exact point on ``plane``, by solving one coordinate."""
    u = as_scalar(rng.randint(-2, 2))
    v = as_scalar(rng.randint(-2, 2))
    if plane.c != 0:
        z = -(plane.a * u + plane.b * v + plane.d) / plane.c
        return Point3D(u, v, z)
    if plane.b != 0:
        y = -(plane.a * u + plane.c * v + plane.d) / plane.b
        return Point3D(u, y, v)
    x = -(plane.b * u + plane.c * v + plane.d) / plane.a
    return Point3D(x, u, v)


def offset_along_normal(point: Point3D, plane: Plane, sign: int) -> Point3D:
    dx, dy, dz = plane.normal
    return point.translated(sign * dx, sign * dy, sign * dz)


def random_intersecting_planes(rng: random.Random) -> tuple[Plane, Plane]:
    for _ in range(200):
        p = random_plane(rng, 1)
        q = random_plane(rng, 2)
        if not normals_parallel(p, q):
            return p, q
    raise RuntimeError("failed to sample intersecting planes")


def random_parallel_planes(rng: random.Random) -> tuple[Plane, Plane]:
    p = random_plane(rng, 1)
    offset = as_scalar(random_int(rng, -3, 3, nonzero=True))
    q = Plane(id=2, a=p.a, b=p.b, c=p.c, d=p.d + offset)
    return p, q


def random_plane_through(rng: random.Random, point: Point3D, plane_id: int) -> Plane:
    for _ in range(200):
        a = random_int(rng)
        b = random_int(rng)
        c = random_int(rng)
        if a == 0 and b == 0 and c == 0:
            continue
        d = -(as_scalar(a) * point.x + as_scalar(b) * point.y + as_scalar(c) * point.z)
        return Plane(id=plane_id, a=a, b=b, c=c, d=d)
    raise RuntimeError("failed to sample a plane through a point")


def random_triple_planes(rng: random.Random) -> tuple[Point3D, Plane, Plane, Plane]:
    for _ in range(200):
        point = random_point(rng)
        p = random_plane_through(rng, point, 1)
        q = random_plane_through(rng, point, 2)
        r = random_plane_through(rng, point, 3)
        if det3((p.normal, q.normal, r.normal)) == 0:
            continue
        return point, p, q, r
    raise RuntimeError("failed to sample three planes through one point")


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


def _has_vertical_intersection(candidate: Plane, planes: Sequence[Plane]) -> bool:
    for other in planes:
        line = intersect_planes(candidate, other)
        if line is PARALLEL:
            continue
        dx, dy, _dz = line.direction
        if dx == 0 and dy == 0:
            return True
    return False


def display_t_range(direction: tuple[Scalar, Scalar, Scalar], target: float = 3.0) -> tuple[float, float]:
    mag = max(abs(to_float(c)) for c in direction)
    if mag == 0:
        return -1.0, 1.0
    span = target / mag
    return -span, span
