"""Exact random planes and points for the review viewer."""

from __future__ import annotations

import random
import time

from vd3d.geometry.linalg import det3
from vd3d.geometry.plane import Plane, normals_parallel
from vd3d.geometry.points import Point3D
from vd3d.geometry.sampling import (
    random_general_position_planes,
    random_int,
    random_parallel_family,
    random_plane,
)
from vd3d.geometry.scalar import Scalar, as_scalar
from vd3d.viz.convert import to_float

__all__ = [
    "choose_seed",
    "display_t_range",
    "format_plane",
    "format_point",
    "offset_along_normal",
    "point_on_plane",
    "random_general_position_planes",
    "random_int",
    "random_intersecting_planes",
    "random_parallel_family",
    "random_parallel_planes",
    "random_plane",
    "random_plane_through",
    "random_point",
    "random_triple_planes",
    "resolve_n",
    "view_lim",
]


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


def display_t_range(direction: tuple[Scalar, Scalar, Scalar], target: float = 3.0) -> tuple[float, float]:
    mag = max(abs(to_float(c)) for c in direction)
    if mag == 0:
        return -1.0, 1.0
    span = target / mag
    return -span, span
