"""Exact intersections among planes and horizontal slices."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.geometry.linalg import cross, solve2, solve3
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane, normals_parallel
from vd3d.geometry.points import Point3D
from vd3d.geometry.scalar import Scalar, as_scalar


@dataclass(frozen=True, slots=True)
class ParallelPlanes:
    """Sentinel: two planes have parallel (or anti-parallel) normals."""


PARALLEL = ParallelPlanes()


def intersect_planes(p: Plane, q: Plane) -> Line3D | ParallelPlanes:
    """Intersect two planes.

    If the normals are parallel, return ``PARALLEL`` (this includes the
    coincident-plane case; we do not distinguish it yet).
    """
    if normals_parallel(p, q):
        return PARALLEL
    direction = cross(p.normal, q.normal)
    point = _point_on_both_planes(p, q, direction)
    return Line3D(point=point, direction=direction, plane_a=p.id, plane_b=q.id)


def intersect_three_planes(p: Plane, q: Plane, r: Plane) -> Point3D | None:
    """Unique intersection of three planes, or ``None`` if the 3x3 is singular."""
    matrix = (
        (p.a, p.b, p.c),
        (q.a, q.b, q.c),
        (r.a, r.b, r.c),
    )
    rhs = (-p.d, -q.d, -r.d)
    solution = solve3(matrix, rhs)
    if solution is None:
        return None
    return Point3D(*solution)


def slice_plane_at_z(plane: Plane, z: int | Scalar | str) -> Line2D | None:
    """Intersect ``plane`` with the horizontal plane at height ``z``.

    Substituting ``z`` into ``a x + b y + c z + d = 0`` yields
    ``a x + b y + (c z + d) = 0``. Returns ``None`` when ``a = b = 0``
    (horizontal plane: the slice is empty or the whole plane).
    """
    z = as_scalar(z)
    a = plane.a
    b = plane.b
    c = plane.c * z + plane.d
    if a == 0 and b == 0:
        return None
    return Line2D(a=a, b=b, c=c, source_plane_id=plane.id)


def build_slice_lines(planes: list[Plane], z: int | Scalar | str) -> list[Line2D]:
    """Slice every plane at ``z`` and drop degenerate (horizontal) slices."""
    lines: list[Line2D] = []
    for plane in planes:
        line = slice_plane_at_z(plane, z)
        if line is not None:
            lines.append(line)
    return lines


def _point_on_both_planes(
    p: Plane,
    q: Plane,
    direction: tuple[Scalar, Scalar, Scalar],
) -> Point3D:
    """One exact point on both planes.

    Set the coordinate matching a nonzero direction component to 0 and
    solve the remaining 2x2. That 2x2 determinant equals plus or minus
    the chosen direction component, so it is nonzero.
    """
    dx, dy, dz = direction
    if dz != 0:
        sol = solve2(((p.a, p.b), (q.a, q.b)), (-p.d, -q.d))
        if sol is None:
            raise RuntimeError("expected nonsingular 2x2 when direction.z != 0")
        return Point3D(sol[0], sol[1], 0)
    if dy != 0:
        sol = solve2(((p.a, p.c), (q.a, q.c)), (-p.d, -q.d))
        if sol is None:
            raise RuntimeError("expected nonsingular 2x2 when direction.y != 0")
        return Point3D(sol[0], 0, sol[1])
    sol = solve2(((p.b, p.c), (q.b, q.c)), (-p.d, -q.d))
    if sol is None:
        raise RuntimeError("expected nonsingular 2x2 when direction.x != 0")
    return Point3D(0, sol[0], sol[1])
