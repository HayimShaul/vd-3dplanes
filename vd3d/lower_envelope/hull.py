"""Exact lower 3D convex hull of dual points (Fraction arithmetic).

A triangle is a lower-hull facet when it can be oriented with outward normal
``n_z < 0`` and every other point lies on the inward side of its plane.
"""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.geometry.linalg import cross
from vd3d.geometry.scalar import Scalar

Vec3 = tuple[Scalar, Scalar, Scalar]


def _sub(u: Vec3, v: Vec3) -> Vec3:
    return (u[0] - v[0], u[1] - v[1], u[2] - v[2])


def _dot(u: Vec3, v: Vec3) -> Scalar:
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def is_lower_hull_facet(
    points: Sequence[Vec3], i: int, j: int, k: int
) -> tuple[int, int, int] | None:
    """If ``(i,j,k)`` is a lower-hull facet, return it oriented with ``n_z < 0``."""
    pi, pj, pk = points[i], points[j], points[k]
    n = cross(_sub(pj, pi), _sub(pk, pi))
    if n == (0, 0, 0) or n[2] == 0:
        return None
    order = (i, j, k)
    if n[2] > 0:
        n = (-n[0], -n[1], -n[2])
        order = (i, k, j)
    # All other points on the inward side: (q - pi) · n <= 0
    for t, q in enumerate(points):
        if t in (i, j, k):
            continue
        if _dot(_sub(q, points[order[0]]), n) > 0:
            return None
    return order


def lower_hull_facets(points: Sequence[Vec3]) -> tuple[tuple[int, int, int], ...]:
    """Enumerate lower-hull triangular facets of ``points``.

    Expected ``O(n^4)``; fine for the small instances this package targets.
    """
    n = len(points)
    if n < 3:
        return ()
    if n == 3:
        facet = is_lower_hull_facet(points, 0, 1, 2)
        return (facet,) if facet is not None else ()
    facets: list[tuple[int, int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                facet = is_lower_hull_facet(points, i, j, k)
                if facet is not None:
                    facets.append(facet)
    return tuple(facets)


def edges_from_facets(
    facets: Sequence[tuple[int, int, int]],
) -> tuple[tuple[int, int], ...]:
    seen: set[tuple[int, int]] = set()
    for a, b, c in facets:
        for u, v in ((a, b), (b, c), (c, a)):
            seen.add((u, v) if u <= v else (v, u))
    return tuple(sorted(seen))
