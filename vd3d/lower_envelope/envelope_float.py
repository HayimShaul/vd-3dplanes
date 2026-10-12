"""Float lower envelope of planes via dual 3D convex hull (scipy).

Each non-vertical plane ``a x + b y + c z + d = 0`` (``c != 0``) is written
``z = A x + B y + C`` with ``A = -a/c``, ``B = -b/c``, ``C = -d/c`` and mapped
to the dual point ``(A, B, C)``. Lower-hull facets of those points are the
envelope vertices (triple intersections); hull vertices are the planes that
appear on the envelope.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.spatial import ConvexHull, QhullError

from vd3d.geometry.plane import Plane

LOWER_NORMAL_EPS = 1e-12


@dataclass(frozen=True, slots=True)
class FloatLowerEnvelope:
    """Lower envelope computed in floating point."""

    planes: tuple[Plane, ...]
    dual_points: np.ndarray  # (n, 3)
    facets: np.ndarray  # (m, 3) plane indices
    vertices: np.ndarray  # (m, 3) triple-intersection points
    edges: tuple[tuple[int, int], ...]
    active: np.ndarray  # unique plane indices on the envelope


def plane_to_dual_abc(plane: Plane) -> tuple[float, float, float]:
    """Map ``a x + b y + c z + d = 0`` to dual ``(A, B, C)`` for ``z = A x + B y + C``."""
    c = float(plane.c)
    if c == 0.0:
        raise ValueError(
            f"plane P{plane.id} is vertical (c = 0); lower-envelope duality requires c != 0"
        )
    return (-float(plane.a) / c, -float(plane.b) / c, -float(plane.d) / c)


def dual_points_from_planes(planes: Sequence[Plane]) -> np.ndarray:
    """Stack dual points; rejects vertical planes."""
    if not planes:
        raise ValueError("need at least one non-vertical plane")
    rows = [plane_to_dual_abc(plane) for plane in planes]
    return np.asarray(rows, dtype=float)


def _intersect_three_graphs(
    dual: np.ndarray, i: int, j: int, k: int
) -> np.ndarray | None:
    """Intersection of ``z = A x + B y + C`` for three dual rows, or ``None`` if parallel."""
    rows = dual[[i, j, k]]
    # A x + B y - z = -C
    matrix = np.array(
        [
            [rows[0, 0], rows[0, 1], -1.0],
            [rows[1, 0], rows[1, 1], -1.0],
            [rows[2, 0], rows[2, 1], -1.0],
        ],
        dtype=float,
    )
    rhs = -rows[:, 2]
    try:
        point = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError:
        return None
    if not np.all(np.isfinite(point)):
        return None
    return point


def _edges_from_facets(facets: np.ndarray) -> tuple[tuple[int, int], ...]:
    seen: set[tuple[int, int]] = set()
    for i, j, k in facets:
        for a, b in ((int(i), int(j)), (int(j), int(k)), (int(k), int(i))):
            seen.add((a, b) if a <= b else (b, a))
    return tuple(sorted(seen))


def _vertices_for_facets(
    dual: np.ndarray, facets: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    verts: list[np.ndarray] = []
    kept: list[list[int]] = []
    for triple in facets:
        i, j, k = (int(triple[0]), int(triple[1]), int(triple[2]))
        point = _intersect_three_graphs(dual, i, j, k)
        if point is None:
            continue
        verts.append(point)
        kept.append([i, j, k])
    if not kept:
        return np.zeros((0, 3), dtype=int), np.zeros((0, 3), dtype=float)
    return np.asarray(kept, dtype=int), np.stack(verts)


def _lower_facets_from_hull(hull: ConvexHull) -> np.ndarray:
    lower = hull.equations[:, 2] < -LOWER_NORMAL_EPS
    return np.asarray(hull.simplices[lower], dtype=int)


def _lower_facet_three_points(dual: np.ndarray) -> np.ndarray:
    """Single triangle on three dual points: keep the orientation with nz < 0."""
    if dual.shape[0] != 3:
        raise ValueError("expected exactly three dual points")
    i, j, k = 0, 1, 2
    a = dual[j] - dual[i]
    b = dual[k] - dual[i]
    normal = np.cross(a, b)
    if abs(float(normal[2])) <= LOWER_NORMAL_EPS:
        return np.zeros((0, 3), dtype=int)
    if float(normal[2]) < 0:
        return np.asarray([[i, j, k]], dtype=int)
    return np.asarray([[i, k, j]], dtype=int)


def lower_envelope_float(planes: Sequence[Plane]) -> FloatLowerEnvelope:
    """Compute the lower envelope of non-vertical planes using scipy's Qhull."""
    packed = tuple(planes)
    dual = dual_points_from_planes(packed)
    n = dual.shape[0]

    if n == 1:
        facets = np.zeros((0, 3), dtype=int)
        vertices = np.zeros((0, 3), dtype=float)
        active = np.array([0], dtype=int)
        edges: tuple[tuple[int, int], ...] = ()
    elif n == 2:
        facets = np.zeros((0, 3), dtype=int)
        vertices = np.zeros((0, 3), dtype=float)
        # Both planes appear on the envelope unless parallel as graphs with one above.
        active = np.array([0, 1], dtype=int)
        edges = ((0, 1),)
    elif n == 3:
        facets = _lower_facet_three_points(dual)
        facets, vertices = _vertices_for_facets(dual, facets)
        active = np.unique(facets) if len(facets) else np.arange(3, dtype=int)
        edges = _edges_from_facets(facets) if len(facets) else ((0, 1), (1, 2), (0, 2))
    else:
        try:
            hull = ConvexHull(dual)
        except QhullError as exc:
            raise ValueError(f"convex hull failed: {exc}") from exc
        facets = _lower_facets_from_hull(hull)
        facets, vertices = _vertices_for_facets(dual, facets)
        if len(facets) == 0:
            raise ValueError("lower hull has no facets (check plane orientations / inputs)")
        active = np.unique(facets)
        edges = _edges_from_facets(facets)

    return FloatLowerEnvelope(
        planes=packed,
        dual_points=dual,
        facets=facets,
        vertices=vertices,
        edges=edges,
        active=active,
    )


def format_float_envelope(env: FloatLowerEnvelope, *, seed: int | None = None) -> str:
    """Human-readable summary of a float lower envelope."""
    lines = [
        f"planes: {len(env.planes)}",
        f"active: {len(env.active)}",
        f"facets: {len(env.facets)}",
        f"edges: {len(env.edges)}",
    ]
    if seed is not None:
        lines.append(f"seed: {seed}")
    lines.extend(["", "input planes:"])
    for plane in env.planes:
        lines.append(
            f"  P{plane.id}: a={float(plane.a):g} b={float(plane.b):g} "
            f"c={float(plane.c):g} d={float(plane.d):g}"
        )
    lines.append("")
    lines.append("active plane indices (0-based into input list):")
    if len(env.active) == 0:
        lines.append("  (none)")
    else:
        for idx in env.active:
            plane = env.planes[int(idx)]
            lines.append(f"  [{int(idx)}] P{plane.id}")
    lines.append("")
    lines.append("envelope vertices (triple intersections):")
    if len(env.vertices) == 0:
        lines.append("  (none)")
    else:
        for facet, vertex in zip(env.facets, env.vertices, strict=True):
            i, j, k = (int(facet[0]), int(facet[1]), int(facet[2]))
            ids = (env.planes[i].id, env.planes[j].id, env.planes[k].id)
            lines.append(
                f"  planes {ids} -> "
                f"({vertex[0]:.6g}, {vertex[1]:.6g}, {vertex[2]:.6g})"
            )
    return "\n".join(lines) + "\n"
