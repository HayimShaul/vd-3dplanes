"""Combinatorial and volume checks for the 3D vertical decomposition."""

from __future__ import annotations

import random
from itertools import combinations

import numpy as np
import pytest
from scipy.spatial import ConvexHull, QhullError
from sympy import Point3D, Plane

import project
import vd

GEOM_EPS = 1e-8
VERTICAL_EPS = 1e-8
SEED = 7


def _f(v) -> float:
    return float(v)


def random_planes(n: int, seed: int) -> list[Plane]:
    """Non-vertical, pairwise non-parallel planes, as required by ``vd.vd``."""
    rng = random.Random(seed)
    planes: list[Plane] = []
    for _ in range(400 * max(n, 1)):
        if len(planes) >= n:
            break
        p1 = Point3D(rng.uniform(-10, 10), rng.uniform(-10, 10), rng.uniform(-10, 10))
        p2 = Point3D(rng.uniform(-10, 10), rng.uniform(-10, 10), rng.uniform(-10, 10))
        p3 = Point3D(rng.uniform(-10, 10), rng.uniform(-10, 10), rng.uniform(-10, 10))
        try:
            h = Plane(p1, p2, p3)
        except Exception:
            continue
        normal = np.array([_f(c) for c in h.normal_vector], dtype=float)
        nrm = float(np.linalg.norm(normal))
        if nrm < 1e-12:
            continue
        if abs(normal[2]) / nrm < 0.08:
            continue
        if any(h.is_parallel(p) for p in planes):
            continue
        planes.append(h)
    if len(planes) < n:
        raise RuntimeError(f"could only sample {len(planes)}/{n} random planes")
    return planes


def cell_interior_point(cell) -> np.ndarray:
    xy = vd.find_center_point(cell)
    z_floor, z_ceil = cell[4], cell[5]
    zs: list[float] = []
    if z_floor is not None:
        zs.append(_f(project.project(xy, z_floor, "z").z))
    if z_ceil is not None:
        zs.append(_f(project.project(xy, z_ceil, "z").z))
    if len(zs) == 2:
        z = 0.5 * (zs[0] + zs[1])
    elif z_floor is not None:
        z = zs[0] + 1.0
    elif z_ceil is not None:
        z = zs[0] - 1.0
    else:
        z = 0.0
    return np.array([_f(xy.x), _f(xy.y), z], dtype=float)


def _orient_plane(normal, point, interior: np.ndarray) -> tuple[np.ndarray, float] | None:
    n = np.asarray(normal, dtype=float)
    nrm = float(np.linalg.norm(n))
    if nrm <= VERTICAL_EPS:
        return None
    n = n / nrm
    p = np.asarray(point, dtype=float)
    if float(np.dot(n, interior - p)) > 0:
        n = -n
    return n, float(np.dot(n, p))


def cell_halfspaces(cell, interior: np.ndarray) -> list[tuple[np.ndarray, float]]:
    """Inward halfspaces ``n · x <= d`` of a 6-tuple prism cell."""
    x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil = cell
    hs: list[tuple[np.ndarray, float]] = []
    if x_floor is not None:
        plane = _orient_plane((-1.0, 0.0, 0.0), (_f(x_floor), 0.0, 0.0), interior)
        if plane is not None:
            hs.append(plane)
    if x_ceil is not None:
        plane = _orient_plane((1.0, 0.0, 0.0), (_f(x_ceil), 0.0, 0.0), interior)
        if plane is not None:
            hs.append(plane)
    if y_floor is not None:
        p1 = np.array([_f(y_floor.p1.x), _f(y_floor.p1.y), 0.0])
        p2 = np.array([_f(y_floor.p2.x), _f(y_floor.p2.y), 0.0])
        plane = _orient_plane(np.cross(p2 - p1, (0.0, 0.0, 1.0)), p1, interior)
        if plane is not None:
            hs.append(plane)
    if y_ceil is not None:
        p1 = np.array([_f(y_ceil.p1.x), _f(y_ceil.p1.y), 0.0])
        p2 = np.array([_f(y_ceil.p2.x), _f(y_ceil.p2.y), 0.0])
        plane = _orient_plane(np.cross(p2 - p1, (0.0, 0.0, 1.0)), p1, interior)
        if plane is not None:
            hs.append(plane)
    if z_floor is not None:
        n = np.array([_f(c) for c in z_floor.normal_vector], dtype=float)
        p = np.array([_f(z_floor.p1.x), _f(z_floor.p1.y), _f(z_floor.p1.z)])
        plane = _orient_plane(n, p, interior)
        if plane is not None:
            hs.append(plane)
    if z_ceil is not None:
        n = np.array([_f(c) for c in z_ceil.normal_vector], dtype=float)
        p = np.array([_f(z_ceil.p1.x), _f(z_ceil.p1.y), _f(z_ceil.p1.z)])
        plane = _orient_plane(n, p, interior)
        if plane is not None:
            hs.append(plane)
    return hs


def bbox_halfspaces(lo: np.ndarray, hi: np.ndarray) -> list[tuple[np.ndarray, float]]:
    return [
        (np.array([-1.0, 0.0, 0.0]), -float(lo[0])),
        (np.array([1.0, 0.0, 0.0]), float(hi[0])),
        (np.array([0.0, -1.0, 0.0]), -float(lo[1])),
        (np.array([0.0, 1.0, 0.0]), float(hi[1])),
        (np.array([0.0, 0.0, -1.0]), -float(lo[2])),
        (np.array([0.0, 0.0, 1.0]), float(hi[2])),
    ]


def unique_points(pts: list[np.ndarray], eps: float = 10 * GEOM_EPS) -> list[np.ndarray]:
    if not pts:
        return []
    arr = np.stack(pts)
    _, idx = np.unique(np.round(arr / eps), axis=0, return_index=True)
    return [arr[i] for i in np.sort(idx)]


def halfspace_vertices(planes_nd, eps: float = 10 * GEOM_EPS) -> list[np.ndarray]:
    pts: list[np.ndarray] = []
    for (n1, d1), (n2, d2), (n3, d3) in combinations(planes_nd, 3):
        N = np.stack([n1, n2, n3])
        try:
            p = np.linalg.solve(N, [d1, d2, d3])
        except np.linalg.LinAlgError:
            continue
        if not np.all(np.isfinite(p)):
            continue
        if all(float(np.dot(n, p)) <= d + eps for n, d in planes_nd):
            pts.append(p)
    return unique_points(pts, eps=eps)


def line_y_at_x(line, x: float) -> float:
    dx = _f(line.p2.x) - _f(line.p1.x)
    if abs(dx) < 1e-6:
        return _f(line.p1.y)
    t = (x - _f(line.p1.x)) / dx
    return _f(line.p1.y) + t * (_f(line.p2.y) - _f(line.p1.y))


def z_on_plane(plane, x: float, y: float) -> float:
    items = {str(k): float(v) for k, v in plane.equation().as_coefficients_dict().items()}
    a = items.get("x", 0.0)
    b = items.get("y", 0.0)
    c = items.get("z", 0.0)
    d = items.get("1", 0.0)
    if abs(c) < VERTICAL_EPS:
        raise ValueError("vertical plane has no z-graph")
    return -(a * x + b * y + d) / c


def combinatorial_corners(cell) -> dict[str, list[np.ndarray]]:
    """Trapezoidal-prism corners, grouped by the supporting bound they lie on.

    A cell is a prism over an xy trapezoid: at most two x-walls, two y-walls,
    one floor and one ceiling. Each of those faces therefore has at most four
    corners (x ∈ {x_floor, x_ceil} × y-line ∈ {y_floor, y_ceil} × z-plane).
    """
    x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil = cell
    xs: list[tuple[str, float]] = []
    if x_floor is not None:
        xs.append(("x_floor", _f(x_floor)))
    if x_ceil is not None:
        xs.append(("x_ceil", _f(x_ceil)))
    ys: list[tuple[str, object]] = []
    if y_floor is not None:
        ys.append(("y_floor", y_floor))
    if y_ceil is not None:
        ys.append(("y_ceil", y_ceil))
    zs: list[tuple[str, object]] = []
    if z_floor is not None:
        zs.append(("z_floor", z_floor))
    if z_ceil is not None:
        zs.append(("z_ceil", z_ceil))

    faces: dict[str, list[np.ndarray]] = {name: [] for name, _ in xs + ys + zs}
    for x_name, x in xs:
        for y_name, y_line in ys:
            y = line_y_at_x(y_line, x)
            for z_name, z_plane in zs:
                z = z_on_plane(z_plane, x, y)
                p = np.array([x, y, z], dtype=float)
                if not np.all(np.isfinite(p)):
                    continue
                faces[x_name].append(p)
                faces[y_name].append(p)
                faces[z_name].append(p)
    return {name: unique_points(pts) for name, pts in faces.items() if pts}


def convex_volume(points: list[np.ndarray]) -> float:
    if len(points) < 4:
        return 0.0
    try:
        return float(ConvexHull(np.stack(points)).volume)
    except (QhullError, ValueError):
        return 0.0


def clipped_cell_volume(cell, lo: np.ndarray, hi: np.ndarray) -> float:
    interior = cell_interior_point(cell)
    hs = cell_halfspaces(cell, interior) + bbox_halfspaces(lo, hi)
    return convex_volume(halfspace_vertices(hs))


def finite_vertices(cells) -> np.ndarray:
    """Finite corners of the decomposition, plus one interior sample per cell."""
    verts: list[np.ndarray] = []
    for cell in cells:
        for pts in combinatorial_corners(cell).values():
            verts.extend(pts)
        try:
            interior = cell_interior_point(cell)
        except Exception:
            continue
        if np.all(np.isfinite(interior)):
            verts.append(interior)
    if not verts:
        raise AssertionError("decomposition has no finite vertices")
    arr = np.stack(verts)
    arr = arr[np.all(np.isfinite(arr), axis=1)]
    # Near-vertical y-walls can evaluate to huge coordinates; drop those outliers.
    if len(arr) == 0:
        raise AssertionError("decomposition has no finite vertices")
    scale = np.percentile(np.abs(arr), 80, axis=0)
    scale = np.maximum(scale, 1.0)
    kept = arr[np.all(np.abs(arr) <= 20.0 * scale, axis=1)]
    return kept if len(kept) else arr


def measurement_bbox(vertices: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lo = vertices.min(axis=0)
    hi = vertices.max(axis=0)
    pad = 0.5
    lo = lo - pad
    hi = hi + pad
    extent = np.maximum(hi - lo, 1e-9)
    return lo, lo + extent


@pytest.fixture(scope="module")
def n_planes(request):
    n = request.config.getoption("--n-planes")
    if n < 1:
        raise pytest.UsageError("--n-planes must be at least 1")
    return n


@pytest.fixture(scope="module")
def decomposition(n_planes):
    planes = random_planes(n_planes, SEED)
    assert len(planes) == n_planes
    cells = vd.vd(planes)
    assert cells, "vertical decomposition produced no cells"
    return cells


def test_cells_have_trapezoidal_prism_faces(decomposition):
    for cell in decomposition:
        x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil = cell
        n_floor = 0 if z_floor is None else 1
        n_ceil = 0 if z_ceil is None else 1
        n_walls = sum(bound is not None for bound in (x_floor, x_ceil, y_floor, y_ceil))
        assert n_floor <= 1
        assert n_ceil <= 1
        assert n_walls <= 4
        for face_pts in combinatorial_corners(cell).values():
            assert len(face_pts) <= 4


def test_clipped_volumes_sum_to_bounding_box(decomposition):
    cells = decomposition
    lo, hi = measurement_bbox(finite_vertices(cells))
    box_volume = float(np.prod(hi - lo))
    clipped_sum = sum(clipped_cell_volume(cell, lo, hi) for cell in cells)
    assert clipped_sum == pytest.approx(box_volume, rel=1e-2, abs=1e-6)
