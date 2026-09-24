"""Clip planes, lines, and 3D VD cells to an axis-aligned viewing box.

Float-only helpers for the interactive GUI. Kernel packages must not import this.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

from vd3d.cells3d.types import Cell3D
from vd3d.events.lines import compute_intersection_lines
from vd3d.events.slice import compute_vd_at_z
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane
from vd3d.geometry.points import Point3D
from vd3d.geometry.scalar import Scalar, as_scalar
from vd3d.sweep.algorithm import SweepResult
from vd3d.sweep.matching import cell_signature
from vd3d.viz.convert import to_float

GEOM_EPS = 1e-8


@dataclass(frozen=True)
class ViewBox:
    """Axis-aligned cube ``[-half, half]^3`` centered at the origin."""

    half: float

    def __post_init__(self) -> None:
        if self.half <= 0:
            raise ValueError("half must be positive")

    @property
    def lim(self) -> float:
        return self.half

    def corners(self) -> np.ndarray:
        h = self.half
        xs = (-h, h)
        ys = (-h, h)
        zs = (-h, h)
        return np.array([[x, y, z] for z in zs for y in ys for x in xs], dtype=float)

    def halfspaces(self) -> list[tuple[np.ndarray, float]]:
        h = self.half
        return [
            (np.array([-1.0, 0.0, 0.0]), h),
            (np.array([1.0, 0.0, 0.0]), h),
            (np.array([0.0, -1.0, 0.0]), h),
            (np.array([0.0, 1.0, 0.0]), h),
            (np.array([0.0, 0.0, -1.0]), h),
            (np.array([0.0, 0.0, 1.0]), h),
        ]


def default_view_half(result: SweepResult, *, minimum: float = 3.0) -> float:
    """A cube large enough to cover finite event heights with a little margin."""
    values = [minimum]
    for event in result.events:
        values.append(abs(to_float(event.z)) + 1.0)
    for cell in result.cells:
        if cell.lower_z is not None:
            values.append(abs(to_float(cell.lower_z)) + 1.0)
        if cell.upper_z is not None:
            values.append(abs(to_float(cell.upper_z)) + 1.0)
    return float(max(values))


def mid_interval_z(lower: Scalar | None, upper: Scalar | None) -> Scalar:
    if lower is None and upper is None:
        return as_scalar(0)
    if lower is None:
        return upper - 1  # type: ignore[operator]
    if upper is None:
        return lower + 1
    return (lower + upper) / 2


def interior_point_of_cell(result: SweepResult, cell_id: int) -> np.ndarray | None:
    """A point in the open 3D cell, or ``None`` if none can be reconstructed."""
    for interval in result.intervals:
        for sig, cid in interval.binding:
            if cid != cell_id:
                continue
            z = mid_interval_z(interval.lower_z, interval.upper_z)
            try:
                vd = compute_vd_at_z(result.planes, z)
            except (ValueError, RuntimeError):
                continue
            for cell2d in vd.cells:
                if cell_signature(vd, cell2d) == sig:
                    return np.array(
                        [
                            to_float(cell2d.representative.x),
                            to_float(cell2d.representative.y),
                            to_float(z),
                        ],
                        dtype=float,
                    )
    return None


def _oriented_halfspace(
    normal: Sequence[float], offset: float, interior: np.ndarray
) -> tuple[np.ndarray, float] | None:
    """``n · x <= d`` containing ``interior``. ``offset`` is ``n · p`` for a point on the plane."""
    n = np.asarray(normal, dtype=float)
    nrm = float(np.linalg.norm(n))
    if nrm < GEOM_EPS:
        return None
    n = n / nrm
    d = offset / nrm
    if float(np.dot(n, interior)) > d + GEOM_EPS:
        n = -n
        d = -d
    return n, d


def plane_halfspace(plane: Plane, interior: np.ndarray) -> tuple[np.ndarray, float] | None:
    a, b, c, d0 = (to_float(plane.a), to_float(plane.b), to_float(plane.c), to_float(plane.d))
    # a x + b y + c z = -d0
    return _oriented_halfspace((a, b, c), -d0, interior)


def cell_halfspaces(cell: Cell3D, interior: np.ndarray) -> list[tuple[np.ndarray, float]]:
    """Oriented halfspaces of the cell (floor, ceiling, walls, finite z-extent)."""
    out: list[tuple[np.ndarray, float]] = []
    for plane in (cell.floor, cell.ceiling, *cell.vertical_walls):
        if plane is None:
            continue
        hs = plane_halfspace(plane, interior)
        if hs is not None:
            out.append(hs)
    if cell.lower_z is not None:
        z = to_float(cell.lower_z)
        hs = _oriented_halfspace((0.0, 0.0, -1.0), -z, interior)
        if hs is not None:
            out.append(hs)
    if cell.upper_z is not None:
        z = to_float(cell.upper_z)
        hs = _oriented_halfspace((0.0, 0.0, 1.0), z, interior)
        if hs is not None:
            out.append(hs)
    return out


def unique_points(pts: list[np.ndarray], eps: float = 10 * GEOM_EPS) -> list[np.ndarray]:
    if not pts:
        return []
    arr = np.stack(pts)
    keys = np.round(arr / eps)
    _, idx = np.unique(keys, axis=0, return_index=True)
    return [arr[i] for i in np.sort(idx)]


def halfspace_vertices(
    planes_nd: Sequence[tuple[np.ndarray, float]], *, eps: float = 10 * GEOM_EPS
) -> list[np.ndarray]:
    """Vertices of a convex polyhedron given as ``n · x <= d`` halfspaces."""
    pts: list[np.ndarray] = []
    for (n1, d1), (n2, d2), (n3, d3) in combinations(planes_nd, 3):
        matrix = np.stack([n1, n2, n3])
        try:
            point = np.linalg.solve(matrix, [d1, d2, d3])
        except np.linalg.LinAlgError:
            continue
        if not np.all(np.isfinite(point)):
            continue
        if all(float(np.dot(n, point)) <= d + eps for n, d in planes_nd):
            pts.append(point)
    return unique_points(pts, eps=eps)


def order_polygon(verts: np.ndarray, normal: np.ndarray) -> np.ndarray:
    if len(verts) <= 3:
        return verts
    n = np.asarray(normal, dtype=float)
    nrm = float(np.linalg.norm(n))
    if nrm < GEOM_EPS:
        return verts
    n = n / nrm
    tangent = np.cross(n, np.array([0.0, 0.0, 1.0]))
    if float(np.linalg.norm(tangent)) < 1e-8:
        tangent = np.cross(n, np.array([0.0, 1.0, 0.0]))
    tangent = tangent / float(np.linalg.norm(tangent))
    bitangent = np.cross(n, tangent)
    rel = verts - verts.mean(axis=0)
    angles = np.arctan2(rel @ bitangent, rel @ tangent)
    return verts[np.argsort(angles)]


def faces_from_halfspaces(
    halfspaces: Sequence[tuple[np.ndarray, float]],
    vertices: Sequence[np.ndarray],
    *,
    eps: float | None = None,
) -> list[np.ndarray]:
    """Convex polygonal faces: vertices lying on each supporting halfspace."""
    if len(vertices) < 3:
        return []
    arr = np.stack(vertices)
    scale = float(np.max(np.abs(arr))) if len(arr) else 1.0
    tol = eps if eps is not None else max(20 * GEOM_EPS, 1e-6 * max(scale, 1.0))
    faces: list[np.ndarray] = []
    for normal, offset in halfspaces:
        on_face = arr[np.abs(arr @ normal - offset) <= tol]
        on_face = unique_points([p for p in on_face], eps=tol)
        if len(on_face) < 3:
            continue
        faces.append(order_polygon(np.stack(on_face), normal))
    return faces


def clipped_cell_faces(
    result: SweepResult, cell: Cell3D, box: ViewBox
) -> list[np.ndarray]:
    """Faces of ``cell ∩ box``, or empty if the cell misses the box."""
    interior = interior_point_of_cell(result, cell.id)
    if interior is None:
        return []
    hs = cell_halfspaces(cell, interior) + box.halfspaces()
    verts = halfspace_vertices(hs)
    if len(verts) < 3:
        return []
    return faces_from_halfspaces(hs, verts)


def plane_polygon_in_box(plane: Plane, box: ViewBox) -> np.ndarray | None:
    """Convex polygon ``plane ∩ box``, or ``None`` if they miss."""
    a, b, c, d0 = (to_float(plane.a), to_float(plane.b), to_float(plane.c), to_float(plane.d))
    n = np.array([a, b, c], dtype=float)
    nrm = float(np.linalg.norm(n))
    if nrm < GEOM_EPS:
        return None
    n = n / nrm
    # n · x = -d0 / nrm
    target = -d0 / nrm

    def signed(pt: np.ndarray) -> float:
        return float(np.dot(n, pt) - target)

    corners = box.corners()
    edges = (
        (0, 1),
        (0, 2),
        (0, 4),
        (1, 3),
        (1, 5),
        (2, 3),
        (2, 6),
        (3, 7),
        (4, 5),
        (4, 6),
        (5, 7),
        (6, 7),
    )
    pts: list[np.ndarray] = []
    for i, j in edges:
        p, q = corners[i], corners[j]
        sp, sq = signed(p), signed(q)
        if abs(sp) <= 10 * GEOM_EPS:
            pts.append(p)
        if abs(sq) <= 10 * GEOM_EPS:
            pts.append(q)
        if sp * sq < -GEOM_EPS * GEOM_EPS:
            t = sp / (sp - sq)
            pts.append(p + t * (q - p))
    pts = unique_points(pts)
    if len(pts) < 3:
        return None
    return order_polygon(np.stack(pts), n)


def clip_line_to_box(line: Line3D, box: ViewBox) -> tuple[np.ndarray, np.ndarray] | None:
    """Segment of ``line`` inside ``box``, or ``None`` if it misses."""
    origin = np.array(
        [to_float(line.point.x), to_float(line.point.y), to_float(line.point.z)],
        dtype=float,
    )
    direction = np.array([to_float(c) for c in line.direction], dtype=float)
    nrm = float(np.linalg.norm(direction))
    if nrm < GEOM_EPS:
        return None
    direction = direction / nrm
    h = box.half
    t_enter, t_exit = -np.inf, np.inf
    for axis in range(3):
        o = origin[axis]
        d = direction[axis]
        if abs(d) < GEOM_EPS:
            if o < -h - GEOM_EPS or o > h + GEOM_EPS:
                return None
            continue
        t0 = (-h - o) / d
        t1 = (h - o) / d
        t_lo, t_hi = (t0, t1) if t0 <= t1 else (t1, t0)
        t_enter = max(t_enter, t_lo)
        t_exit = min(t_exit, t_hi)
        if t_enter > t_exit:
            return None
    if not np.isfinite(t_enter) or not np.isfinite(t_exit):
        return None
    return origin + t_enter * direction, origin + t_exit * direction


def draw_planes_grey(ax, planes: Sequence[Plane], box: ViewBox) -> None:
    for plane in planes:
        poly = plane_polygon_in_box(plane, box)
        if poly is None:
            continue
        ax.add_collection3d(
            Poly3DCollection(
                [poly],
                facecolors=(0.55, 0.55, 0.55, 0.18),
                edgecolors=(0.35, 0.35, 0.35, 0.35),
                linewidths=0.4,
            )
        )


def draw_intersection_lines(ax, planes: Sequence[Plane], box: ViewBox) -> None:
    segments = []
    for line in compute_intersection_lines(planes):
        clipped = clip_line_to_box(line, box)
        if clipped is None:
            continue
        segments.append(clipped)
    if not segments:
        return
    ax.add_collection3d(
        Line3DCollection(segments, colors=(0.2, 0.2, 0.2, 1.0), linewidths=1.8)
    )


def draw_cell_red(ax, result: SweepResult, cell: Cell3D, box: ViewBox) -> None:
    faces = clipped_cell_faces(result, cell, box)
    if not faces:
        return
    ax.add_collection3d(
        Poly3DCollection(
            faces,
            facecolors=(0.85, 0.15, 0.12, 0.35),
            edgecolors=(0.55, 0.05, 0.05, 0.85),
            linewidths=1.2,
        )
    )


def draw_bbox_wire(ax, box: ViewBox) -> None:
    corners = box.corners()
    edges = (
        (0, 1),
        (0, 2),
        (0, 4),
        (1, 3),
        (1, 5),
        (2, 3),
        (2, 6),
        (3, 7),
        (4, 5),
        (4, 6),
        (5, 7),
        (6, 7),
    )
    segments = [(corners[i], corners[j]) for i, j in edges]
    ax.add_collection3d(
        Line3DCollection(segments, colors=(0.7, 0.7, 0.7, 0.6), linewidths=0.8, linestyles="--")
    )


def point3d_array(point: Point3D) -> np.ndarray:
    return np.array([to_float(point.x), to_float(point.y), to_float(point.z)], dtype=float)
