"""Interactive GUI for the vertical decomposition of planes in R^3.

Uses the main-branch API: ``vd.vd(planes)`` returns prism cells as
6-tuples ``(x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil)``.
Unbounded prism cells are clipped to a viewing cube for display (the
arrangement itself is not boxed). The left pane is the 3D prisms; the
right pane is the xy trapezoidal map that the decomposition is stacked
over.

Run from the repo root:

    python example_vertical_decomposition_gui.py
    python example_vertical_decomposition_gui.py 4

The positional argument is the number of input planes (default 3).
All planes are random (reproducible with ``--seed``). Requires matplotlib
and scipy (``pip install matplotlib scipy``).

Drag to rotate the 3D view. Click the xy-map to highlight the prism stack
under that point.
"""

from __future__ import annotations

import argparse
import colorsys
import random
import sys
import tkinter as tk
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from tkinter import ttk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.figure import Figure
from matplotlib.patches import Polygon as MplPolygon
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
from scipy.spatial import ConvexHull, QhullError
from sympy import Point3D, Plane

import project
import vd

GEOM_EPS = 1e-8
VERTICAL_EPS = 1e-8


# ---------------------------------------------------------------------------
# Input
# ---------------------------------------------------------------------------

def _random_planes(n: int, seed: int) -> list[Plane]:
    """Random non-vertical, non-parallel planes, in the style of ``vd.test_vd``."""
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
        normal = np.array([float(c) for c in h.normal_vector], dtype=float)
        nrm = float(np.linalg.norm(normal))
        if nrm < 1e-12:
            continue
        # ``vd`` projects along z, so skip near-vertical planes.
        if abs(normal[2]) / nrm < 0.08:
            continue
        if any(h.is_parallel(p) for p in planes):
            continue
        planes.append(h)
    if len(planes) < n:
        raise RuntimeError(f"could only sample {len(planes)}/{n} random planes")
    return planes


def example_planes(n: int, seed: int = 2) -> list[Plane]:
    """``n`` random input planes."""
    if n < 1:
        raise ValueError("n must be at least 1")
    return _random_planes(n, seed)


# ---------------------------------------------------------------------------
# Geometry for drawing (clip unbounded cells to a finite viewing box)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ViewBox:
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    zmin: float
    zmax: float

    def corners(self) -> np.ndarray:
        xs = (self.xmin, self.xmax)
        ys = (self.ymin, self.ymax)
        zs = (self.zmin, self.zmax)
        return np.array([[x, y, z] for z in zs for y in ys for x in xs], dtype=float)


def _f(v) -> float:
    return float(v)


def line_y_at_x(line, x: float) -> float:
    """y-coordinate of an xy-line (``Line3D``) at the given x."""
    dx = _f(line.p2.x) - _f(line.p1.x)
    if abs(dx) < VERTICAL_EPS:
        return _f(line.p1.y)
    t = (x - _f(line.p1.x)) / dx
    return _f(line.p1.y) + t * (_f(line.p2.y) - _f(line.p1.y))


def _cell_interior_point(cell) -> np.ndarray:
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


def _sample_cell_points(cell) -> list[np.ndarray]:
    """A few points on the cell, used only to size the viewing cube."""
    pts: list[np.ndarray] = []
    try:
        xy = vd.find_center_point(cell)
    except Exception:
        return pts
    x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil = cell
    xs = [_f(xy.x)]
    if x_floor is not None:
        xs.append(_f(x_floor))
    if x_ceil is not None:
        xs.append(_f(x_ceil))
    for x in xs:
        ys = [_f(xy.y)]
        if y_floor is not None:
            try:
                ys.append(line_y_at_x(y_floor, x))
            except Exception:
                pass
        if y_ceil is not None:
            try:
                ys.append(line_y_at_x(y_ceil, x))
            except Exception:
                pass
        for y in ys:
            p_xy = Point3D(x, y, 0)
            if z_floor is not None:
                try:
                    p = project.project(p_xy, z_floor, "z")
                    pts.append(np.array([_f(p.x), _f(p.y), _f(p.z)], dtype=float))
                except Exception:
                    pass
            if z_ceil is not None:
                try:
                    p = project.project(p_xy, z_ceil, "z")
                    pts.append(np.array([_f(p.x), _f(p.y), _f(p.z)], dtype=float))
                except Exception:
                    pass
    return pts


def viewing_bbox(cells, scale: float = 1.0) -> ViewBox:
    """Cubic box around sample points of the cells."""
    pts: list[np.ndarray] = []
    for cell in cells:
        for p in _sample_cell_points(cell):
            if np.all(np.isfinite(p)) and np.max(np.abs(p)) < 1e6:
                pts.append(p)
    if not pts:
        h = 10.0 * max(scale, 0.25)
        return ViewBox(-h, h, -h, h, -h, h)
    arr = np.stack(pts)
    lo = arr.min(axis=0) - 0.75
    hi = arr.max(axis=0) + 0.75
    center = 0.5 * (lo + hi)
    half = 0.5 * float(np.max(hi - lo)) * max(scale, 0.25)
    half = max(half, 1.0)
    lo_c = center - half
    hi_c = center + half
    return ViewBox(float(lo_c[0]), float(hi_c[0]), float(lo_c[1]), float(hi_c[1]), float(lo_c[2]), float(hi_c[2]))


def _order_polygon(verts: np.ndarray, normal: np.ndarray) -> np.ndarray:
    if len(verts) <= 3:
        return verts
    n = np.asarray(normal, dtype=float)
    nrm = float(np.linalg.norm(n))
    if nrm < VERTICAL_EPS:
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


def _add_oriented_plane(planes_nd: list, normal, point, interior: np.ndarray) -> None:
    n = np.asarray(normal, dtype=float)
    nrm = float(np.linalg.norm(n))
    if nrm <= VERTICAL_EPS:
        return
    n = n / nrm
    p = np.asarray(point, dtype=float)
    # Halfspace n · x <= n · p contains ``interior``.
    if float(np.dot(n, interior - p)) > 0:
        n = -n
    planes_nd.append((n, float(np.dot(n, p))))


def _cell_halfspaces(cell, interior: np.ndarray) -> list[tuple[np.ndarray, float]]:
    x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil = cell
    hs: list[tuple[np.ndarray, float]] = []
    if x_floor is not None:
        _add_oriented_plane(hs, (-1.0, 0.0, 0.0), (_f(x_floor), 0.0, 0.0), interior)
    if x_ceil is not None:
        _add_oriented_plane(hs, (1.0, 0.0, 0.0), (_f(x_ceil), 0.0, 0.0), interior)
    if y_floor is not None:
        p1 = np.array([_f(y_floor.p1.x), _f(y_floor.p1.y), 0.0])
        p2 = np.array([_f(y_floor.p2.x), _f(y_floor.p2.y), 0.0])
        _add_oriented_plane(hs, np.cross(p2 - p1, (0.0, 0.0, 1.0)), p1, interior)
    if y_ceil is not None:
        p1 = np.array([_f(y_ceil.p1.x), _f(y_ceil.p1.y), 0.0])
        p2 = np.array([_f(y_ceil.p2.x), _f(y_ceil.p2.y), 0.0])
        _add_oriented_plane(hs, np.cross(p2 - p1, (0.0, 0.0, 1.0)), p1, interior)
    if z_floor is not None:
        n = np.array([_f(c) for c in z_floor.normal_vector], dtype=float)
        p = np.array([_f(z_floor.p1.x), _f(z_floor.p1.y), _f(z_floor.p1.z)])
        _add_oriented_plane(hs, n, p, interior)
    if z_ceil is not None:
        n = np.array([_f(c) for c in z_ceil.normal_vector], dtype=float)
        p = np.array([_f(z_ceil.p1.x), _f(z_ceil.p1.y), _f(z_ceil.p1.z)])
        _add_oriented_plane(hs, n, p, interior)
    return hs


def _bbox_halfspaces(box: ViewBox) -> list[tuple[np.ndarray, float]]:
    return [
        (np.array([-1.0, 0.0, 0.0]), -box.xmin),
        (np.array([1.0, 0.0, 0.0]), box.xmax),
        (np.array([0.0, -1.0, 0.0]), -box.ymin),
        (np.array([0.0, 1.0, 0.0]), box.ymax),
        (np.array([0.0, 0.0, -1.0]), -box.zmin),
        (np.array([0.0, 0.0, 1.0]), box.zmax),
    ]


def unique_points(pts: list[np.ndarray], eps: float = 10 * GEOM_EPS) -> list[np.ndarray]:
    if not pts:
        return []
    arr = np.stack(pts)
    keys = np.round(arr / eps)
    _, idx = np.unique(keys, axis=0, return_index=True)
    return [arr[i] for i in np.sort(idx)]


def _halfspace_vertices(planes_nd, eps: float = 10 * GEOM_EPS):
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


def _faces_from_points(points: list[np.ndarray]):
    """Convex polygonal faces and boundary edges of a point set."""
    arr = np.stack(points)
    if len(arr) < 3:
        return [], []
    if len(arr) == 3:
        return [arr], [(arr[0], arr[1]), (arr[1], arr[2]), (arr[2], arr[0])]
    try:
        hull = ConvexHull(arr)
    except (QhullError, ValueError):
        try:
            hull = ConvexHull(arr, qhull_options="QJ")
        except (QhullError, ValueError):
            return [], []

    groups: dict[tuple, set[int]] = {}
    for simplex in hull.simplices:
        p0, p1, p2 = arr[simplex[0]], arr[simplex[1]], arr[simplex[2]]
        n = np.cross(p1 - p0, p2 - p0)
        nrm = float(np.linalg.norm(n))
        if nrm < GEOM_EPS:
            continue
        n = n / nrm
        d = float(np.dot(n, p0))
        for component in n:
            if abs(component) > 1e-8:
                if component < 0:
                    n = -n
                    d = -d
                break
        key = tuple(np.round(n, 5)) + (round(d, 5),)
        groups.setdefault(key, set()).update(int(i) for i in simplex)

    faces: list[np.ndarray] = []
    edges: list[tuple[np.ndarray, np.ndarray]] = []
    for key, idxs in groups.items():
        ordered = _order_polygon(arr[list(idxs)], np.asarray(key[:3], dtype=float))
        if len(ordered) < 3:
            continue
        faces.append(ordered)
        for i in range(len(ordered)):
            edges.append((ordered[i], ordered[(i + 1) % len(ordered)]))
    return faces, edges


def clip_cell_mesh(cell, bbox: ViewBox):
    """Faces and edges of ``cell ∩ bbox`` (empty if the clip has no volume)."""
    try:
        interior = _cell_interior_point(cell)
    except Exception:
        return [], []
    hs = _cell_halfspaces(cell, interior) + _bbox_halfspaces(bbox)
    pts = _halfspace_vertices(hs)
    if len(pts) < 4:
        return [], []
    return _faces_from_points(pts)


def clip_trap_to_rect(cell, xmin: float, xmax: float, ymin: float, ymax: float):
    """xy-polygon of a cell's trapezoid clipped to an axis-aligned rectangle."""
    x_floor, x_ceil, y_floor, y_ceil = cell[0], cell[1], cell[2], cell[3]
    x_left = xmin if x_floor is None else max(xmin, _f(x_floor))
    x_right = xmax if x_ceil is None else min(xmax, _f(x_ceil))
    if x_right - x_left <= 10 * GEOM_EPS:
        return None

    def y_span(x: float) -> tuple[float, float] | None:
        y_lo = ymin if y_floor is None else line_y_at_x(y_floor, x)
        y_hi = ymax if y_ceil is None else line_y_at_x(y_ceil, x)
        y_lo = max(y_lo, ymin)
        y_hi = min(y_hi, ymax)
        if y_hi - y_lo <= 10 * GEOM_EPS:
            return None
        return y_lo, y_hi

    left = y_span(x_left)
    right = y_span(x_right)
    if left is None or right is None:
        return None
    return np.array(
        [
            [x_left, left[0]],
            [x_right, right[0]],
            [x_right, right[1]],
            [x_left, left[1]],
        ],
        dtype=float,
    )


def plane_polygon_in_box(plane: Plane, bbox: ViewBox) -> np.ndarray | None:
    """Convex polygon = plane ∩ bbox, or None if they miss."""
    n = np.array([_f(c) for c in plane.normal_vector], dtype=float)
    nrm = float(np.linalg.norm(n))
    if nrm < VERTICAL_EPS:
        return None
    n = n / nrm
    p0 = np.array([_f(plane.p1.x), _f(plane.p1.y), _f(plane.p1.z)])

    def signed_distance(pt: np.ndarray) -> float:
        return float(np.dot(n, pt - p0))

    corners = bbox.corners()
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
        a, c = corners[i], corners[j]
        da, dc = signed_distance(a), signed_distance(c)
        if abs(da) <= 10 * GEOM_EPS:
            pts.append(a)
        if abs(dc) <= 10 * GEOM_EPS:
            pts.append(c)
        if da * dc < -GEOM_EPS * GEOM_EPS:
            t = da / (da - dc)
            pts.append(a + t * (c - a))
    pts = unique_points(pts, eps=10 * GEOM_EPS)
    if len(pts) < 3:
        return None
    return _order_polygon(np.stack(pts), n)


def _cell_color(index: int, _n: int) -> tuple[float, float, float]:
    hue = (index * 0.61803398875) % 1.0
    return colorsys.hsv_to_rgb(hue, 0.55, 0.88)


def _box_edges(bbox: ViewBox) -> list[tuple[np.ndarray, np.ndarray]]:
    c = bbox.corners()
    pairs = (
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
    return [(c[i], c[j]) for i, j in pairs]


def trap_key(cell) -> tuple:
    x_floor, x_ceil, y_floor, y_ceil = cell[0], cell[1], cell[2], cell[3]

    def line_key(line) -> tuple | None:
        if line is None:
            return None
        return (
            round(_f(line.p1.x), 8),
            round(_f(line.p1.y), 8),
            round(_f(line.p2.x), 8),
            round(_f(line.p2.y), 8),
        )

    return (
        None if x_floor is None else round(_f(x_floor), 8),
        None if x_ceil is None else round(_f(x_ceil), 8),
        line_key(y_floor),
        line_key(y_ceil),
    )


def cell_contains_xy(cell, x: float, y: float) -> bool:
    x_floor, x_ceil, y_floor, y_ceil = cell[0], cell[1], cell[2], cell[3]
    if x_floor is not None and x < _f(x_floor) - 10 * GEOM_EPS:
        return False
    if x_ceil is not None and x > _f(x_ceil) + 10 * GEOM_EPS:
        return False
    if y_floor is not None and y < line_y_at_x(y_floor, x) - 10 * GEOM_EPS:
        return False
    if y_ceil is not None and y > line_y_at_x(y_ceil, x) + 10 * GEOM_EPS:
        return False
    return True


def z_on_plane(plane: Plane, x: float, y: float) -> float:
    return _f(project.project(Point3D(x, y, 0), plane, "z").z)


def cell_contains_z(cell, x: float, y: float, z: float) -> bool:
    z_floor, z_ceil = cell[4], cell[5]
    if z_floor is not None and z < z_on_plane(z_floor, x, y) - 10 * GEOM_EPS:
        return False
    if z_ceil is not None and z > z_on_plane(z_ceil, x, y) + 10 * GEOM_EPS:
        return False
    return True


def cell_is_unbounded(cell) -> bool:
    return any(bound is None for bound in cell)


def clipped_volume(faces) -> float:
    if not faces:
        return 0.0
    pts = unique_points([v for face in faces for v in face])
    if len(pts) < 4:
        return 0.0
    try:
        return float(ConvexHull(np.stack(pts)).volume)
    except (QhullError, ValueError):
        return 0.0


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

class VDViewer(tk.Tk):
    def __init__(self, planes: list[Plane]):
        super().__init__()
        self.title(f"Vertical decomposition of {len(planes)} planes in R³")
        self.geometry("1280x780")
        self.minsize(900, 560)

        self.planes = planes
        self.cells: list | None = None
        self.bbox: ViewBox | None = None
        self.meshes: list[tuple[list, list]] = []
        self.trap_repr: list = []
        self.trap_stacks: dict = {}
        self.selected: int = -1  # -1 = all cells
        self._draw_planes = tk.BooleanVar(value=True)
        self._draw_edges = tk.BooleanVar(value=True)
        self._building = False
        self._scale_job: str | None = None
        self._elev = 22.0
        self._azim = -60.0
        self._view_initialized = False

        self._build_chrome()
        self._build_figure()
        self.bind("<Left>", lambda _e: self._step_cell(-1))
        self.bind("<Right>", lambda _e: self._step_cell(+1))
        self.bind("a", lambda _e: self._select_cell(-1))
        self.protocol("WM_DELETE_WINDOW", self.destroy)
        self.after(50, self._compute)

    def _build_chrome(self) -> None:
        bar = ttk.Frame(self, padding=(8, 6))
        bar.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(bar, text="Cell").pack(side=tk.LEFT)
        ttk.Button(bar, text="All", width=4, command=lambda: self._select_cell(-1)).pack(
            side=tk.LEFT, padx=(4, 2)
        )
        ttk.Button(bar, text="◀", width=3, command=lambda: self._step_cell(-1)).pack(side=tk.LEFT)
        ttk.Button(bar, text="▶", width=3, command=lambda: self._step_cell(+1)).pack(
            side=tk.LEFT, padx=(0, 8)
        )

        ttk.Label(bar, text="Opacity").pack(side=tk.LEFT)
        self.opacity = tk.DoubleVar(value=0.35)
        ttk.Scale(
            bar,
            from_=0.05,
            to=0.9,
            variable=self.opacity,
            command=lambda _v: self._redraw(),
            length=120,
        ).pack(side=tk.LEFT, padx=(4, 12))

        ttk.Label(bar, text="View").pack(side=tk.LEFT)
        self.view_scale = tk.DoubleVar(value=1.0)
        ttk.Scale(
            bar,
            from_=0.5,
            to=2.0,
            variable=self.view_scale,
            command=self._on_view_scale,
            length=100,
        ).pack(side=tk.LEFT, padx=(4, 12))

        ttk.Checkbutton(
            bar, text="Input planes", variable=self._draw_planes, command=self._redraw
        ).pack(side=tk.LEFT, padx=4)
        ttk.Checkbutton(
            bar, text="Edges", variable=self._draw_edges, command=self._redraw
        ).pack(side=tk.LEFT, padx=4)

        self.status = ttk.Label(self, text="Building…", padding=(8, 4))
        self.status.pack(side=tk.BOTTOM, fill=tk.X)

        hint = ttk.Label(
            self,
            text="Drag to rotate · click the xy-map to pick a prism · ←/→ cycle cells · A shows all.  "
            "Cells recede to infinity; they are clipped to the viewing cube for display.",
            padding=(8, 0),
        )
        hint.pack(side=tk.BOTTOM, fill=tk.X)

    def _build_figure(self) -> None:
        holder = ttk.Frame(self)
        holder.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.fig = Figure(figsize=(12.0, 6.2), dpi=100, layout="constrained")
        self.ax3d = self.fig.add_subplot(1, 2, 1, projection="3d")
        self.ax2d = self.fig.add_subplot(1, 2, 2)

        self.canvas = FigureCanvasTkAgg(self.fig, master=holder)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas, holder, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)

        self.canvas.mpl_connect("button_press_event", self._on_click)

    def _compute(self) -> None:
        self.status.configure(text="Computing vertical decomposition…")
        self.update_idletasks()
        try:
            self.cells = vd.vd(self.planes)
        except Exception as exc:
            self.status.configure(text=f"vd() failed: {exc}")
            return
        self.selected = -1
        stacks: dict[tuple, list[int]] = defaultdict(list)
        order: list[tuple] = []
        reprs: dict[tuple, object] = {}
        for i, cell in enumerate(self.cells):
            key = trap_key(cell)
            if key not in stacks:
                order.append(key)
                reprs[key] = cell
            stacks[key].append(i)
        self.trap_stacks = stacks
        self.trap_repr = [(key, reprs[key]) for key in order]
        self._rebuild_meshes()

    def _on_view_scale(self, _value: str) -> None:
        if self._scale_job is not None:
            self.after_cancel(self._scale_job)
        self._scale_job = self.after(180, self._rebuild_meshes)

    def _rebuild_meshes(self, *_args) -> None:
        self._scale_job = None
        if self.cells is None or self._building:
            return
        self._building = True
        try:
            self.bbox = viewing_bbox(self.cells, scale=float(self.view_scale.get()))
            self.meshes = [clip_cell_mesh(cell, self.bbox) for cell in self.cells]
            self._redraw()
        finally:
            self._building = False

    def _select_cell(self, index: int) -> None:
        n = len(self.cells) if self.cells is not None else 0
        if n == 0:
            self.selected = -1
        elif index < 0:
            self.selected = -1
        else:
            self.selected = int(index) % n
        self._redraw()

    def _step_cell(self, delta: int) -> None:
        n = len(self.cells) if self.cells is not None else 0
        if n == 0:
            return
        if self.selected < 0:
            self.selected = 0 if delta > 0 else n - 1
        else:
            self.selected = (self.selected + delta) % n
        self._redraw()

    def _on_click(self, event) -> None:
        if event.inaxes is not self.ax2d or event.xdata is None or event.ydata is None:
            return
        if self.cells is None or self.bbox is None:
            return
        x, y = float(event.xdata), float(event.ydata)
        hits = [i for i, cell in enumerate(self.cells) if cell_contains_xy(cell, x, y)]
        if not hits:
            return
        z_samples = [
            0.5 * (self.bbox.zmin + self.bbox.zmax),
            *(
                (1 - t) * self.bbox.zmin + t * self.bbox.zmax
                for t in (0.15, 0.35, 0.65, 0.85)
            ),
        ]
        for z in z_samples:
            for i in hits:
                try:
                    if cell_contains_z(self.cells[i], x, y, z):
                        self._select_cell(i)
                        return
                except Exception:
                    continue
        self._select_cell(hits[0])

    def _redraw(self) -> None:
        if self.cells is None or self.bbox is None:
            return
        if self._view_initialized:
            self._elev = getattr(self.ax3d, "elev", self._elev)
            self._azim = getattr(self.ax3d, "azim", self._azim)
        self.ax3d.cla()
        self.ax2d.cla()
        self._draw_3d()
        self._draw_2d()
        self._update_status()
        self.canvas.draw_idle()

    def _draw_3d(self) -> None:
        ax = self.ax3d
        box = self.bbox
        n = len(self.cells)
        alpha = float(self.opacity.get())
        show_all = self.selected < 0

        face_verts: list[np.ndarray] = []
        face_colors: list[tuple] = []
        edge_segs: list = []
        edge_colors: list[tuple] = []
        edge_widths: list[float] = []

        for i, (faces, edges) in enumerate(self.meshes):
            rgb = _cell_color(i, n)
            highlighted = show_all or i == self.selected
            a = alpha if show_all else (0.82 if i == self.selected else 0.06)
            if not highlighted and not show_all:
                rgb_draw = (0.55, 0.55, 0.58)
            else:
                rgb_draw = rgb
            for face in faces:
                face_verts.append(face)
                face_colors.append((*rgb_draw, a))
            if self._draw_edges.get():
                ew = 1.15 if (show_all or i == self.selected) else 0.25
                ec = (0.12, 0.12, 0.14, 0.85) if (show_all or i == self.selected) else (
                    0.4,
                    0.4,
                    0.42,
                    0.25,
                )
                for p, q in edges:
                    edge_segs.append(np.stack([p, q]))
                    edge_colors.append(ec)
                    edge_widths.append(ew)

        if face_verts:
            coll = Poly3DCollection(
                face_verts, facecolors=face_colors, linewidths=0, shade=False
            )
            ax.add_collection3d(coll)
        if edge_segs:
            ax.add_collection3d(
                Line3DCollection(
                    edge_segs, colors=edge_colors, linewidths=edge_widths, zorder=3
                )
            )

        if self._draw_planes.get():
            plane_polys = []
            for h in self.planes:
                poly = plane_polygon_in_box(h, box)
                if poly is not None:
                    plane_polys.append(poly)
            if plane_polys:
                ax.add_collection3d(
                    Poly3DCollection(
                        plane_polys,
                        facecolors=(0.15, 0.15, 0.18, 0.07),
                        edgecolors=(0.05, 0.05, 0.08, 0.7),
                        linewidths=1.1,
                        linestyles="--",
                        shade=False,
                    )
                )

        ax.add_collection3d(
            Line3DCollection(
                [np.stack(seg) for seg in _box_edges(box)],
                colors=(0.35, 0.35, 0.38, 0.5),
                linewidths=0.6,
                linestyles=":",
            )
        )

        ax.set_xlim(box.xmin, box.xmax)
        ax.set_ylim(box.ymin, box.ymax)
        ax.set_zlim(box.zmin, box.zmax)
        try:
            ax.set_box_aspect(
                (
                    box.xmax - box.xmin,
                    box.ymax - box.ymin,
                    box.zmax - box.zmin,
                )
            )
        except (AttributeError, ValueError):
            pass
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_title("3D prisms (clipped to viewing cube)")
        ax.view_init(elev=self._elev, azim=self._azim)
        self._view_initialized = True

    def _draw_2d(self) -> None:
        ax = self.ax2d
        box = self.bbox
        n = len(self.cells)
        patches: list[MplPolygon] = []
        colors: list[tuple] = []
        outlines: list = []

        selected_key = None
        if 0 <= self.selected < n:
            selected_key = trap_key(self.cells[self.selected])

        for key, trap_cell in self.trap_repr:
            poly = clip_trap_to_rect(trap_cell, box.xmin, box.xmax, box.ymin, box.ymax)
            if poly is None or len(poly) < 3:
                continue
            outlines.append(np.vstack([poly, poly[0]]))
            stack = self.trap_stacks.get(key, [])
            if selected_key is not None and key == selected_key:
                rgb = _cell_color(self.selected, n)
                colors.append((*rgb, 0.75))
            elif self.selected < 0 and stack:
                rgb = _cell_color(stack[0], n)
                colors.append((*rgb, 0.28))
            else:
                colors.append((0.82, 0.84, 0.86, 0.35))
            patches.append(MplPolygon(poly, closed=True))

        if patches:
            ax.add_collection(PatchCollection(patches, facecolors=colors, edgecolors="none"))
        if outlines:
            ax.add_collection(
                LineCollection(outlines, colors=(0.15, 0.15, 0.18, 0.8), linewidths=0.7)
            )

        ax.set_xlim(box.xmin, box.xmax)
        ax.set_ylim(box.ymin, box.ymax)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("xy trapezoidal map")
        n_traps = len(self.trap_repr)
        ax.text(
            0.02,
            0.98,
            f"{n_traps} trapezoids  ·  {n} prisms",
            transform=ax.transAxes,
            va="top",
            fontsize=9,
            color="0.25",
        )

    def _update_status(self) -> None:
        n = len(self.cells)
        n_planes = len(self.planes)
        n_traps = len(self.trap_repr)
        if self.selected < 0:
            self.status.configure(
                text=f"{n_planes} planes  ·  {n} prism cells  ·  {n_traps} xy-trapezoids  ·  "
                "showing all cells (clipped)"
            )
            return
        cell = self.cells[self.selected]
        kind = "unbounded" if cell_is_unbounded(cell) else "bounded"
        faces, _edges = self.meshes[self.selected]
        n_verts = len(unique_points([v for face in faces for v in face])) if faces else 0
        vol = clipped_volume(faces)
        names = ("x−", "x+", "y−", "y+", "z−", "z+")
        missing = [name for name, bound in zip(names, cell) if bound is None]
        bounds = "all six bounds" if not missing else ("open " + ",".join(missing))
        self.status.configure(
            text=(
                f"Cell {self.selected}/{n}  ·  {kind}  ·  {bounds}  ·  "
                f"{n_verts} clipped vertices  ·  {len(faces)} faces  ·  "
                f"clipped volume {vol:.4f}"
            )
        )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Interactive 3D/2D viewer for the vertical decomposition of planes."
    )
    parser.add_argument(
        "n",
        nargs="?",
        type=int,
        default=3,
        help="number of input planes (default: 3)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2,
        help="RNG seed for the random planes (default: 2)",
    )
    args = parser.parse_args()
    if args.n < 1:
        parser.error("n must be at least 1")
    if args.n > 6:
        print(
            f"warning: n = {args.n} may be slow (this demo is intended for n ≲ 6)",
            file=sys.stderr,
        )
    planes = example_planes(args.n, seed=args.seed)
    app = VDViewer(planes)
    app.mainloop()


if __name__ == "__main__":
    main()
