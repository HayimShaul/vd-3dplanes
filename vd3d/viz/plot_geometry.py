"""Float-only drawing helpers for Phase 1 figures."""

from __future__ import annotations

import numpy as np

from vd3d.geometry.line2d import Line2D
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane
from vd3d.geometry.points import Point3D
from vd3d.viz.convert import to_float


def plane_mesh(plane: Plane, *, lim: float = 3.0, n: int = 16) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample ``plane`` as a surface over a cube of half-width ``lim``."""
    a, b, c, d = (to_float(plane.a), to_float(plane.b), to_float(plane.c), to_float(plane.d))
    grid = np.linspace(-lim, lim, n)
    aa, bb, cc = abs(a), abs(b), abs(c)
    if cc >= aa and cc >= bb:
        x, y = np.meshgrid(grid, grid)
        z = -(a * x + b * y + d) / c
        return x, y, z
    if bb >= aa:
        x, z = np.meshgrid(grid, grid)
        y = -(a * x + c * z + d) / b
        return x, y, z
    y, z = np.meshgrid(grid, grid)
    x = -(b * y + c * z + d) / a
    return x, y, z


def draw_plane(ax, plane: Plane, *, lim: float = 3.0, color: str, alpha: float = 0.35, label: str | None = None) -> None:
    x, y, z = plane_mesh(plane, lim=lim)
    ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0, antialiased=True)
    if label:
        # A single proxy scatter so the surface can appear in the legend.
        ax.scatter([], [], [], color=color, label=label)


def draw_line3d(
    ax,
    line: Line3D,
    *,
    t_min: float = -3.0,
    t_max: float = 3.0,
    color: str = "cyan",
    linewidth: float = 3.0,
    label: str | None = None,
) -> None:
    p0 = line.point_at(0)
    dx, dy, dz = (to_float(c) for c in line.direction)
    origin = np.array([to_float(p0.x), to_float(p0.y), to_float(p0.z)], dtype=float)
    direction = np.array([dx, dy, dz], dtype=float)
    start = origin + t_min * direction
    end = origin + t_max * direction
    ax.plot(
        [start[0], end[0]],
        [start[1], end[1]],
        [start[2], end[2]],
        color=color,
        linewidth=linewidth,
        label=label,
    )


def draw_point3d(
    ax,
    point: Point3D,
    *,
    color: str,
    size: float = 60,
    label: str | None = None,
    edgecolor: str = "black",
) -> None:
    ax.scatter(
        [to_float(point.x)],
        [to_float(point.y)],
        [to_float(point.z)],
        c=color,
        s=size,
        edgecolors=edgecolor,
        depthshade=False,
        label=label,
    )


def draw_lifted_line2d(
    ax,
    line: Line2D,
    z: int | float,
    *,
    extent: float = 4.0,
    color: str = "red",
    label: str | None = None,
) -> None:
    """Draw ``line`` at height ``z`` as a 3D segment."""
    p0, p1 = line.sample_points()
    dx = to_float(p1.x - p0.x)
    dy = to_float(p1.y - p0.y)
    length = (dx * dx + dy * dy) ** 0.5
    if length == 0:
        raise ValueError("cannot draw a degenerate 2D line")
    ux, uy = dx / length, dy / length
    x0, y0 = to_float(p0.x), to_float(p0.y)
    zf = float(z)
    ax.plot(
        [x0 - extent * ux, x0 + extent * ux],
        [y0 - extent * uy, y0 + extent * uy],
        [zf, zf],
        color=color,
        linewidth=3,
        label=label,
    )


def draw_line2d(ax, line: Line2D, *, lim: float = 5.0, color: str = "black", label: str | None = None) -> None:
    a, b, c = to_float(line.a), to_float(line.b), to_float(line.c)
    if abs(b) >= abs(a):
        xs = np.array([-lim, lim], dtype=float)
        ys = -(a * xs + c) / b
    else:
        ys = np.array([-lim, lim], dtype=float)
        xs = -(b * ys + c) / a
    ax.plot(xs, ys, color=color, linewidth=2, label=label)


def set_equal_3d(ax, lim: float = 3.0) -> None:
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
