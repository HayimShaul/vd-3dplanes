"""Float-only drawing helpers for 2D vertical decompositions."""

from __future__ import annotations

from matplotlib.patches import Polygon

from vd3d.geometry.line2d import Line2D
from vd3d.vertical_decomposition.types import Hit, VDCell2D, VerticalDecomposition, VerticalRay
from vd3d.viz.convert import to_float
from vd3d.viz.plot_arrangement import (
    face_color,
    viewing_limit,
)


def vd_viewing_limit(vd: VerticalDecomposition, *, margin: float = 1.8, minimum: float = 2.8) -> int:
    return viewing_limit(vd.arrangement, margin=margin, minimum=minimum)


def draw_cells_filled(ax, vd: VerticalDecomposition, *, lim: float) -> None:
    n = len(vd.cells)
    for cell in vd.cells:
        poly = cell_clip_polygon(vd, cell, lim)
        if len(poly) < 3:
            continue
        ax.add_patch(
            Polygon(
                poly,
                closed=True,
                facecolor=face_color(cell.id, n),
                edgecolor="none",
                alpha=0.9,
                zorder=1,
            )
        )
        cx = sum(p[0] for p in poly) / len(poly)
        cy = sum(p[1] for p in poly) / len(poly)
        ax.text(
            cx,
            cy,
            str(cell.id),
            ha="center",
            va="center",
            fontsize=9,
            zorder=4,
        )


def draw_vertical_walls(ax, vd: VerticalDecomposition, *, lim: float) -> None:
    for wall in vd.walls:
        x = to_float(wall.x)
        y0 = -lim if wall.y_min is None else max(-lim, to_float(wall.y_min))
        y1 = lim if wall.y_max is None else min(lim, to_float(wall.y_max))
        if y1 <= y0:
            continue
        ax.plot(
            [x, x],
            [y0, y1],
            color="crimson",
            linestyle="--",
            linewidth=1.8,
            zorder=3,
        )


def draw_vertical_rays(
    ax,
    rays_and_hits: list[tuple[VerticalRay, Hit]],
    *,
    lim: float,
) -> None:
    for ray, hit in rays_and_hits:
        x = to_float(ray.origin.x)
        y0 = to_float(ray.origin.y)
        if hit.unbounded:
            y1 = lim if ray.direction > 0 else -lim
        else:
            y1 = to_float(hit.point.y)  # type: ignore[union-attr]
            ax.scatter([x], [y1], c="red", s=55, zorder=6, edgecolors="black", linewidths=0.4)
        ax.plot([x, x], [y0, y1], color="cyan", linewidth=2.2, zorder=4)
        ax.annotate(
            "+y" if ray.direction > 0 else "-y",
            (x, (y0 + y1) / 2),
            textcoords="offset points",
            xytext=(6, 0),
            fontsize=8,
            color="teal",
            zorder=5,
        )


def cell_clip_polygon(
    vd: VerticalDecomposition, cell: VDCell2D, lim: float
) -> list[tuple[float, float]]:
    """``cell ∩ [-lim, lim]^2`` as a convex polygon, in CCW order."""
    square = [(-lim, -lim), (lim, -lim), (lim, lim), (-lim, lim)]
    poly = square
    if cell.left_x is not None:
        poly = _clip_halfplane(poly, 1.0, 0.0, -to_float(cell.left_x))
    if cell.right_x is not None:
        poly = _clip_halfplane(poly, -1.0, 0.0, to_float(cell.right_x))
    lines = vd.arrangement.lines
    if cell.lower_line is not None:
        poly = _clip_above_line(poly, lines[cell.lower_line])
    if cell.upper_line is not None:
        poly = _clip_below_line(poly, lines[cell.upper_line])
    return poly


def _clip_above_line(poly: list[tuple[float, float]], line: Line2D) -> list[tuple[float, float]]:
    a, b, c = to_float(line.a), to_float(line.b), to_float(line.c)
    # y >= y_line(x)  iff  (a x + b y + c) * b >= 0
    return _clip_halfplane(poly, a * b, b * b, c * b)


def _clip_below_line(poly: list[tuple[float, float]], line: Line2D) -> list[tuple[float, float]]:
    a, b, c = to_float(line.a), to_float(line.b), to_float(line.c)
    return _clip_halfplane(poly, -a * b, -b * b, -c * b)


def _clip_halfplane(
    poly: list[tuple[float, float]], a: float, b: float, c: float
) -> list[tuple[float, float]]:
    """Keep ``a x + b y + c >= 0`` (Sutherland–Hodgman)."""

    def inside(p: tuple[float, float]) -> bool:
        return a * p[0] + b * p[1] + c >= -1e-12

    def intersect(
        p: tuple[float, float], q: tuple[float, float]
    ) -> tuple[float, float]:
        dp = a * p[0] + b * p[1] + c
        dq = a * q[0] + b * q[1] + c
        t = dp / (dp - dq)
        return (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))

    if not poly:
        return []
    out: list[tuple[float, float]] = []
    prev = poly[-1]
    prev_in = inside(prev)
    for cur in poly:
        cur_in = inside(cur)
        if cur_in:
            if not prev_in:
                out.append(intersect(prev, cur))
            out.append(cur)
        elif prev_in:
            out.append(intersect(prev, cur))
        prev, prev_in = cur, cur_in
    return out
