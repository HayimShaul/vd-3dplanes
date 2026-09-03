"""Float-only drawing helpers for 2D line arrangements."""

from __future__ import annotations

import colorsys
import math

from matplotlib.patches import FancyArrowPatch, Polygon

from vd3d.arrangement2d.invariants import point_in_face
from vd3d.arrangement2d.predicates import line_direction
from vd3d.arrangement2d.types import Arrangement2D, EdgePiece, Face, HalfEdge
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.viz.convert import to_float
from vd3d.viz.plot_geometry import draw_line2d

_PIECE_COLORS = [
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    "#a65628",
    "#f781bf",
    "#999999",
    "#66c2a5",
    "#fc8d62",
    "#8da0cb",
    "#e78ac3",
]


def format_line2d(line: Line2D) -> str:
    return f"{line.a}x+{line.b}y+{line.c}=0"


def face_color(index: int, count: int) -> tuple[float, float, float]:
    hue = index / max(count, 1)
    return colorsys.hsv_to_rgb(hue, 0.45, 0.93)


def viewing_limit(arrangement: Arrangement2D, *, margin: float = 1.5, minimum: float = 2.5) -> int:
    """An integer half-width so window corners are exact ``Point2D`` values."""
    lim = minimum
    for vertex in arrangement.vertices:
        lim = max(lim, abs(to_float(vertex.point.x)) + margin, abs(to_float(vertex.point.y)) + margin)
    if not arrangement.vertices:
        lim = max(lim, 3.0)
    return max(2, math.ceil(lim))


def setup_axes_2d(ax, *, lim: float, title: str) -> None:
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    ax.axhline(0, color="0.7", linewidth=0.6)
    ax.axvline(0, color="0.7", linewidth=0.6)


def draw_arrangement_lines(ax, arrangement: Arrangement2D, *, lim: float, color: str = "0.2") -> None:
    for line in arrangement.lines:
        draw_line2d(ax, line, lim=lim, color=color)


def draw_vertices_numbered(ax, arrangement: Arrangement2D) -> None:
    for vertex in arrangement.vertices:
        x, y = to_float(vertex.point.x), to_float(vertex.point.y)
        ax.scatter([x], [y], c="black", s=40, zorder=5)
        ax.annotate(
            str(vertex.id),
            (x, y),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=9,
            fontweight="bold",
            zorder=6,
        )


def draw_edge_pieces(ax, arrangement: Arrangement2D, *, lim: float, arrows: bool = True) -> None:
    for edge in arrangement.edges:
        color = _PIECE_COLORS[edge.id % len(_PIECE_COLORS)]
        pts = clip_edge_to_window(arrangement, edge, lim)
        if len(pts) >= 2:
            xs, ys = zip(*pts)
            ax.plot(xs, ys, color=color, linewidth=2.4, solid_capstyle="round", zorder=3)
        if arrows and edge.kind == "ray":
            _draw_ray_arrow(ax, arrangement, edge, lim, color)


def clip_edge_to_window(
    arrangement: Arrangement2D,
    edge: EdgePiece,
    lim: float,
) -> list[tuple[float, float]]:
    dx, dy = (to_float(c) for c in line_direction(arrangement.lines[edge.line_index]))
    mag = math.hypot(dx, dy)
    if mag == 0:
        return []
    dx, dy = dx / mag, dy / mag

    if edge.kind == "segment":
        if edge.start_vertex is None or edge.end_vertex is None:
            return []
        a = arrangement.vertices[edge.start_vertex].point
        b = arrangement.vertices[edge.end_vertex].point
        return _clip_segment(to_float(a.x), to_float(a.y), to_float(b.x), to_float(b.y), lim)

    if edge.kind == "ray":
        finite = edge.end_vertex if edge.start_vertex is None else edge.start_vertex
        if finite is None:
            return []
        p = arrangement.vertices[finite].point
        px, py = to_float(p.x), to_float(p.y)
        if edge.start_vertex is None:
            rdx, rdy = -dx, -dy
        else:
            rdx, rdy = dx, dy
        return _clip_ray(px, py, rdx, rdy, lim)

    sample, _ = arrangement.lines[edge.line_index].sample_points()
    return _clip_infinite_line(to_float(sample.x), to_float(sample.y), dx, dy, lim)


def _draw_ray_arrow(
    ax,
    arrangement: Arrangement2D,
    edge: EdgePiece,
    lim: float,
    color: str,
) -> None:
    finite = edge.end_vertex if edge.start_vertex is None else edge.start_vertex
    if finite is None:
        return
    p = arrangement.vertices[finite].point
    px, py = to_float(p.x), to_float(p.y)
    dx, dy = (to_float(c) for c in line_direction(arrangement.lines[edge.line_index]))
    mag = math.hypot(dx, dy)
    if mag == 0:
        return
    ux, uy = dx / mag, dy / mag
    if edge.start_vertex is None:
        ux, uy = -ux, -uy
    span = 0.55 * min(1.2, lim / 3)
    start = (px + 0.15 * ux, py + 0.15 * uy)
    end = (px + span * ux, py + span * uy)
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            color=color,
            linewidth=1.6,
            zorder=4,
        )
    )


def draw_outgoing_half_edges(
    ax,
    arrangement: Arrangement2D,
    vertex_id: int = 0,
    *,
    scale: float = 1.1,
) -> None:
    vertex = arrangement.vertices[vertex_id]
    ox, oy = to_float(vertex.point.x), to_float(vertex.point.y)
    ax.scatter([ox], [oy], c="black", s=70, zorder=5)
    for i, hid in enumerate(vertex.outgoing):
        he = arrangement.half_edges[hid]
        dx, dy = to_float(he.direction[0]), to_float(he.direction[1])
        mag = math.hypot(dx, dy)
        ux, uy = dx / mag, dy / mag
        end = (ox + scale * ux, oy + scale * uy)
        ax.add_patch(
            FancyArrowPatch(
                (ox, oy),
                end,
                arrowstyle="-|>",
                mutation_scale=14,
                color="crimson",
                linewidth=2,
                zorder=4,
            )
        )
        ax.text(
            ox + (scale + 0.18) * ux,
            oy + (scale + 0.18) * uy,
            str(i),
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
            color="crimson",
            zorder=6,
        )


def draw_faces_filled(ax, arrangement: Arrangement2D, *, lim: float) -> None:
    n = len(arrangement.faces)
    for face in arrangement.faces:
        poly = face_clip_polygon(arrangement, face, lim)
        if len(poly) < 3:
            continue
        ax.add_patch(
            Polygon(
                poly,
                closed=True,
                facecolor=face_color(face.id, n),
                edgecolor="none",
                alpha=0.9,
                zorder=1,
            )
        )
        cx = sum(p[0] for p in poly) / len(poly)
        cy = sum(p[1] for p in poly) / len(poly)
        label = f"U{face.id}" if face.unbounded else str(face.id)
        ax.text(cx, cy, label, ha="center", va="center", fontsize=9, zorder=4)


def face_clip_polygon(
    arrangement: Arrangement2D,
    face: Face,
    lim: float,
) -> list[tuple[float, float]]:
    """Vertices of ``face ∩ [-lim, lim]^2``, sorted CCW. Faces are convex."""
    pts: list[tuple[float, float]] = []
    for he in arrangement.cycle(face):
        if he.origin_id is not None:
            p = arrangement.vertices[he.origin_id].point
            x, y = to_float(p.x), to_float(p.y)
            if _inside_box(x, y, lim):
                pts.append((x, y))
        pts.extend(_half_edge_box_hits(arrangement, he, lim))

    ilim = int(lim)
    for x, y in ((-ilim, -ilim), (ilim, -ilim), (ilim, ilim), (-ilim, ilim)):
        if point_in_face(arrangement, face, Point2D(x, y), closed=True):
            pts.append((float(x), float(y)))

    pts = _unique_points(pts)
    return _sort_ccw(pts)


def _half_edge_box_hits(
    arrangement: Arrangement2D,
    he: HalfEdge,
    lim: float,
) -> list[tuple[float, float]]:
    edge = arrangement.edges[he.edge_id]
    pts = clip_edge_to_window(arrangement, edge, lim)
    return [p for p in pts if _on_box_boundary(p[0], p[1], lim)]


def _on_box_boundary(x: float, y: float, lim: float, eps: float = 1e-9) -> bool:
    on_x = abs(abs(x) - lim) <= eps and abs(y) <= lim + eps
    on_y = abs(abs(y) - lim) <= eps and abs(x) <= lim + eps
    return on_x or on_y


def _inside_box(x: float, y: float, lim: float, eps: float = 1e-12) -> bool:
    return -lim - eps <= x <= lim + eps and -lim - eps <= y <= lim + eps


def _unique_points(pts: list[tuple[float, float]], ndigits: int = 9) -> list[tuple[float, float]]:
    seen: dict[tuple[float, float], tuple[float, float]] = {}
    for x, y in pts:
        seen[(round(x, ndigits), round(y, ndigits))] = (x, y)
    return list(seen.values())


def _sort_ccw(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
    if len(pts) < 3:
        return pts
    cx = sum(p[0] for p in pts) / len(pts)
    cy = sum(p[1] for p in pts) / len(pts)
    return sorted(pts, key=lambda p: math.atan2(p[1] - cy, p[0] - cx))


def _clip_segment(
    x0: float, y0: float, x1: float, y1: float, lim: float
) -> list[tuple[float, float]]:
    hits = [(x0, y0, 0.0), (x1, y1, 1.0)]
    dx, dy = x1 - x0, y1 - y0
    hits.extend(_axis_hits(x0, y0, dx, dy, lim))
    inside = [(x, y) for x, y, t in hits if 0 <= t <= 1 and _inside_box(x, y, lim)]
    inside = _unique_points(inside)
    inside.sort(key=lambda p: (p[0] - x0) * dx + (p[1] - y0) * dy)
    return inside


def _clip_ray(px: float, py: float, dx: float, dy: float, lim: float) -> list[tuple[float, float]]:
    hits = _axis_hits(px, py, dx, dy, lim)
    if _inside_box(px, py, lim):
        hits.append((px, py, 0.0))
    positive = [(x, y, t) for x, y, t in hits if t >= -1e-12 and _inside_box(x, y, lim)]
    positive.sort(key=lambda item: item[2])
    pts = _unique_points([(x, y) for x, y, _ in positive])
    if len(pts) >= 2:
        return [pts[0], pts[-1]]
    return pts


def _clip_infinite_line(
    px: float, py: float, dx: float, dy: float, lim: float
) -> list[tuple[float, float]]:
    hits = _axis_hits(px, py, dx, dy, lim)
    pts = _unique_points([(x, y) for x, y, _ in hits if _inside_box(x, y, lim)])
    if len(pts) >= 2:
        return [pts[0], pts[-1]]
    return pts


def _axis_hits(
    px: float, py: float, dx: float, dy: float, lim: float
) -> list[tuple[float, float, float]]:
    hits: list[tuple[float, float, float]] = []
    if dx != 0:
        for xedge in (-lim, lim):
            t = (xedge - px) / dx
            y = py + t * dy
            if abs(y) <= lim + 1e-9:
                hits.append((xedge, y, t))
    if dy != 0:
        for yedge in (-lim, lim):
            t = (yedge - py) / dy
            x = px + t * dx
            if abs(x) <= lim + 1e-9:
                hits.append((x, yedge, t))
    return hits
