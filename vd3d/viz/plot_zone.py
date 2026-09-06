"""Float-only drawing helpers for query-line zones."""

from __future__ import annotations

import colorsys
import math
from collections.abc import Callable, Sequence

import numpy as np
from matplotlib.patches import FancyArrowPatch, Polygon

from vd3d.arrangement2d.types import Arrangement2D, Face
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.viz.convert import to_float
from vd3d.viz.plot_arrangement import face_clip_polygon, viewing_limit
from vd3d.zone.types import Crossing, SupportingLineZone, Zone


def zone_viewing_limit(
    arrangement: Arrangement2D,
    *,
    query: Line2D | None = None,
    crossings: tuple[Crossing, ...] = (),
    extra_points: tuple[Point2D, ...] = (),
    margin: float = 1.5,
    minimum: float = 2.5,
) -> int:
    lim = float(viewing_limit(arrangement, margin=margin, minimum=minimum))
    for crossing in crossings:
        lim = max(lim, abs(to_float(crossing.point.x)) + margin, abs(to_float(crossing.point.y)) + margin)
    for point in extra_points:
        lim = max(lim, abs(to_float(point.x)) + margin, abs(to_float(point.y)) + margin)
    if query is not None:
        p0, p1 = query.sample_points()
        for point in (p0, p1):
            lim = max(lim, abs(to_float(point.x)) + 0.5, abs(to_float(point.y)) + 0.5)
    return max(2, math.ceil(lim))


def draw_query_line(ax, query: Line2D, *, lim: float, color: str = "red") -> None:
    a, b, c = to_float(query.a), to_float(query.b), to_float(query.c)
    if abs(b) >= abs(a):
        xs = np.array([-lim, lim], dtype=float)
        ys = -(a * xs + c) / b
    else:
        ys = np.array([-lim, lim], dtype=float)
        xs = -(b * ys + c) / a
    ax.plot(xs, ys, color=color, linewidth=4.2, zorder=4, solid_capstyle="round")


def draw_crossings_numbered(ax, crossings: tuple[Crossing, ...]) -> None:
    for i, crossing in enumerate(crossings):
        x, y = to_float(crossing.point.x), to_float(crossing.point.y)
        ax.scatter([x], [y], c="crimson", s=80, zorder=7, edgecolors="white", linewidths=1.2)
        ax.annotate(
            str(i),
            (x, y),
            textcoords="offset points",
            xytext=(7, 7),
            fontsize=11,
            fontweight="bold",
            color="crimson",
            zorder=8,
        )


def zone_face_color(index: int, count: int) -> tuple[float, float, float]:
    """Saturated distinct fill; unlike arrangement ``face_color`` this is only for the zone."""
    hue = (index + 0.5) / max(count, 1)
    return colorsys.hsv_to_rgb(hue % 1.0, 0.72, 0.96)


def paint_faces(
    ax,
    arrangement: Arrangement2D,
    faces: Sequence[Face],
    *,
    lim: float,
    label: Callable[[int, Face], str] | None = None,
) -> None:
    """Flood-fill ``faces`` (the zone) and leave every other face light gray."""
    ordered: list[Face] = []
    seen: set[int] = set()
    for face in faces:
        if face.id in seen:
            continue
        seen.add(face.id)
        ordered.append(face)

    for face in arrangement.faces:
        if face.id in seen:
            continue
        poly = face_clip_polygon(arrangement, face, lim)
        if len(poly) < 3:
            continue
        ax.add_patch(
            Polygon(
                poly,
                closed=True,
                facecolor="#e8e8e8",
                edgecolor="none",
                alpha=1.0,
                zorder=1,
            )
        )

    n = max(len(ordered), 1)
    for i, face in enumerate(ordered):
        poly = face_clip_polygon(arrangement, face, lim)
        if len(poly) < 3:
            continue
        ax.add_patch(
            Polygon(
                poly,
                closed=True,
                facecolor=zone_face_color(i, n),
                edgecolor="#222222",
                linewidth=0.9,
                alpha=0.95,
                zorder=2,
            )
        )
        if label is None:
            continue
        rx, ry = to_float(face.representative.x), to_float(face.representative.y)
        ax.text(
            rx,
            ry,
            label(i, face),
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            zorder=5,
        )


def draw_zone_faces(ax, arrangement: Arrangement2D, zone: Zone, *, lim: float) -> None:
    def _label(index: int, face: Face) -> str:
        return f"U{index}" if face.unbounded else str(index)

    paint_faces(ax, arrangement, zone.faces, lim=lim, label=_label)


def draw_supporting_faces(
    ax, arrangement: Arrangement2D, zone: SupportingLineZone, *, lim: float
) -> None:
    paint_faces(
        ax,
        arrangement,
        zone.faces,
        lim=lim,
        label=lambda _i, face: f"U{face.id}" if face.unbounded else str(face.id),
    )


def draw_zone_walk_arrows(ax, zone: Zone) -> None:
    for left, right in zip(zone.faces, zone.faces[1:]):
        x0, y0 = to_float(left.representative.x), to_float(left.representative.y)
        x1, y1 = to_float(right.representative.x), to_float(right.representative.y)
        ax.add_patch(
            FancyArrowPatch(
                (x0, y0),
                (x1, y1),
                arrowstyle="-|>",
                mutation_scale=14,
                color="black",
                linewidth=1.6,
                zorder=6,
            )
        )


def draw_supporting_vertices(ax, zone: SupportingLineZone) -> None:
    for vertex in zone.vertices_on_line:
        x, y = to_float(vertex.point.x), to_float(vertex.point.y)
        ax.scatter([x], [y], c="red", s=90, zorder=7, edgecolors="black", linewidths=0.5)
        ax.annotate(
            "on L",
            (x, y),
            textcoords="offset points",
            xytext=(8, 8),
            fontsize=9,
            color="red",
            fontweight="bold",
            zorder=8,
        )
    for vertex in zone.opposite_vertices:
        x, y = to_float(vertex.point.x), to_float(vertex.point.y)
        ax.scatter([x], [y], c="limegreen", s=90, zorder=7, edgecolors="black", linewidths=0.5)
        ax.annotate(
            "opposite",
            (x, y),
            textcoords="offset points",
            xytext=(8, -14),
            fontsize=9,
            color="green",
            fontweight="bold",
            zorder=8,
        )


__all__ = [
    "draw_crossings_numbered",
    "draw_query_line",
    "draw_supporting_faces",
    "draw_supporting_vertices",
    "draw_zone_faces",
    "draw_zone_walk_arrows",
    "paint_faces",
    "zone_face_color",
    "zone_viewing_limit",
]
