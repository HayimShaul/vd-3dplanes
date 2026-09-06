"""Float-only drawing helpers for Phase 6 alignment figures."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events import (
    AlignmentGeometry,
    Event,
    alignment_event_for_opposite_vertex,
    compute_intersection_lines,
    lift_wall_point,
    wall_frame_point,
)
from vd3d.events.alignment import WallLineZone
from vd3d.geometry.plane import Plane
from vd3d.viz.convert import to_float
from vd3d.viz.plot_arrangement import draw_arrangement_lines
from vd3d.viz.plot_zone import draw_query_line, draw_supporting_faces, zone_viewing_limit


def setup_wall_axes(ax, *, lim: float, title: str) -> None:
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.set_xlabel("z  (sweep)")
    ax.set_ylabel("y  (2D vertical)")
    ax.set_title(title)
    ax.axhline(0, color="0.7", linewidth=0.6)
    ax.axvline(0, color="0.7", linewidth=0.6)


def draw_dashed_vertical(ax, x, *, y_min, y_max, color: str = "black") -> None:
    xf = to_float(x)
    ax.plot(
        [xf, xf],
        [to_float(y_min), to_float(y_max)],
        color=color,
        linestyle="--",
        linewidth=1.8,
        zorder=4,
    )


def draw_alignment_pair(ax, event: Event) -> None:
    data = event.geometric_data
    if not isinstance(data, AlignmentGeometry):
        return
    for point, color in ((data.point_a, "crimson"), (data.point_b, "crimson")):
        ax.scatter(
            [to_float(point.x)],
            [to_float(point.y)],
            c=color,
            s=110,
            zorder=8,
            edgecolors="black",
        )
    draw_dashed_vertical(
        ax,
        data.point_a.x,
        y_min=min(data.point_a.y, data.point_b.y) - 1,
        y_max=max(data.point_a.y, data.point_b.y) + 1,
    )


def draw_wall_zone_vertices(ax, zone: WallLineZone, planes: Sequence[Plane]) -> None:
    """Red = on L (triples). Green = visible alignment. Orange = opposite but blocked."""
    lines = compute_intersection_lines(planes)
    for vertex in zone.supporting.vertices_on_line:
        x, y = to_float(vertex.point.x), to_float(vertex.point.y)
        ax.scatter([x], [y], c="red", s=90, zorder=7, edgecolors="black")
        lifted = lift_wall_point(vertex.point, zone.wall)
        ax.annotate(
            f"triple z={lifted.z}",
            (x, y),
            textcoords="offset points",
            xytext=(8, 8),
            fontsize=9,
            color="red",
            fontweight="bold",
        )
    for vertex in zone.supporting.opposite_vertices:
        x, y = to_float(vertex.point.x), to_float(vertex.point.y)
        event = alignment_event_for_opposite_vertex(zone, vertex, planes, lines)
        if event is None:
            ax.scatter([x], [y], c="darkorange", s=70, zorder=7, edgecolors="black")
            ax.annotate(
                "blocked",
                (x, y),
                textcoords="offset points",
                xytext=(8, -14),
                fontsize=9,
                color="darkorange",
                fontweight="bold",
            )
            continue
        ax.scatter([x], [y], c="limegreen", s=90, zorder=7, edgecolors="black")
        on_L = wall_frame_point(zone.line, event.z)
        draw_dashed_vertical(ax, vertex.point.x, y_min=on_L.y, y_max=vertex.point.y, color="limegreen")
        ax.annotate(
            f"z={event.z}",
            (x, y),
            textcoords="offset points",
            xytext=(8, -14),
            fontsize=9,
            color="green",
            fontweight="bold",
        )


def oracle_table_text(events: Sequence[Event]) -> str:
    if not events:
        return "oracle: (none)"
    rows = []
    for event in events:
        rows.append(f"z={event.z}  lines {event.line_ids}")
    return "oracle:\n" + "\n".join(rows)


def draw_wall_zone_scene(ax, zone: WallLineZone, planes: Sequence[Plane], *, title: str) -> None:
    extra = tuple(
        v.point
        for v in (*zone.supporting.vertices_on_line, *zone.supporting.opposite_vertices)
    )
    lim = zone_viewing_limit(zone.arrangement, extra_points=extra)
    draw_supporting_faces(ax, zone.arrangement, zone.supporting, lim=lim)
    draw_arrangement_lines(ax, zone.arrangement, lim=lim, color="0.25")
    draw_query_line(ax, zone.arrangement.lines[zone.query_index], lim=lim)
    draw_wall_zone_vertices(ax, zone, planes)
    setup_wall_axes(ax, lim=lim, title=title)
