"""Float-only drawing helpers for Phase 5 event figures."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events.types import Event
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane
from vd3d.viz.convert import to_float
from vd3d.viz.plot_geometry import draw_line3d, draw_plane, draw_point3d
from vd3d.viz.random_geom import display_t_range, format_plane

PLANE_COLORS = ("steelblue", "orange", "seagreen", "mediumpurple", "tomato", "cadetblue")
LINE_COLORS = ("magenta", "cyan", "gold", "orangered", "slateblue", "limegreen")


def draw_planes(ax, planes: Sequence[Plane], *, lim: float) -> None:
    for i, plane in enumerate(planes):
        draw_plane(
            ax,
            plane,
            lim=lim,
            color=PLANE_COLORS[i % len(PLANE_COLORS)],
            label=format_plane(plane),
        )


def draw_intersection_lines(ax, lines: Sequence[Line3D], *, target: float = 3.0) -> None:
    for i, line in enumerate(lines):
        t_min, t_max = display_t_range(line.direction, target=target)
        draw_line3d(
            ax,
            line,
            t_min=t_min,
            t_max=t_max,
            color=LINE_COLORS[i % len(LINE_COLORS)],
            linewidth=2.8,
            label=f"L{line.id}  P{line.plane_a}∩P{line.plane_b}",
        )


def draw_triple_markers(ax, events: Sequence[Event]) -> None:
    for event in events:
        point = event.geometric_data
        label = f"({point.x}, {point.y}, {point.z})"
        draw_point3d(ax, point, color="red", size=130, label=label)
        ax.text(
            to_float(point.x),
            to_float(point.y),
            to_float(point.z) + 0.2,
            f"z={event.z}",
            color="red",
        )
