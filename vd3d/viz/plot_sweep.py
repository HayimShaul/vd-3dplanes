"""Float-only drawing helpers for Phase 8 matching, snapshots, and 3D cells."""

from __future__ import annotations

from collections.abc import Sequence

from matplotlib.patches import Polygon

from vd3d.geometry.points import Point3D
from vd3d.sweep.matching import CellMatch, match_cells
from vd3d.sweep.types import SweepSnapshot
from vd3d.vertical_decomposition.types import VerticalDecomposition
from vd3d.viz.convert import to_float
from vd3d.viz.plot_arrangement import draw_arrangement_lines, draw_vertices_numbered, face_color, setup_axes_2d
from vd3d.viz.plot_event_list import shared_vd_limit
from vd3d.viz.plot_geometry import draw_plane, draw_point3d, set_equal_3d
from vd3d.viz.plot_vd import cell_clip_polygon, draw_vertical_walls


UNMATCHED_COLOR = (0.12, 0.12, 0.12)


def match_color_maps(
    vd_before: VerticalDecomposition,
    vd_after: VerticalDecomposition,
    matches: Sequence[CellMatch],
) -> tuple[dict[int, tuple[float, float, float]], dict[int, tuple[float, float, float]]]:
    """Matched cells share a hue; unmatched cells are black."""
    n = max(len(matches), 1)
    before: dict[int, tuple[float, float, float]] = {
        cell.id: UNMATCHED_COLOR for cell in vd_before.cells
    }
    after: dict[int, tuple[float, float, float]] = {
        cell.id: UNMATCHED_COLOR for cell in vd_after.cells
    }
    for index, pair in enumerate(matches):
        color = face_color(index, n)
        before[pair.before.id] = color
        after[pair.after.id] = color
    return before, after


def draw_cells_colored(
    ax,
    vd: VerticalDecomposition,
    colors: dict[int, tuple[float, float, float]],
    *,
    lim: float,
) -> None:
    for cell in vd.cells:
        poly = cell_clip_polygon(vd, cell, lim)
        if len(poly) < 3:
            continue
        ax.add_patch(
            Polygon(
                poly,
                closed=True,
                facecolor=colors.get(cell.id, UNMATCHED_COLOR),
                edgecolor="none",
                alpha=0.9,
                zorder=1,
            )
        )
        cx = sum(p[0] for p in poly) / len(poly)
        cy = sum(p[1] for p in poly) / len(poly)
        ax.text(cx, cy, str(cell.id), ha="center", va="center", fontsize=9, zorder=4)


def draw_matched_vd_pair(
    fig,
    vd_before: VerticalDecomposition,
    vd_after: VerticalDecomposition,
    *,
    title: str,
    left_title: str,
    right_title: str,
) -> tuple[CellMatch, ...]:
    matches = match_cells(vd_before, vd_after)
    colors_b, colors_a = match_color_maps(vd_before, vd_after, matches)
    lim = shared_vd_limit(vd_before, vd_after)
    ax_l = fig.add_subplot(121)
    ax_r = fig.add_subplot(122)
    _draw_vd_colored(ax_l, vd_before, colors_b, lim=lim, title=left_title)
    _draw_vd_colored(ax_r, vd_after, colors_a, lim=lim, title=right_title)
    n_un_b = len(vd_before.cells) - len(matches)
    n_un_a = len(vd_after.cells) - len(matches)
    fig.suptitle(
        f"{title}  —  {len(matches)} matched, {n_un_b} ended (black), {n_un_a} started (black)"
    )
    return matches


def draw_snapshot_caption(snapshot: SweepSnapshot) -> str:
    return (
        f"{snapshot.event.type.name} z={snapshot.event.z}  "
        f"before={snapshot.n_before} after={snapshot.n_after}  "
        f"continued={len(snapshot.matches)} ended={len(snapshot.ended)} "
        f"started={len(snapshot.started)} active={snapshot.n_active}"
    )


def draw_query_points_3d(
    ax,
    result,
    points: Sequence[Point3D],
    cell_ids: Sequence[int | None],
    *,
    lim: float,
) -> None:
    plane_colors = ("steelblue", "orange", "seagreen", "mediumpurple", "goldenrod")
    for i, plane in enumerate(result.planes):
        draw_plane(
            ax,
            plane,
            lim=lim,
            color=plane_colors[i % len(plane_colors)],
            alpha=0.18,
            label=f"P{plane.id}",
        )
    n = max(len(result.cells), 1)
    for point, cell_id in zip(points, cell_ids):
        color = UNMATCHED_COLOR if cell_id is None else face_color(cell_id, n)
        draw_point3d(ax, point, color=color, size=50)
        if cell_id is not None:
            ax.text(
                to_float(point.x),
                to_float(point.y),
                to_float(point.z),
                str(cell_id),
                fontsize=7,
            )
    set_equal_3d(ax, lim=lim)
    ax.legend(loc="upper left", fontsize=8)


def _draw_vd_colored(ax, vd, colors, *, lim: float, title: str) -> None:
    draw_cells_colored(ax, vd, colors, lim=lim)
    draw_arrangement_lines(ax, vd.arrangement, lim=lim, color="0.15")
    draw_vertical_walls(ax, vd, lim=lim)
    draw_vertices_numbered(ax, vd.arrangement)
    setup_axes_2d(
        ax,
        lim=lim,
        title=f"{title}  —  {len(vd.arrangement.vertices)} vertices, {len(vd.cells)} cells",
    )
