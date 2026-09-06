"""Float-only drawing helpers for Phase 7 event-list and slice-VD figures."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events.types import Event, EventType
from vd3d.vertical_decomposition import VerticalDecomposition
from vd3d.viz.convert import to_float
from vd3d.viz.plot_arrangement import draw_arrangement_lines, draw_vertices_numbered, setup_axes_2d
from vd3d.viz.plot_vd import draw_cells_filled, draw_vertical_walls, vd_viewing_limit

TYPE_COLOR = {
    EventType.TRIPLE_INTERSECTION: "crimson",
    EventType.VERTICAL_ALIGNMENT: "seagreen",
}
TYPE_LABEL = {
    EventType.TRIPLE_INTERSECTION: "triple",
    EventType.VERTICAL_ALIGNMENT: "align",
}


def draw_event_timeline(ax, events: Sequence[Event], *, title: str) -> None:
    """1D z-axis with type-colored ticks."""
    zs = [to_float(event.z) for event in events]
    if zs:
        lo, hi = min(zs) - 1.5, max(zs) + 1.5
    else:
        lo, hi = -2.0, 2.0
    ax.plot([lo, hi], [0.0, 0.0], color="0.35", linewidth=2.0, zorder=1)
    seen: dict[float, int] = {}
    for event in events:
        z = to_float(event.z)
        stack = seen.get(z, 0)
        seen[z] = stack + 1
        y = 0.18 + 0.28 * stack
        color = TYPE_COLOR[event.type]
        ax.plot([z, z], [0.0, y], color=color, linewidth=2.2, zorder=2)
        ax.scatter([z], [y], c=color, s=90, zorder=3, edgecolors="black", linewidths=0.4)
        ax.annotate(
            f"{TYPE_LABEL[event.type]}  z={event.z}",
            (z, y),
            textcoords="offset points",
            xytext=(8, 4),
            fontsize=8,
            color=color,
        )
    ax.set_xlim(lo, hi)
    ax.set_ylim(-0.55, 1.55)
    ax.set_yticks([])
    ax.set_xlabel("z")
    ax.set_title(title)
    ax.axvline(0, color="0.75", linewidth=0.6, linestyle=":")
    handles = [
        ax.scatter([], [], c=TYPE_COLOR[kind], s=60, edgecolors="black", label=TYPE_LABEL[kind])
        for kind in (EventType.TRIPLE_INTERSECTION, EventType.VERTICAL_ALIGNMENT)
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.set_aspect("auto")


def draw_vd_slice(ax, vd: VerticalDecomposition, *, lim: float, title: str) -> None:
    draw_cells_filled(ax, vd, lim=lim)
    draw_arrangement_lines(ax, vd.arrangement, lim=lim, color="0.15")
    draw_vertical_walls(ax, vd, lim=lim)
    draw_vertices_numbered(ax, vd.arrangement)
    setup_axes_2d(
        ax,
        lim=lim,
        title=f"{title}  —  {len(vd.arrangement.vertices)} vertices, {len(vd.cells)} cells",
    )


def shared_vd_limit(*vds: VerticalDecomposition, minimum: float = 3.5) -> float:
    lim = minimum
    for vd in vds:
        lim = max(lim, float(vd_viewing_limit(vd, minimum=minimum)))
    return lim
