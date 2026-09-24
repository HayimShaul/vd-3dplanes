"""Draw a 2D vertical decomposition on a movable axis-aligned window."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from matplotlib.patches import Polygon, Rectangle

from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar
from vd3d.events.types import Event, EventType
from vd3d.geometry.plane import Plane
from vd3d.sweep.algorithm import SweepResult
from vd3d.sweep.locate import interval_containing
from vd3d.sweep.matching import cell_signature
from vd3d.sweep.update import event_vertex_plane_keys, vertex_plane_key
from vd3d.vertical_decomposition.types import VDCell2D, VerticalDecomposition
from vd3d.viz.convert import to_float
from vd3d.viz.plot_vd import _clip_above_line, _clip_below_line, _clip_halfplane


@dataclass
class ViewWindow2D:
    """Axis-aligned square window in the sweep plane."""

    cx: float = 0.0
    cy: float = 0.0
    half: float = 4.0

    def __post_init__(self) -> None:
        if self.half <= 0:
            raise ValueError("half must be positive")

    @property
    def xmin(self) -> float:
        return self.cx - self.half

    @property
    def xmax(self) -> float:
        return self.cx + self.half

    @property
    def ymin(self) -> float:
        return self.cy - self.half

    @property
    def ymax(self) -> float:
        return self.cy + self.half

    def contains(self, x: float, y: float, *, eps: float = 1e-9) -> bool:
        return (
            self.xmin - eps <= x <= self.xmax + eps
            and self.ymin - eps <= y <= self.ymax + eps
        )

    def corners(self) -> list[tuple[float, float]]:
        return [
            (self.xmin, self.ymin),
            (self.xmax, self.ymin),
            (self.xmax, self.ymax),
            (self.xmin, self.ymax),
        ]


# Eight visually distinct fills for adjacency colouring of sweep-plane cells.
_SWEEP_FACE_PALETTE: tuple[tuple[float, float, float], ...] = (
    (0.90, 0.30, 0.30),  # red
    (0.20, 0.55, 0.85),  # blue
    (0.30, 0.75, 0.35),  # green
    (0.95, 0.70, 0.15),  # amber
    (0.65, 0.35, 0.80),  # purple
    (0.15, 0.75, 0.75),  # teal
    (0.95, 0.45, 0.70),  # pink
    (0.55, 0.35, 0.15),  # brown
)


def adjacent_cell_colors(vd: VerticalDecomposition) -> dict[int, tuple[float, float, float]]:
    """Colour cells from an 8-colour palette; adjacent cells differ.

    Among colours not used by neighbours, prefer the globally least-used
    palette entry so the drawing spreads across all eight hues.
    """
    n_colors = len(_SWEEP_FACE_PALETTE)
    order = sorted(vd.cells, key=lambda cell: (-len(cell.neighbors), cell.id))
    color_index: dict[int, int] = {}
    usage = [0] * n_colors
    for cell in order:
        used = {color_index[nid] for nid in cell.neighbors if nid in color_index}
        free = [i for i in range(n_colors) if i not in used]
        if not free:
            # Degenerate: more than 8 mutually adjacent cells; wrap.
            free = list(range(n_colors))
        preferred = (cell.id * 3) % n_colors
        free.sort(key=lambda i: (usage[i], (i - preferred) % n_colors, i))
        chosen = free[0]
        color_index[cell.id] = chosen
        usage[chosen] += 1
    return {
        cell_id: _SWEEP_FACE_PALETTE[index] for cell_id, index in color_index.items()
    }


def cell_label(
    result: SweepResult | None,
    vd: VerticalDecomposition,
    cell: VDCell2D,
    z,
) -> int:
    """Prefer the active 3D cell id; fall back to the 2D cell id."""
    if result is None:
        return cell.id
    interval = interval_containing(result, z)
    if interval is None:
        return cell.id
    return interval.lookup().get(cell_signature(vd, cell), cell.id)


def cell_clip_polygon_window(
    vd: VerticalDecomposition, cell: VDCell2D, window: ViewWindow2D
) -> list[tuple[float, float]]:
    """``cell ∩ window`` as a convex polygon."""
    poly = window.corners()
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


def clip_line_to_window(line: Line2D, window: ViewWindow2D) -> list[tuple[float, float]]:
    """Segment of an infinite 2D line inside ``window``."""
    a, b, c = to_float(line.a), to_float(line.b), to_float(line.c)
    corners = window.corners()
    edges = (
        (corners[0], corners[1]),
        (corners[1], corners[2]),
        (corners[2], corners[3]),
        (corners[3], corners[0]),
    )

    def signed(p: tuple[float, float]) -> float:
        return a * p[0] + b * p[1] + c

    pts: list[tuple[float, float]] = []
    for p, q in edges:
        sp, sq = signed(p), signed(q)
        if abs(sp) <= 1e-12:
            pts.append(p)
        if abs(sq) <= 1e-12:
            pts.append(q)
        if sp * sq < -1e-24:
            t = sp / (sp - sq)
            pts.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    # Dedup
    uniq: list[tuple[float, float]] = []
    for x, y in pts:
        if all(abs(x - u) > 1e-9 or abs(y - v) > 1e-9 for u, v in uniq):
            uniq.append((x, y))
    if len(uniq) < 2:
        return []
    # Extreme pair along the line direction
    if abs(b) >= abs(a):
        uniq.sort(key=lambda p: p[0])
    else:
        uniq.sort(key=lambda p: p[1])
    return [uniq[0], uniq[-1]]


def event_highlights(
    vd: VerticalDecomposition,
    events: Sequence[Event],
    planes: Sequence[Plane],
) -> tuple[frozenset[int], frozenset[Scalar]]:
    """Plane ids (triple lines) and wall ``x`` values (alignment walls) to paint red."""
    plane_ids: set[int] = set()
    wall_xs: set[Scalar] = set()
    for event in events:
        if event.type is EventType.TRIPLE_INTERSECTION:
            plane_ids.update(event.plane_ids)
        elif event.type is EventType.VERTICAL_ALIGNMENT:
            keys = event_vertex_plane_keys(event, planes)
            for vertex in vd.arrangement.vertices:
                if vertex_plane_key(vd.arrangement, vertex) in keys:
                    wall_xs.add(vertex.point.x)
    return frozenset(plane_ids), frozenset(wall_xs)


def draw_sweep_slice(
    ax,
    vd: VerticalDecomposition,
    *,
    window: ViewWindow2D,
    z,
    title: str,
    result: SweepResult | None = None,
    event_points: Sequence[Point2D] = (),
    highlight_plane_ids: frozenset[int] | None = None,
    highlight_wall_xs: frozenset[Scalar] | None = None,
) -> None:
    """One sweep-plane panel: coloured cells, solid lines, dashed walls, event dots.

    Triple events: arrangement lines of the three planes are drawn solid red.
    Alignment events: the two Steiner walls at the aligning vertices are dashed red.
    """
    ax.clear()
    hot_planes = highlight_plane_ids or frozenset()
    hot_xs = highlight_wall_xs or frozenset()
    colors = adjacent_cell_colors(vd)
    for cell in vd.cells:
        poly = cell_clip_polygon_window(vd, cell, window)
        if len(poly) < 3:
            continue
        ax.add_patch(
            Polygon(
                poly,
                closed=True,
                facecolor=colors[cell.id],
                edgecolor="none",
                alpha=0.85,
                zorder=1,
            )
        )
        cx = sum(p[0] for p in poly) / len(poly)
        cy = sum(p[1] for p in poly) / len(poly)
        if window.contains(cx, cy):
            ax.text(
                cx,
                cy,
                str(cell_label(result, vd, cell, z)),
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
                zorder=4,
            )

    for line in vd.arrangement.lines:
        seg = clip_line_to_window(line, window)
        if len(seg) < 2:
            continue
        hot = line.source_plane_id is not None and line.source_plane_id in hot_planes
        ax.plot(
            [seg[0][0], seg[1][0]],
            [seg[0][1], seg[1][1]],
            color="red" if hot else "black",
            linestyle="-",
            linewidth=2.6 if hot else 1.8,
            zorder=4 if hot else 3,
        )

    for wall in vd.walls:
        x = to_float(wall.x)
        if x < window.xmin - 1e-9 or x > window.xmax + 1e-9:
            continue
        y0 = window.ymin if wall.y_min is None else max(window.ymin, to_float(wall.y_min))
        y1 = window.ymax if wall.y_max is None else min(window.ymax, to_float(wall.y_max))
        if y1 <= y0:
            continue
        hot = wall.x in hot_xs
        ax.plot(
            [x, x],
            [y0, y1],
            color="red" if hot else "black",
            linestyle="--",
            linewidth=2.4 if hot else 1.6,
            zorder=4 if hot else 3,
        )

    for point in event_points:
        x, y = to_float(point.x), to_float(point.y)
        if window.contains(x, y, eps=window.half * 0.02):
            ax.scatter(
                [x],
                [y],
                c="red",
                s=70,
                zorder=6,
                edgecolors="darkred",
                linewidths=0.6,
            )

    ax.add_patch(
        Rectangle(
            (window.xmin, window.ymin),
            2 * window.half,
            2 * window.half,
            fill=False,
            edgecolor="0.4",
            linestyle=":",
            linewidth=1.0,
            zorder=5,
        )
    )
    ax.set_xlim(window.xmin, window.xmax)
    ax.set_ylim(window.ymin, window.ymax)
    ax.set_aspect("equal")
    ax.grid(True, linestyle=":", alpha=0.35)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
