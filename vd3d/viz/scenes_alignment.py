"""Phase 6 review scenes: y-parallel walls and alignment events."""

from __future__ import annotations

import random

from matplotlib.figure import Figure

from vd3d.arrangement2d import build_line_arrangement
from vd3d.events import (
    VerticalIntersectionLine,
    build_vertical_wall,
    compute_intersection_lines,
    compute_line_zone_on_wall,
    enumerate_alignment_pairs,
    generate_alignment_events,
    point_on_line_at_z,
)
from vd3d.events.samples import planes_alignment_at_z2
from vd3d.geometry import Plane, Point3D, build_slice_lines
from vd3d.viz.convert import to_float
from vd3d.viz.plot_alignment import (
    draw_alignment_pair,
    draw_wall_zone_scene,
    oracle_table_text,
)
from vd3d.viz.plot_arrangement import draw_arrangement_lines, draw_vertices_numbered, setup_axes_2d, viewing_limit
from vd3d.viz.plot_geometry import draw_line3d, draw_plane, set_equal_3d
from vd3d.viz.random_geom import (
    choose_seed,
    display_t_range,
    random_general_position_planes,
    resolve_n,
    view_lim,
)
from vd3d.viz.scene import Scene


def _l12(planes: list[Plane]):
    for line in compute_intersection_lines(planes):
        if {line.plane_a, line.plane_b} == {1, 2}:
            return line
    raise RuntimeError("missing L12")


def _draw_wall_3d(fig: Figure, planes: list[Plane], title: str) -> None:
    line = _l12(planes)
    wall = build_vertical_wall(line)
    t_min, t_max = display_t_range(line.direction, target=3.0)
    lim = view_lim(line.point_at(-1), line.point, line.point_at(1), Point3D(0, 3, 0), minimum=3.5)
    ax = fig.add_subplot(111, projection="3d")
    draw_plane(ax, wall, lim=lim, color="mediumpurple", alpha=0.35, label="wall  x = a z + b")
    draw_line3d(ax, line, t_min=t_min, t_max=t_max, color="magenta", linewidth=3.2, label="L")
    set_equal_3d(ax, lim=lim)
    ax.set_title(title)
    ax.legend(loc="upper left")


def _draw_slice(fig: Figure, planes: list[Plane], z: int, title: str, *, mark_align: bool) -> None:
    arr = build_line_arrangement(build_slice_lines(planes, z))
    events = enumerate_alignment_pairs(planes)
    lim = viewing_limit(arr, minimum=5)
    ax = fig.add_subplot(111)
    draw_arrangement_lines(ax, arr, lim=lim, color="0.3")
    draw_vertices_numbered(ax, arr)
    if mark_align and events:
        draw_alignment_pair(ax, events[0])
    else:
        lines = compute_intersection_lines(planes)
        pts = []
        for line in lines:
            if {line.plane_a, line.plane_b} in ({1, 2}, {3, 4}):
                try:
                    pts.append(point_on_line_at_z(line, z))
                except VerticalIntersectionLine:
                    continue
        for point in pts:
            ax.scatter(
                [to_float(point.x)],
                [to_float(point.y)],
                c="steelblue",
                s=90,
                zorder=7,
                edgecolors="black",
            )
        if len(pts) == 2 and pts[0].x != pts[1].x:
            ax.annotate(
                "not aligned",
                (
                    (to_float(pts[0].x) + to_float(pts[1].x)) / 2,
                    (to_float(pts[0].y) + to_float(pts[1].y)) / 2,
                ),
                fontsize=10,
                color="0.2",
            )
    setup_axes_2d(ax, lim=lim, title=title)


def _draw_wall_arr(fig: Figure, planes: list[Plane], title: str) -> None:
    line = _l12(planes)
    wall = build_vertical_wall(line)
    zone = compute_line_zone_on_wall(line, wall, planes)
    ax = fig.add_subplot(111)
    draw_wall_zone_scene(ax, zone, planes, title=title)


def _draw_zone_oracle(fig: Figure, planes: list[Plane], title: str) -> None:
    line = _l12(planes)
    wall = build_vertical_wall(line)
    zone = compute_line_zone_on_wall(line, wall, planes)
    oracle = enumerate_alignment_pairs(planes)
    ax = fig.add_subplot(111)
    draw_wall_zone_scene(ax, zone, planes, title=title)
    ax.text(
        0.02,
        0.02,
        oracle_table_text(oracle),
        transform=ax.transAxes,
        va="bottom",
        ha="left",
        fontsize=9,
        family="monospace",
        bbox={"facecolor": "white", "edgecolor": "0.6", "alpha": 0.9},
    )


PHASE6_SCENES: tuple[Scene, ...] = (
    Scene(
        name="wall_of_L",
        title="Step 6.1: y-parallel wall",
        caption=(
            "What you must see: a magenta intersection line L lying in a "
            "purple translucent plane. That plane does not change when you "
            "move along y — it is x = a z + b. Rotate: L stays in the mesh "
            "and the mesh is independent of y."
        ),
        figsize=(8, 7),
        draw=lambda fig: _draw_wall_3d(
            fig, planes_alignment_at_z2(), "Step 6.1: wall of L  (x = z)"
        ),
    ),
    Scene(
        name="slice_aligned_z2",
        title="Step 6.2: slice at z=2",
        caption=(
            "What you must see: the xy-slice at z=2. Two crimson vertices "
            "sit on one dashed vertical line (same x, different y). That is "
            "the alignment. No other vertex lies on the open segment between them."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_slice(
            fig, planes_alignment_at_z2(), 2, "Step 6.2: z=2  — aligned", mark_align=True
        ),
    ),
    Scene(
        name="slice_not_aligned_z3",
        title="Step 6.2: slice at z=3",
        caption=(
            "What you must see: the same two vertices at z=3. They no longer "
            "share an x — there is no dashed vertical through both. Alignment "
            "holds only at the predicted z=2."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_slice(
            fig, planes_alignment_at_z2(), 3, "Step 6.2: z=3  — not aligned", mark_align=False
        ),
    ),
    Scene(
        name="wall_arrangement",
        title="Step 6.3: wall arrangement",
        caption=(
            "What you must see: a 2D line arrangement in the (z, y) frame of "
            "the wall (z to the right, y up). The bold red line is L. A "
            "visible alignment is a vertical dashed segment (same z) from L "
            "to a green point. This is an ordinary line arrangement; L is "
            "one of the lines."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_wall_arr(
            fig, planes_alignment_at_z2(), "Step 6.3: wall arrangement, L bold red"
        ),
    ),
    Scene(
        name="zone_vs_oracle",
        title="Step 6.4: zone vs oracle",
        caption=(
            "What you must see: green points are visible alignments — each "
            "has a dashed vertical (same z) to L. Orange points are opposite "
            "vertices with a line crossing that vertical, so they are not "
            "events. Red points are on L (triples). Green z values match "
            "the oracle table."
        ),
        figsize=(8, 7.5),
        draw=lambda fig: _draw_zone_oracle(
            fig, planes_alignment_at_z2(), "Step 6.4: opposite (green) vs oracle"
        ),
    ),
)


def _scene_wall_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (4, 5))
    planes = random_general_position_planes(rng, n)
    lines = compute_intersection_lines(planes)
    line = lines[0]
    wall = build_vertical_wall(line)

    def draw(fig: Figure) -> None:
        t_min, t_max = display_t_range(line.direction)
        lim = view_lim(line.point_at(-1), line.point, line.point_at(1), minimum=3.5)
        ax = fig.add_subplot(111, projection="3d")
        draw_plane(ax, wall, lim=lim, color="mediumpurple", alpha=0.35, label="wall")
        draw_line3d(ax, line, t_min=t_min, t_max=t_max, color="magenta", linewidth=3.2, label="L")
        set_equal_3d(ax, lim=lim)
        ax.set_title(f"Step 6.1: n={n} wall of L")
        ax.legend(loc="upper left")

    return Scene(
        name="wall_random",
        title="Step 6.1: y-parallel wall",
        caption=(
            f"seed={seed} n={n}. Magenta L lies in the purple y-parallel wall. "
            "Rotate: L stays in the mesh; the wall does not depend on y."
        ),
        figsize=(8, 7),
        draw=draw,
    )


def _scene_zone_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (4, 5))
    planes = random_general_position_planes(rng, n)
    events = generate_alignment_events(planes)
    zs = ", ".join(str(e.z) for e in events) or "none"

    def draw(fig: Figure) -> None:
        line = compute_intersection_lines(planes)[0]
        wall = build_vertical_wall(line)
        zone = compute_line_zone_on_wall(line, wall, planes)
        ax = fig.add_subplot(111)
        draw_wall_zone_scene(ax, zone, planes, title=f"Step 6.4: n={n} wall zone")
        ax.text(
            0.02,
            0.02,
            oracle_table_text(enumerate_alignment_pairs(planes)),
            transform=ax.transAxes,
            va="bottom",
            ha="left",
            fontsize=8,
            family="monospace",
            bbox={"facecolor": "white", "edgecolor": "0.6", "alpha": 0.9},
        )

    return Scene(
        name="zone_oracle_random",
        title="Step 6.4: zone vs oracle",
        caption=(
            f"seed={seed} n={n}. Green = visible alignment (dashed vertical "
            f"to L). Orange = opposite but blocked. Red = on L. Green z "
            f"labels must match the oracle table. All alignments: z = {zs}."
        ),
        figsize=(8, 7.5),
        draw=draw,
    )


def make_phase6_scenes(seed: int, n: int | None = None) -> tuple[Scene, ...]:
    rng = random.Random(seed)
    return (
        _scene_wall_random(rng, seed, n),
        _scene_zone_random(rng, seed, n),
    )


def phase6_scenes(
    *, seed: int | None = None, fixtures: bool = False, n: int | None = None
) -> tuple[Scene, ...]:
    if fixtures:
        return PHASE6_SCENES
    return make_phase6_scenes(choose_seed(seed), n=n)
