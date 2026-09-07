"""Phase 10 review scenes: zone degeneracies and simultaneous groups."""

from __future__ import annotations

import random

from matplotlib.figure import Figure

from vd3d.arrangement2d import build_line_arrangement
from vd3d.events import generate_all_events
from vd3d.events.samples import (
    planes_alignment_at_z2,
    planes_four_through_123,
    planes_vertical_and_slanted,
)
from vd3d.geometry import Line2D
from vd3d.geometry.points import Point3D
from vd3d.sweep import (
    compute_vd_around_event,
    group_events_by_z,
    locate_cell3d,
    vertical_decomposition_3d,
)
from vd3d.viz.plot_arrangement import draw_arrangement_lines, draw_vertices_numbered, setup_axes_2d
from vd3d.viz.plot_sweep import draw_matched_vd_pair, draw_query_points_3d
from vd3d.viz.plot_zone import (
    draw_crossings_numbered,
    draw_query_line,
    draw_zone_faces,
    draw_zone_walk_arrows,
    zone_viewing_limit,
)
from vd3d.viz.random_geom import choose_seed, random_general_position_planes, resolve_n
from vd3d.viz.scene import Scene
from vd3d.zone import compute_zone


def _triangle() -> list[Line2D]:
    return [
        Line2D(a=1, b=0, c=0, id=0),
        Line2D(a=0, b=1, c=0, id=1),
        Line2D(a=1, b=1, c=-1, id=2),
    ]


def _draw_zone(fig: Figure, lines: list[Line2D], query: Line2D, title: str) -> None:
    arr = build_line_arrangement(lines)
    zone = compute_zone(arr, query)
    extra = tuple(v.point for v in zone.vertices)
    lim = zone_viewing_limit(arr, query=query, crossings=zone.crossings, extra_points=extra)
    ax = fig.add_subplot(111)
    draw_zone_faces(ax, arr, zone, lim=lim)
    draw_arrangement_lines(ax, arr, lim=lim, color="0.25")
    draw_query_line(ax, query, lim=lim)
    draw_vertices_numbered(ax, arr)
    draw_crossings_numbered(ax, zone.crossings)
    draw_zone_walk_arrows(ax, zone)
    setup_axes_2d(
        ax,
        lim=lim,
        title=(
            f"{title}  —  {len(zone.crossings)} crossings, "
            f"{len(zone.vertices)} vertices, {len(zone.faces)} zone faces"
        ),
    )


def _draw_group_match(fig: Figure, planes, title: str) -> None:
    events = generate_all_events(planes)
    groups = group_events_by_z(events)
    group = next((item for item in groups if len(item.events) > 1), groups[0])
    event = group.representative
    vd_before, vd_after, z_minus, z_plus = compute_vd_around_event(
        planes, event, events
    )
    draw_matched_vd_pair(
        fig,
        vd_before,
        vd_after,
        title=f"{title}  {len(group.events)} events at z={group.z}",
        left_title=f"z− = {z_minus}",
        right_title=f"z+ = {z_plus}",
    )


def _draw_query_scene(fig: Figure, planes, title: str, *, seed: int = 10) -> None:
    result = vertical_decomposition_3d(planes)
    rng = random.Random(seed)
    points = [
        Point3D(rng.randint(-2, 2), rng.randint(-2, 2), rng.randint(-2, 2))
        for _ in range(24)
    ]
    cell_ids = [locate_cell3d(result, point) for point in points]
    ax = fig.add_subplot(111, projection="3d")
    n_cells = len(result.cells)
    n_groups = len(group_events_by_z(result.events))
    draw_query_points_3d(ax, result, points, cell_ids, lim=3.5)
    ax.set_title(
        f"{title}  —  {n_cells} 3D cells, {n_groups} event groups "
        f"(drag to rotate)"
    )


def _draw_through_vertex(fig: Figure) -> None:
    _draw_zone(
        fig,
        _triangle(),
        Line2D(a=1, b=-1, c=0),
        "Step 10.1: y=x through the origin",
    )


def _draw_overlap(fig: Figure) -> None:
    _draw_zone(
        fig,
        _triangle(),
        Line2D(a=0, b=1, c=0),
        "Step 10.2: query overlaps y=0",
    )


def _draw_alignment_group(fig: Figure) -> None:
    _draw_group_match(fig, planes_alignment_at_z2(), "Step 10.3: simultaneous group")


def _draw_four_planes(fig: Figure) -> None:
    _draw_group_match(fig, planes_four_through_123(), "Step 10.4: four triples at z=3")


def _draw_vertical(fig: Figure) -> None:
    _draw_query_scene(fig, planes_vertical_and_slanted(), "Step 10.4: vertical planes")


PHASE10_SCENES: tuple[Scene, ...] = (
    Scene(
        name="zone_through_vertex",
        title="Step 10.1: through a vertex",
        caption=(
            "What you must see: bold red query y=x through the origin. "
            "Crossing 0 sits on the numbered vertex at (0,0). The walk is "
            "unbounded → triangle → unbounded. The triangle is painted; "
            "wedges that only touch L at the origin stay gray."
        ),
        figsize=(7.5, 7.5),
        draw=_draw_through_vertex,
    ),
    Scene(
        name="zone_overlap",
        title="Step 10.2: overlapping query",
        caption=(
            "What you must see: bold red query on y=0, an arrangement line. "
            "Every face incident to y=0 is painted (including the triangle). "
            "Faces that do not touch y=0 are gray. Numbered crossings are "
            "the on-line vertices (0,0) and (1,0)."
        ),
        figsize=(7.5, 7.5),
        draw=_draw_overlap,
    ),
    Scene(
        name="simultaneous_group",
        title="Step 10.3: simultaneous events",
        caption=(
            "What you must see: one before/after pair for a whole z-group "
            "(several events at that height). Matched cells share a colour; "
            "unmatched cells are black only near the event strip."
        ),
        figsize=(11, 5.6),
        draw=_draw_alignment_group,
    ),
    Scene(
        name="four_planes_group",
        title="Step 10.4: four planes at a point",
        caption=(
            "What you must see: four triples at z=3 processed as one group. "
            "Left is z−, right is z+. Far cells keep their colour; the "
            "concurrent-line neighbourhood is the only recolor."
        ),
        figsize=(11, 5.6),
        draw=_draw_four_planes,
    ),
    Scene(
        name="vertical_planes",
        title="Step 10.4: vertical input planes",
        caption=(
            "What you must see: query points coloured by 3D cell for a "
            "scene with vertical planes x=0 and y=0. Rotate: a point on "
            "the wrong side of a plane is a matching bug. Boundary points "
            "are black."
        ),
        figsize=(8, 7),
        draw=_draw_vertical,
    ),
)


def _scene_group_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (3, 4, 5))
    planes = random_general_position_planes(rng, n)
    events = generate_all_events(planes)
    groups = group_events_by_z(events)
    n_groups = len(groups)
    n_multi = sum(1 for group in groups if len(group.events) > 1)

    def draw(fig: Figure) -> None:
        if not events:
            ax = fig.add_subplot(111)
            ax.set_title(f"Step 10: n={n}  —  no events")
            ax.set_axis_off()
            return
        _draw_group_match(fig, planes, f"Step 10: n={n}")

    if events:
        caption = (
            f"seed={seed} n={n}. {n_groups} event groups "
            f"({n_multi} simultaneous). Matched cells share a colour; "
            "unmatched (local) cells are black."
        )
    else:
        caption = f"seed={seed} n={n}. No events."

    return Scene(
        name="group_random",
        title="Step 10: event group",
        caption=caption,
        figsize=(11, 5.6),
        draw=draw,
    )


def make_phase10_scenes(seed: int, n: int | None = None) -> tuple[Scene, ...]:
    rng = random.Random(seed)
    return (_scene_group_random(rng, seed, n),)


def phase10_scenes(
    *, seed: int | None = None, fixtures: bool = False, n: int | None = None
) -> tuple[Scene, ...]:
    if fixtures:
        return PHASE10_SCENES
    return make_phase10_scenes(choose_seed(seed), n=n)
