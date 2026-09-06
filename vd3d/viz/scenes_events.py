"""Phase 5 review scenes: intersection lines and triple events."""

from __future__ import annotations

import random

from matplotlib.figure import Figure

from vd3d.events import compute_intersection_lines, generate_triple_events
from vd3d.geometry import Plane, Point3D
from vd3d.viz.plot_events import draw_intersection_lines, draw_planes, draw_triple_markers
from vd3d.viz.plot_geometry import draw_height_guide, draw_z_axis, set_equal_3d
from vd3d.viz.random_geom import (
    choose_seed,
    random_general_position_planes,
    random_parallel_family,
    resolve_n,
    view_lim,
)
from vd3d.viz.scene import Scene


def planes_through_123() -> list[Plane]:
    """Design Test 15: three planes meeting at ``(1, 2, 3)``."""
    return [
        Plane(id=1, a=1, b=0, c=1, d=-4),
        Plane(id=2, a=0, b=1, c=1, d=-5),
        Plane(id=3, a=1, b=1, c=1, d=-6),
    ]


def parallel_family_fixture() -> list[Plane]:
    """Design Test 16: pairwise-parallel planes."""
    return [
        Plane(id=1, a=1, b=1, c=1, d=0),
        Plane(id=2, a=1, b=1, c=1, d=-1),
        Plane(id=3, a=1, b=1, c=1, d=-2),
    ]


def _view_lim_for(planes: list[Plane], extra: tuple[Point3D, ...] = ()) -> float:
    lines = compute_intersection_lines(planes)
    points = list(extra)
    for line in lines:
        points.append(line.point)
        points.append(line.point_at(1))
        points.append(line.point_at(-1))
    if not points:
        points.append(Point3D(0, 0, 0))
        if planes:
            # A point on the first plane so parallel families still frame.
            p0 = planes[0]
            if p0.c != 0:
                points.append(Point3D(0, 0, -p0.d / p0.c))
            elif p0.b != 0:
                points.append(Point3D(0, -p0.d / p0.b, 0))
            else:
                points.append(Point3D(-p0.d / p0.a, 0, 0))
    return view_lim(*points, minimum=3.5)


def _draw_lines_scene(fig: Figure, planes: list[Plane], title: str) -> None:
    lines = compute_intersection_lines(planes)
    lim = _view_lim_for(planes)
    ax = fig.add_subplot(111, projection="3d")
    draw_planes(ax, planes, lim=lim)
    draw_intersection_lines(ax, lines, target=lim)
    set_equal_3d(ax, lim=lim)
    ax.set_title(f"{title}  —  {len(lines)} intersection lines")
    ax.legend(loc="upper left", fontsize=8)


def _draw_triples_scene(fig: Figure, planes: list[Plane], title: str) -> None:
    lines = compute_intersection_lines(planes)
    events = generate_triple_events(planes)
    extra = tuple(event.geometric_data for event in events)
    lim = _view_lim_for(planes, extra)
    ax = fig.add_subplot(111, projection="3d")
    draw_planes(ax, planes, lim=lim)
    draw_intersection_lines(ax, lines, target=lim)
    if events:
        for event in events:
            draw_z_axis(
                ax,
                lim=lim,
                tick_z=float(event.z),
                tick_label=f"z={event.z}",
            )
            draw_height_guide(ax, event.geometric_data)
    else:
        draw_z_axis(ax, lim=lim)
    draw_triple_markers(ax, events)
    set_equal_3d(ax, lim=lim)
    ax.set_title(f"{title}  —  {len(events)} triple event(s)")
    ax.legend(loc="upper left", fontsize=8)


def _draw_parallel_scene(fig: Figure, planes: list[Plane], title: str) -> None:
    lines = compute_intersection_lines(planes)
    events = generate_triple_events(planes)
    lim = _view_lim_for(planes)
    ax = fig.add_subplot(111, projection="3d")
    draw_planes(ax, planes, lim=lim)
    draw_z_axis(ax, lim=lim)
    set_equal_3d(ax, lim=lim)
    ax.set_title(f"{title}  —  {len(lines)} lines, {len(events)} triples")
    ax.legend(loc="upper left", fontsize=8)


PHASE5_SCENES: tuple[Scene, ...] = (
    Scene(
        name="intersection_lines_three_planes",
        title="Step 5.1: intersection lines",
        caption=(
            "What you must see: three translucent planes and three colored "
            "intersection lines (one crease per pair). The three lines meet "
            "at a single point. Rotate: each colored line must stay the "
            "crease of its two planes; there is no fourth line."
        ),
        figsize=(8, 7),
        draw=lambda fig: _draw_lines_scene(
            fig, planes_through_123(), "Step 5.1: 3 planes, 3 lines"
        ),
    ),
    Scene(
        name="triple_event_z3",
        title="Step 5.2: triple event",
        caption=(
            "What you must see: one fat red marker at (1, 2, 3). A tick on "
            "the z-axis at height 3, labeled z=3, and a dashed guide at that "
            "same height from the axis to the marker. One obvious point, one "
            "obvious z. Rotate until the marker sits on all three meshes."
        ),
        figsize=(8, 7),
        draw=lambda fig: _draw_triples_scene(
            fig, planes_through_123(), "Step 5.2: triple at (1,2,3), z=3"
        ),
    ),
    Scene(
        name="parallel_family_no_triples",
        title="Step 5.2: parallel family",
        caption=(
            "What you must see: three parallel planes (same normal, shifted). "
            "No intersection lines. No red triple marker. No z-tick other "
            "than the bare axis. Rotate to confirm they never meet."
        ),
        figsize=(8, 7),
        draw=lambda fig: _draw_parallel_scene(
            fig, parallel_family_fixture(), "Step 5.2: parallel family"
        ),
    ),
)


def _scene_lines_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (3, 4))
    planes = random_general_position_planes(rng, n)
    lines = compute_intersection_lines(planes)

    def draw(fig: Figure) -> None:
        _draw_lines_scene(fig, planes, f"Step 5.1: n={n} intersection lines")

    return Scene(
        name="intersection_lines_random",
        title="Step 5.1: intersection lines",
        caption=(
            f"seed={seed} n={n}. {len(lines)} colored intersection lines "
            f"(expected {n * (n - 1) // 2}). Each line is the crease of two "
            "planes. Rotate: no line should float off its pair."
        ),
        figsize=(8, 7),
        draw=draw,
    )


def _scene_triples_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (3, 4))
    planes = random_general_position_planes(rng, n)
    events = generate_triple_events(planes)
    zs = ", ".join(str(event.z) for event in events) or "none"

    def draw(fig: Figure) -> None:
        _draw_triples_scene(fig, planes, f"Step 5.2: n={n} triples")

    return Scene(
        name="triple_events_random",
        title="Step 5.2: triple events",
        caption=(
            f"seed={seed} n={n}. {len(events)} red triple marker(s) "
            f"(expected {n * (n - 1) * (n - 2) // 6}) at z = {zs}. Each "
            "marker sits on three planes; each has a z-axis tick at the same "
            "height. Rotate to check."
        ),
        figsize=(8, 7),
        draw=draw,
    )


def _scene_parallel_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (3, 4))
    planes = random_parallel_family(rng, n)

    def draw(fig: Figure) -> None:
        _draw_parallel_scene(fig, planes, f"Step 5.2: n={n} parallel family")

    return Scene(
        name="parallel_family_random",
        title="Step 5.2: parallel family",
        caption=(
            f"seed={seed} n={n}. {n} parallel planes. No intersection lines "
            "and no triple markers. Rotate to confirm they never meet."
        ),
        figsize=(8, 7),
        draw=draw,
    )


def make_phase5_scenes(seed: int, n: int | None = None) -> tuple[Scene, ...]:
    rng = random.Random(seed)
    return (
        _scene_lines_random(rng, seed, n),
        _scene_triples_random(rng, seed, n),
        _scene_parallel_random(rng, seed, n),
    )


def phase5_scenes(
    *, seed: int | None = None, fixtures: bool = False, n: int | None = None
) -> tuple[Scene, ...]:
    if fixtures:
        return PHASE5_SCENES
    return make_phase5_scenes(choose_seed(seed), n=n)
