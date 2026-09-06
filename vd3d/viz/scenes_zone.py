"""Phase 4 review scenes: zone of a query line."""

from __future__ import annotations

import random

from matplotlib.figure import Figure

from vd3d.arrangement2d import build_line_arrangement
from vd3d.arrangement2d.sampling import random_query_missing_vertices, random_simple_lines
from vd3d.geometry import Line2D
from vd3d.viz.plot_arrangement import draw_arrangement_lines, draw_vertices_numbered, setup_axes_2d
from vd3d.viz.plot_zone import (
    draw_crossings_numbered,
    draw_query_line,
    draw_supporting_faces,
    draw_supporting_vertices,
    draw_zone_faces,
    draw_zone_walk_arrows,
    zone_viewing_limit,
)
from vd3d.viz.random_geom import choose_seed, resolve_n
from vd3d.viz.scene import Scene
from vd3d.zone import compute_supporting_line_zone, compute_zone


def _triangle() -> list[Line2D]:
    return [
        Line2D(a=1, b=0, c=0, id=0),
        Line2D(a=0, b=1, c=0, id=1),
        Line2D(a=1, b=1, c=-1, id=2),
    ]


def _two_axes() -> list[Line2D]:
    return [Line2D(a=1, b=0, c=0, id=0), Line2D(a=0, b=1, c=0, id=1)]


def _query_y(value: str | int) -> Line2D:
    return Line2D(a=0, b=-1, c=value)


def _draw_crossings(fig: Figure, lines: list[Line2D], query: Line2D, title: str) -> None:
    arr = build_line_arrangement(lines)
    zone = compute_zone(arr, query)
    lim = zone_viewing_limit(arr, query=query, crossings=zone.crossings)
    ax = fig.add_subplot(111)
    draw_zone_faces(ax, arr, zone, lim=lim)
    draw_arrangement_lines(ax, arr, lim=lim, color="0.25")
    draw_query_line(ax, query, lim=lim)
    draw_vertices_numbered(ax, arr)
    draw_crossings_numbered(ax, zone.crossings)
    setup_axes_2d(
        ax,
        lim=lim,
        title=f"{title}  —  {len(zone.crossings)} crossings, {len(zone.faces)} zone faces",
    )


def _draw_walk(fig: Figure, lines: list[Line2D], query: Line2D, title: str) -> None:
    arr = build_line_arrangement(lines)
    zone = compute_zone(arr, query)
    lim = zone_viewing_limit(arr, query=query, crossings=zone.crossings)
    ax = fig.add_subplot(111)
    draw_zone_faces(ax, arr, zone, lim=lim)
    draw_arrangement_lines(ax, arr, lim=lim, color="0.25")
    draw_query_line(ax, query, lim=lim)
    draw_crossings_numbered(ax, zone.crossings)
    draw_zone_walk_arrows(ax, zone)
    labels = []
    for i, face in enumerate(zone.faces):
        labels.append("unbounded" if face.unbounded else "interior")
    setup_axes_2d(ax, lim=lim, title=f"{title}  —  {' → '.join(labels)}")


def _draw_supporting(fig: Figure, lines: list[Line2D], line_index: int, title: str) -> None:
    arr = build_line_arrangement(lines)
    zone = compute_supporting_line_zone(arr, line_index)
    query = arr.lines[line_index]
    extra = tuple(v.point for v in (*zone.vertices_on_line, *zone.opposite_vertices))
    lim = zone_viewing_limit(arr, query=query, extra_points=extra)
    ax = fig.add_subplot(111)
    draw_supporting_faces(ax, arr, zone, lim=lim)
    draw_arrangement_lines(ax, arr, lim=lim, color="0.25")
    draw_query_line(ax, query, lim=lim)
    draw_supporting_vertices(ax, zone)
    setup_axes_2d(
        ax,
        lim=lim,
        title=(
            f"{title}  —  {len(zone.vertices_on_line)} on L (red), "
            f"{len(zone.opposite_vertices)} opposite (green)"
        ),
    )


PHASE4_SCENES: tuple[Scene, ...] = (
    Scene(
        name="crossings_triangle",
        title="Step 4.1: crossings",
        caption=(
            "What you must see: a bold red query y=1/2. Colored fills are its "
            "zone faces (left unbounded, the triangle, right unbounded). Gray "
            "faces are not in the zone. Two crimson dots: 0 at (0, 1/2) then 1 "
            "at (1/2, 1/2), increasing left→right along L. No mark at a vertex."
        ),
        figsize=(7, 7),
        draw=lambda fig: _draw_crossings(
            fig, _triangle(), _query_y("1/2"), "Step 4.1: crossings of y=1/2"
        ),
    ),
    Scene(
        name="zone_walk_triangle",
        title="Step 4.2: face walk",
        caption=(
            "What you must see: a bold red query y=1/2. The three zone faces "
            "are painted, not gray. U0 left unbounded, 1 the bounded triangle "
            "(filled), U2 right unbounded. Arrows 0→1→2. Sequence unbounded → "
            "interior → unbounded."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_walk(
            fig, _triangle(), _query_y("1/2"), "Step 4.2: zone walk y=1/2"
        ),
    ),
    Scene(
        name="zone_far_from_vertices",
        title="Step 4.3: far from vertices",
        caption=(
            "What you must see: a bold red query y=-2. Painted faces are its "
            "zone, all unbounded, all below the x-axis. The triangle is gray — "
            "L never enters it. Numbered crossings sit on x=0 and on x+y=1, "
            "far from (0,0), (1,0), (0,1)."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_walk(
            fig, _triangle(), _query_y(-2), "Step 4.3: y=-2 far from vertices"
        ),
    ),
    Scene(
        name="zone_near_miss",
        title="Step 4.3: near-miss",
        caption=(
            "What you must see: a bold red query y=1/10. The triangle is "
            "painted as a zone face. L misses (0,0) and (1,0) by 1/10. Walk is "
            "still unbounded → triangle → unbounded. Gray faces are not in "
            "the zone."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_walk(
            fig, _triangle(), _query_y("1/10"), "Step 4.3: y=1/10 near-miss"
        ),
    ),
    Scene(
        name="supporting_two_lines",
        title="Step 4.4: supporting line, two lines",
        caption=(
            "What you must see: all four quadrants are painted — they are the "
            "zone of supporting line y=0, drawn bold red. One red point at the "
            "origin (on L). No green opposite vertices."
        ),
        figsize=(7, 7),
        draw=lambda fig: _draw_supporting(
            fig, _two_axes(), 1, "Step 4.4: supporting y=0, two lines"
        ),
    ),
    Scene(
        name="supporting_triangle",
        title="Step 4.4: supporting line, triangle",
        caption=(
            "What you must see: bold red L is y=0. Painted faces are those "
            "incident to L (including the triangle). Faces that do not touch L "
            "are gray. Red on-line vertices at (0,0) and (1,0). Green opposite "
            "vertex at (0,1), the far corner of faces touching L."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_supporting(
            fig, _triangle(), 1, "Step 4.4: supporting y=0 of the triangle"
        ),
    ),
    Scene(
        name="zone_five_lines",
        title="Step 4.2: five lines (Test 13)",
        caption=(
            "What you must see: five lines, a bold red query, and the zone "
            "faces painted in walk order. Gray faces are not crossed by L. "
            "Numbered dots sit on the red query and on arrangement edges."
        ),
        figsize=(8, 8),
        draw=lambda fig: _draw_five(fig),
    ),
)


def _draw_five(fig: Figure) -> None:
    rng = random.Random(13)
    lines = random_simple_lines(rng, 5)
    arr = build_line_arrangement(lines)
    query = random_query_missing_vertices(rng, arr)
    _draw_walk(fig, lines, query, "Step 4.2: random n=5, seed=13")


def _scene_crossings_random(
    rng: random.Random, seed: int, n: int | None = None
) -> Scene:
    n = resolve_n(n, rng, (3, 4, 5))
    lines = random_simple_lines(rng, n)
    arr = build_line_arrangement(lines)
    query = random_query_missing_vertices(rng, arr)

    def draw(fig: Figure) -> None:
        _draw_crossings(fig, lines, query, f"Step 4.1: n={n} crossings")

    return Scene(
        name="crossings_random",
        title="Step 4.1: crossings",
        caption=(
            f"seed={seed} n={n}. Bold red query L. Painted faces are its zone. Gray "
            "faces are not in the zone. Numbered crossings increase along L "
            "and sit on arrangement edges, none at a vertex."
        ),
        figsize=(7.5, 7.5),
        draw=draw,
    )


def _scene_walk_random(rng: random.Random, seed: int, n: int | None = None) -> Scene:
    n = resolve_n(n, rng, (3, 4, 5))
    lines = random_simple_lines(rng, n)
    arr = build_line_arrangement(lines)
    query = random_query_missing_vertices(rng, arr)
    zone = compute_zone(arr, query)

    def draw(fig: Figure) -> None:
        _draw_walk(fig, lines, query, f"Step 4.2: n={n} walk")

    return Scene(
        name="zone_walk_random",
        title="Step 4.2: face walk",
        caption=(
            f"seed={seed} n={n}. Bold red query. {len(zone.faces)} painted zone "
            "faces in walk order, arrows between consecutive sample points. "
            "Gray faces are not crossed by L."
        ),
        figsize=(8, 8),
        draw=draw,
    )


def _scene_supporting_random(
    rng: random.Random, seed: int, n: int | None = None
) -> Scene:
    n = resolve_n(n, rng, (2, 3, 4))
    lines = random_simple_lines(rng, n)
    line_index = rng.randrange(n)

    def draw(fig: Figure) -> None:
        _draw_supporting(fig, lines, line_index, f"Step 4.4: supporting line {line_index}")

    zone = compute_supporting_line_zone(build_line_arrangement(lines), line_index)
    return Scene(
        name="supporting_random",
        title="Step 4.4: supporting line",
        caption=(
            f"seed={seed} n={n}. Painted faces are incident to the bold red "
            f"supporting line. Gray faces do not touch L. Red = vertices on L. "
            f"Green = opposite vertices ({len(zone.opposite_vertices)})."
        ),
        figsize=(7.5, 7.5),
        draw=draw,
    )


def make_phase4_scenes(seed: int, n: int | None = None) -> tuple[Scene, ...]:
    rng = random.Random(seed)
    return (
        _scene_crossings_random(rng, seed, n=n),
        _scene_walk_random(rng, seed, n=n),
        _scene_supporting_random(rng, seed, n=n),
    )


def phase4_scenes(
    *, seed: int | None = None, fixtures: bool = False, n: int | None = None
) -> tuple[Scene, ...]:
    if fixtures:
        return PHASE4_SCENES
    return make_phase4_scenes(choose_seed(seed), n=n)
