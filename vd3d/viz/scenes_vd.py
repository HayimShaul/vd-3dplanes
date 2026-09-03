"""Phase 3 review scenes: 2D vertical decomposition."""

from __future__ import annotations

import random

from matplotlib.figure import Figure

from vd3d.arrangement2d import build_line_arrangement
from vd3d.arrangement2d.sampling import random_simple_lines
from vd3d.geometry import Line2D
from vd3d.vertical_decomposition import (
    compute_vertical_decomposition,
    first_hit,
    vertical_rays_from,
)
from vd3d.viz.plot_arrangement import (
    draw_arrangement_lines,
    draw_vertices_numbered,
    setup_axes_2d,
    viewing_limit,
)
from vd3d.viz.plot_vd import (
    draw_cells_filled,
    draw_vertical_rays,
    draw_vertical_walls,
    vd_viewing_limit,
)
from vd3d.viz.random_geom import choose_seed
from vd3d.viz.scene import Scene


def _two_diagonals() -> list[Line2D]:
    return [Line2D(a=1, b=-1, c=0, id=0), Line2D(a=1, b=1, c=0, id=1)]


def _triangle() -> list[Line2D]:
    return [
        Line2D(a=1, b=0, c=0, id=0),
        Line2D(a=0, b=1, c=0, id=1),
        Line2D(a=1, b=1, c=-1, id=2),
    ]


def _rays_and_hits(arrangement):
    out = []
    for vertex in arrangement.vertices:
        for ray in vertical_rays_from(vertex):
            out.append((ray, first_hit(arrangement, ray)))
    return out


def _draw_rays(fig: Figure, lines: list[Line2D], title: str) -> None:
    arr = build_line_arrangement(lines)
    lim = viewing_limit(arr, margin=1.8, minimum=2.8)
    ax = fig.add_subplot(111)
    draw_arrangement_lines(ax, arr, lim=lim, color="0.15")
    rays = _rays_and_hits(arr)
    draw_vertical_rays(ax, rays, lim=lim)
    draw_vertices_numbered(ax, arr)
    n_hit = sum(1 for _ray, hit in rays if not hit.unbounded)
    n_miss = sum(1 for _ray, hit in rays if hit.unbounded)
    setup_axes_2d(
        ax,
        lim=lim,
        title=f"{title}  —  {n_hit} hits (red), {n_miss} unbounded (to window)",
    )


def _draw_cells(fig: Figure, lines: list[Line2D], title: str) -> None:
    arr = build_line_arrangement(lines)
    vd = compute_vertical_decomposition(arr)
    lim = vd_viewing_limit(vd)
    ax = fig.add_subplot(111)
    draw_cells_filled(ax, vd, lim=lim)
    draw_arrangement_lines(ax, arr, lim=lim, color="0.15")
    draw_vertical_walls(ax, vd, lim=lim)
    draw_vertices_numbered(ax, arr)
    setup_axes_2d(
        ax,
        lim=lim,
        title=f"{title}  —  {len(vd.cells)} cells, {len(vd.walls)} dashed vertical walls",
    )


PHASE3_SCENES: tuple[Scene, ...] = (
    Scene(
        name="rays_two_lines",
        title="Step 3.1: rays, two lines",
        caption=(
            "What you must see: y=x and y=-x cross at vertex 0. Two cyan vertical "
            "rays leave the origin, +y and -y, and reach the window edge. No red "
            "hit mark: the two input lines go through the vertex, they are not "
            "hits above or below. No cyan ray crosses a black line without a hit."
        ),
        figsize=(7, 7),
        draw=lambda fig: _draw_rays(fig, _two_diagonals(), "Step 3.1: ±y rays from the crossing"),
    ),
    Scene(
        name="rays_triangle",
        title="Step 3.1: rays, triangle",
        caption=(
            "What you must see: cyan vertical rays from the three vertices. From "
            "(0,0) the +y ray hits (0,1) with a red mark; from (0,1) the -y ray "
            "hits (0,0). From (1,0) both rays miss and go to the window edge. "
            "No ray crosses a black line without a red hit."
        ),
        figsize=(7, 7),
        draw=lambda fig: _draw_rays(fig, _triangle(), "Step 3.1: ±y rays on the triangle"),
    ),
    Scene(
        name="cells_two_lines",
        title="Step 3.2: cells, two lines",
        caption=(
            "What you must see: six colored trapezoids. Dashed red vertical walls "
            "through the origin split the top and bottom wedges; the left and "
            "right wedges are unsplit. Walls are strictly vertical (x=const). "
            "Title count is 6 cells, 2 walls. No visible gap or overlap."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_cells(fig, _two_diagonals(), "Step 3.2: VD of y=x and y=-x"),
    ),
    Scene(
        name="cells_triangle",
        title="Step 3.2–3.3: triangle VD",
        caption=(
            "What you must see: nine colored cells. The bounded triangle is one "
            "cell (a degenerate trapezoid with a vertical left side at x=0 and a "
            "point at (1,0)). Dashed red walls at x=0 (existing) and x=1 (new, "
            "through (1,0) to the window). Trapezoids look like trapezoids; walls "
            "are strictly vertical; no gaps or overlaps. This picture is the "
            "definition of correct for the checked-in triangle JSON fixture."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_cells(fig, _triangle(), "Step 3.3: triangle vertical decomposition"),
    ),
    Scene(
        name="vd_random_n4",
        title="Step 3.3: random n=4",
        caption=(
            "What you must see: a simple arrangement of 4 lines, decomposed into "
            "trapezoids with dashed red vertical walls. Every cell is a trapezoid "
            "(or an unbounded / degenerate one). No leftover sliver, no missing "
            "wedge, walls strictly vertical."
        ),
        figsize=(8, 8),
        draw=lambda fig: _draw_cells(
            fig, random_simple_lines(random.Random(0), 4), "Step 3.3: random n=4, seed=0"
        ),
    ),
    Scene(
        name="vd_random_n6",
        title="Step 3.3: random n=6",
        caption=(
            "What you must see: a simple arrangement of 6 lines with vertical "
            "walls dropped from every vertex. Colored trapezoids fill the window, "
            "dashed red walls are vertical, cell ids sit inside their cells. No "
            "gap or overlap you can see."
        ),
        figsize=(8, 8),
        draw=lambda fig: _draw_cells(
            fig, random_simple_lines(random.Random(1), 6), "Step 3.3: random n=6, seed=1"
        ),
    ),
)


def _scene_rays_random(rng: random.Random, seed: int) -> Scene:
    n = rng.choice([2, 3])
    lines = random_simple_lines(rng, n)

    def draw(fig: Figure) -> None:
        _draw_rays(fig, lines, f"Step 3.1: n={n} rays")

    return Scene(
        name="rays_random",
        title="Step 3.1: rays",
        caption=(
            f"seed={seed}. Cyan ±y rays from every vertex. A red mark is the first "
            "hit; otherwise the ray goes to the window edge. No ray may cross a "
            "black line without a hit."
        ),
        figsize=(7, 7),
        draw=draw,
    )


def _scene_cells_random(rng: random.Random, seed: int) -> Scene:
    n = rng.choice([3, 4, 5])
    lines = random_simple_lines(rng, n)
    arr = build_line_arrangement(lines)
    vd = compute_vertical_decomposition(arr)

    def draw(fig: Figure) -> None:
        _draw_cells(fig, lines, f"Step 3.2: n={n} cells")

    return Scene(
        name="cells_random",
        title="Step 3.2: cells",
        caption=(
            f"seed={seed}. {len(vd.cells)} cells, {len(vd.walls)} dashed vertical "
            "walls. Trapezoids, strictly vertical walls, no visible gap or overlap."
        ),
        figsize=(8, 8),
        draw=draw,
    )


def make_phase3_scenes(seed: int) -> tuple[Scene, ...]:
    rng = random.Random(seed)
    return (
        _scene_rays_random(rng, seed),
        _scene_cells_random(rng, seed),
    )


def phase3_scenes(*, seed: int | None = None, fixtures: bool = False) -> tuple[Scene, ...]:
    if fixtures:
        return PHASE3_SCENES
    return make_phase3_scenes(choose_seed(seed))
