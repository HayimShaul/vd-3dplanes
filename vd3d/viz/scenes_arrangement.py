"""Phase 2 review scenes: 2D line arrangements."""

from __future__ import annotations

import random

from matplotlib.figure import Figure

from vd3d.arrangement2d import (
    COINCIDENT,
    PARALLEL,
    build_line_arrangement,
    intersect_lines_2d,
)
from vd3d.arrangement2d.sampling import random_simple_lines
from vd3d.geometry import Line2D, Point2D
from vd3d.viz.convert import to_float
from vd3d.viz.plot_arrangement import (
    draw_edge_pieces,
    draw_faces_filled,
    draw_outgoing_half_edges,
    draw_vertices_numbered,
    format_line2d,
    setup_axes_2d,
    viewing_limit,
)
from vd3d.viz.plot_geometry import draw_line2d
from vd3d.viz.random_geom import choose_seed
from vd3d.viz.scene import Scene


def _axes() -> list[Line2D]:
    return [Line2D(a=1, b=0, c=0, id=0), Line2D(a=0, b=1, c=0, id=1)]


def _triangle() -> list[Line2D]:
    return [
        Line2D(a=1, b=0, c=0, id=0),
        Line2D(a=0, b=1, c=0, id=1),
        Line2D(a=1, b=1, c=-1, id=2),
    ]


def _draw_intersect_pairs(fig: Figure) -> None:
    cases = [
        (
            "crossing",
            Line2D(a=1, b=0, c=0),
            Line2D(a=0, b=1, c=0),
            "x=0 ∩ y=0",
        ),
        (
            "generic",
            Line2D(a=1, b=1, c=-3),
            Line2D(a=1, b=-1, c=-1),
            "x+y=3 ∩ x-y=1",
        ),
        (
            "parallel",
            Line2D(a=1, b=0, c=0),
            Line2D(a=1, b=0, c=-2),
            "x=0 ∥ x=2",
        ),
        (
            "coincident",
            Line2D(a=1, b=0, c=0),
            Line2D(a=2, b=0, c=0),
            "x=0 and 2x=0",
        ),
    ]
    for i, (_name, left, right, title) in enumerate(cases):
        ax = fig.add_subplot(2, 2, i + 1)
        draw_line2d(ax, left, lim=4, color="steelblue", label=format_line2d(left))
        draw_line2d(ax, right, lim=4, color="orange", label=format_line2d(right))
        result = intersect_lines_2d(left, right)
        if result is PARALLEL:
            ax.set_title(f"{title}\nPARALLEL — no mark")
        elif result is COINCIDENT:
            ax.set_title(f"{title}\nCOINCIDENT — no extra mark")
        else:
            ax.scatter(
                [to_float(result.x)],
                [to_float(result.y)],
                c="red",
                s=80,
                zorder=5,
                label=f"({result.x}, {result.y})",
            )
            ax.set_title(f"{title}\nintersection ({result.x}, {result.y})")
        ax.set_xlim(-4, 4)
        ax.set_ylim(-4, 4)
        ax.set_aspect("equal")
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.legend(loc="upper right", fontsize=7)
    fig.suptitle("Step 2.1: line–line intersection (marks sit on both lines; parallels have none)")


def _draw_pieces(fig: Figure, lines: list[Line2D], title: str) -> None:
    arr = build_line_arrangement(lines)
    lim = viewing_limit(arr)
    ax = fig.add_subplot(111)
    draw_edge_pieces(ax, arr, lim=lim, arrows=True)
    draw_vertices_numbered(ax, arr)
    rays = sum(1 for e in arr.edges if e.kind == "ray")
    segs = sum(1 for e in arr.edges if e.kind == "segment")
    n_v = len(arr.vertices)
    v_word = "vertex" if n_v == 1 else "vertices"
    setup_axes_2d(
        ax,
        lim=lim,
        title=f"{title}  —  {n_v} {v_word}, {rays} rays, {segs} segments",
    )


def _draw_circulation(fig: Figure) -> None:
    arr = build_line_arrangement(_axes())
    ax = fig.add_subplot(111)
    draw_outgoing_half_edges(ax, arr, 0, scale=1.3)
    setup_axes_2d(
        ax,
        lim=2.2,
        title="Step 2.3: outgoing half-edges at (0,0), numbered CCW from +x",
    )
    ax.annotate(
        "+x",
        (1.6, 0.05),
        fontsize=9,
        color="0.3",
    )


def _draw_faces(fig: Figure, lines: list[Line2D], title: str) -> None:
    arr = build_line_arrangement(lines)
    lim = viewing_limit(arr, margin=1.8, minimum=2.8)
    ax = fig.add_subplot(111)
    draw_faces_filled(ax, arr, lim=lim)
    for line in arr.lines:
        draw_line2d(ax, line, lim=lim, color="0.15")
    draw_vertices_numbered(ax, arr)
    n_u = len(arr.unbounded_faces)
    n_b = len(arr.bounded_faces)
    setup_axes_2d(
        ax,
        lim=lim,
        title=f"{title}  —  {n_b} bounded, {n_u} unbounded (U), {len(arr.faces)} faces total",
    )


def _draw_random_faces(fig: Figure, seed: int = 0, n: int = 6) -> None:
    rng = random.Random(seed)
    lines = random_simple_lines(rng, n)
    _draw_faces(fig, lines, f"Step 2.5: random simple arrangement n={n}, seed={seed}")


PHASE2_SCENES: tuple[Scene, ...] = (
    Scene(
        name="intersect_pairs",
        title="Step 2.1: intersections",
        caption=(
            "What you must see: four pairs. Top-left: red mark at (0,0) on both axes. "
            "Top-right: red mark at (2,1) on both lines. Bottom-left: two vertical "
            "parallels, no intersection mark. Bottom-right: coincident vertical line "
            "drawn twice, no extra intersection mark."
        ),
        figsize=(9, 8),
        draw=_draw_intersect_pairs,
    ),
    Scene(
        name="pieces_two_lines",
        title="Step 2.2: two lines",
        caption=(
            "What you must see: two lines crossing at vertex 0. Four rays in four "
            "colors, each with an arrow pointing away from the vertex. No interior "
            "segments. Title count is 1 vertex, 4 rays, 0 segments."
        ),
        figsize=(7, 7),
        draw=lambda fig: _draw_pieces(fig, _axes(), "Step 2.2: two intersecting lines"),
    ),
    Scene(
        name="pieces_triangle",
        title="Step 2.2: three lines",
        caption=(
            "What you must see: a triangle with vertices 0,1,2. Three interior "
            "segments (the triangle sides) plus six rays with arrows leaving the "
            "triangle. Title count is 3 vertices, 6 rays, 3 segments."
        ),
        figsize=(7, 7),
        draw=lambda fig: _draw_pieces(fig, _triangle(), "Step 2.2: three lines in general position"),
    ),
    Scene(
        name="circulation_vertex",
        title="Step 2.3: CCW circulation",
        caption=(
            "What you must see: the origin with four outgoing arrows numbered 0,1,2,3. "
            "Numbers increase counterclockwise, starting at +x (arrow 0 points right, "
            "then up, then left, then down). Face-on-the-left follows this CCW order."
        ),
        figsize=(6.5, 6.5),
        draw=_draw_circulation,
    ),
    Scene(
        name="faces_triangle",
        title="Step 2.4: faces",
        caption=(
            "What you must see: one bounded triangular cell with a numeric id at its "
            "interior sample, plus six unbounded outer cells labeled U. Seven colors, "
            "no leftover sliver, no missing wedge. The triangle is the only bounded face."
        ),
        figsize=(7.5, 7.5),
        draw=lambda fig: _draw_faces(fig, _triangle(), "Step 2.4: triangle arrangement"),
    ),
    Scene(
        name="faces_random_n6",
        title="Step 2.5: random n=6",
        caption=(
            "What you must see: a simple arrangement of 6 lines, every wedge filled, "
            "no leftover sliver, no missing cell. Bounded faces have numeric ids; "
            "unbounded faces are labeled U. Vertex/edge/face counts match the formula "
            "V=15, E=36, F=22."
        ),
        figsize=(8, 8),
        draw=lambda fig: _draw_random_faces(fig, seed=0, n=6),
    ),
)


def _scene_intersect_random(rng: random.Random, seed: int) -> Scene:
    kind = rng.choice(["crossing", "parallel"])
    if kind == "parallel":
        a, b = rng.randint(-3, 3), rng.randint(-3, 3)
        if a == 0 and b == 0:
            a = 1
        c1 = rng.randint(-3, 3)
        c2 = c1 + rng.choice([-3, -2, -1, 1, 2, 3])
        left, right = Line2D(a=a, b=b, c=c1), Line2D(a=a, b=b, c=c2)
    else:
        lines = random_simple_lines(rng, 2)
        left, right = lines[0], lines[1]

    def draw(fig: Figure) -> None:
        ax = fig.add_subplot(111)
        draw_line2d(ax, left, lim=5, color="steelblue", label=format_line2d(left))
        draw_line2d(ax, right, lim=5, color="orange", label=format_line2d(right))
        result = intersect_lines_2d(left, right)
        if isinstance(result, Point2D):
            ax.scatter([to_float(result.x)], [to_float(result.y)], c="red", s=90, zorder=5)
            ax.set_title(f"Step 2.1: intersection ({result.x}, {result.y})")
        else:
            ax.set_title("Step 2.1: PARALLEL — no intersection mark")
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_aspect("equal")
        ax.grid(True, linestyle=":", alpha=0.5)
        ax.legend(loc="upper right")

    return Scene(
        name="intersect_random",
        title="Step 2.1: intersection",
        caption=(
            f"seed={seed}. Blue {format_line2d(left)} and orange {format_line2d(right)}. "
            "A red mark must sit on both lines if they cross; parallels have no mark."
        ),
        figsize=(7, 7),
        draw=draw,
    )


def _scene_pieces_random(rng: random.Random, seed: int) -> Scene:
    n = rng.choice([2, 3])
    lines = random_simple_lines(rng, n)
    arr = build_line_arrangement(lines)
    rays = sum(1 for e in arr.edges if e.kind == "ray")
    segs = sum(1 for e in arr.edges if e.kind == "segment")

    def draw(fig: Figure) -> None:
        _draw_pieces(fig, lines, f"Step 2.2: n={n} pieces")

    return Scene(
        name="pieces_random",
        title="Step 2.2: pieces",
        caption=(
            f"seed={seed}. {len(arr.vertices)} vertices, {rays} rays (arrows), "
            f"{segs} segments. Count the colored pieces; they must match the title."
        ),
        figsize=(7, 7),
        draw=draw,
    )


def _scene_circulation_random(rng: random.Random, seed: int) -> Scene:
    lines = random_simple_lines(rng, 2)
    arr = build_line_arrangement(lines)

    def draw(fig: Figure) -> None:
        ax = fig.add_subplot(111)
        draw_outgoing_half_edges(ax, arr, 0, scale=1.3)
        setup_axes_2d(ax, lim=2.4, title="Step 2.3: outgoing CCW at the crossing")

    return Scene(
        name="circulation_random",
        title="Step 2.3: circulation",
        caption=(
            f"seed={seed}. Four outgoing arrows at the crossing, numbered 0–3 "
            "increasing counterclockwise from the +x side."
        ),
        figsize=(6.5, 6.5),
        draw=draw,
    )


def _scene_faces_random(rng: random.Random, seed: int) -> Scene:
    n = rng.choice([3, 4, 6])
    lines = random_simple_lines(rng, n)
    arr = build_line_arrangement(lines)

    def draw(fig: Figure) -> None:
        _draw_faces(fig, lines, f"Step 2.4–2.5: n={n} faces")

    return Scene(
        name="faces_random",
        title="Step 2.4: faces",
        caption=(
            f"seed={seed}. {len(arr.bounded_faces)} bounded + {len(arr.unbounded_faces)} "
            f"unbounded = {len(arr.faces)} faces. Every wedge filled, no leftover sliver."
        ),
        figsize=(8, 8),
        draw=draw,
    )


def make_phase2_scenes(seed: int) -> tuple[Scene, ...]:
    rng = random.Random(seed)
    return (
        _scene_intersect_random(rng, seed),
        _scene_pieces_random(rng, seed),
        _scene_circulation_random(rng, seed),
        _scene_faces_random(rng, seed),
    )


def phase2_scenes(*, seed: int | None = None, fixtures: bool = False) -> tuple[Scene, ...]:
    if fixtures:
        return PHASE2_SCENES
    return make_phase2_scenes(choose_seed(seed))
