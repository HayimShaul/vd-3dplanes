"""Phase 8 review scenes: cell matching, snapshots, and 3D cells."""

from __future__ import annotations

import random

from matplotlib.figure import Figure

from vd3d.events import events_have_unique_z, generate_all_events
from vd3d.events.samples import (
    planes_alignment_at_z2,
    planes_one,
    planes_through_123,
    planes_two,
)
from vd3d.events.types import EventType
from vd3d.geometry.points import Point3D
from vd3d.sweep import (
    SimultaneousEvents,
    compute_vd_around_event,
    locate_cell3d,
    vertical_decomposition_3d,
)
from vd3d.viz.plot_sweep import (
    draw_matched_vd_pair,
    draw_query_points_3d,
    draw_snapshot_caption,
)
from vd3d.viz.random_geom import choose_seed, random_general_position_planes, resolve_n
from vd3d.viz.scene import Scene


def _planes_unique_z(rng: random.Random, n: int):
    for _ in range(40):
        planes = random_general_position_planes(rng, n)
        if events_have_unique_z(generate_all_events(planes)):
            return planes
    return random_general_position_planes(rng, n)


def _draw_match_around_first_event(fig: Figure, planes, title: str, *, isolate: bool = False) -> None:
    events = generate_all_events(planes)
    if not events:
        raise RuntimeError("expected at least one event")
    if isolate:
        event = events[0]
        around = (event,)
        if any(e.type is EventType.VERTICAL_ALIGNMENT for e in events):
            event = next(e for e in events if e.type is EventType.VERTICAL_ALIGNMENT)
            around = (event,)
    else:
        event = events[0]
        around = events
    vd_before, vd_after, z_minus, z_plus = compute_vd_around_event(planes, event, around)
    draw_matched_vd_pair(
        fig,
        vd_before,
        vd_after,
        title=f"{title}  {event.type.name} z={event.z}",
        left_title=f"z− = {z_minus}",
        right_title=f"z+ = {z_plus}",
    )


def _draw_triple_match(fig: Figure) -> None:
    _draw_match_around_first_event(fig, planes_through_123(), "Step 8.1: triple match")


def _draw_alignment_match(fig: Figure) -> None:
    _draw_match_around_first_event(
        fig, planes_alignment_at_z2(), "Step 8.1: alignment match", isolate=True
    )


def _draw_triple_snapshot(fig: Figure) -> None:
    planes = planes_through_123()
    result = vertical_decomposition_3d(planes)
    event = result.events[0]
    vd_before, vd_after, z_minus, z_plus = compute_vd_around_event(
        planes, event, result.events
    )
    draw_matched_vd_pair(
        fig,
        vd_before,
        vd_after,
        title=f"Step 8.3: {draw_snapshot_caption(result.snapshots[0])}",
        left_title=f"z− = {z_minus}",
        right_title=f"z+ = {z_plus}",
    )


def _draw_query_scene(fig: Figure, planes, title: str, *, seed: int = 8) -> None:
    try:
        result = vertical_decomposition_3d(planes)
    except SimultaneousEvents:
        ax = fig.add_subplot(111)
        ax.set_title(f"{title}  —  simultaneous events (not drawn)")
        ax.set_axis_off()
        return
    rng = random.Random(seed)
    points = [
        Point3D(rng.randint(-2, 2), rng.randint(-2, 2), rng.randint(-2, 2))
        for _ in range(24)
    ]
    cell_ids = [locate_cell3d(result, point) for point in points]
    ax = fig.add_subplot(111, projection="3d")
    n_cells = len(result.cells)
    floors = sum(cell.floor is not None for cell in result.cells)
    ceils = sum(cell.ceiling is not None for cell in result.cells)
    walls = max((len(cell.vertical_walls) for cell in result.cells), default=0)
    draw_query_points_3d(ax, result, points, cell_ids, lim=3.5)
    ax.set_title(
        f"{title}  —  {n_cells} 3D cells, {floors} floors, {ceils} ceilings, "
        f"max walls={walls} (drag to rotate)"
    )


PHASE8_SCENES: tuple[Scene, ...] = (
    Scene(
        name="match_triple",
        title="Step 8.1: match across a triple",
        caption=(
            "What you must see: side-by-side 2D VDs of the Test 21 triple, "
            "just below and just above z=3. Matched cells share a colour. "
            "Unmatched cells (the event neighbourhood) are black. Away from "
            "the small triangle the colours agree; only the local strip recolors."
        ),
        figsize=(11, 5.6),
        draw=_draw_triple_match,
    ),
    Scene(
        name="match_alignment",
        title="Step 8.1: match across an alignment",
        caption=(
            "What you must see: side-by-side VDs around the z=2 alignment. "
            "Matched cells share a colour; any local change is black. Far from "
            "the shared vertical wall the colours agree."
        ),
        figsize=(11, 5.6),
        draw=_draw_alignment_match,
    ),
    Scene(
        name="snapshot_triple",
        title="Step 8.3: event snapshot",
        caption=(
            "What you must see: the same before/after pair as the triple match, "
            "with the snapshot counts in the title (continued / ended / started). "
            "Active 3D cells equal the after 2D cell count. Combinatorics change "
            "only in the black neighbourhood."
        ),
        figsize=(11, 5.6),
        draw=_draw_triple_snapshot,
    ),
    Scene(
        name="cells3d_one_plane",
        title="Step 8.4: one plane, two 3D cells",
        caption=(
            "What you must see: the plane y+z=0 and two colours of query "
            "points, one on each side. No vertical walls. Rotate: a point's "
            "colour is the 3D cell below or above the plane."
        ),
        figsize=(8, 7),
        draw=lambda fig: _draw_query_scene(fig, planes_one(), "Step 8.4: Test 19"),
    ),
    Scene(
        name="cells3d_two_planes",
        title="Step 8.4: two planes, four 3D cells",
        caption=(
            "What you must see: two planes and query points in four colours "
            "(the four quadrants of the crossing). Each cell has at most one "
            "vertical wall. Rotate until the crease of the two planes is obvious."
        ),
        figsize=(8, 7),
        draw=lambda fig: _draw_query_scene(fig, planes_two(), "Step 8.4: Test 20"),
    ),
    Scene(
        name="query_points_triple",
        title="Step 8.5: query points in 3D cells",
        caption=(
            "What you must see: the Test 21 planes and a handful of query "
            "points. Each interior point is coloured by its 3D cell. A wrong "
            "match would put a point on the wrong colour. Boundary points "
            "(on a plane or at an event z) are black."
        ),
        figsize=(8, 7),
        draw=lambda fig: _draw_query_scene(fig, planes_through_123(), "Step 8.5: queries"),
    ),
)


def _scene_match_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (3, 4))
    planes = _planes_unique_z(rng, n)
    events = generate_all_events(planes)

    def draw(fig: Figure) -> None:
        if not events:
            _draw_query_scene(fig, planes, f"Step 8.4: n={n} (no events)", seed=seed)
            return
        try:
            _draw_match_around_first_event(fig, planes, f"Step 8.1: n={n}")
        except (SimultaneousEvents, RuntimeError):
            _draw_query_scene(fig, planes, f"Step 8.4: n={n}", seed=seed)

    if events:
        caption = (
            f"seed={seed} n={n}. Before/after the first event at z={events[0].z}. "
            "Matched cells share a colour; unmatched (local) cells are black."
        )
    else:
        caption = f"seed={seed} n={n}. No events. Query points coloured by 3D cell."

    return Scene(
        name="match_random",
        title="Step 8.1: match across an event",
        caption=caption,
        figsize=(11, 5.6),
        draw=draw,
    )


def _scene_query_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (3, 4))
    planes = _planes_unique_z(rng, n)

    def draw(fig: Figure) -> None:
        _draw_query_scene(fig, planes, f"Step 8.5: n={n} queries", seed=seed)

    return Scene(
        name="query_random",
        title="Step 8.5: query points",
        caption=(
            f"seed={seed} n={n}. Query points coloured by 3D cell. Rotate: a "
            "point on the wrong side of a plane is a matching bug."
        ),
        figsize=(8, 7),
        draw=draw,
    )


def make_phase8_scenes(seed: int, n: int | None = None) -> tuple[Scene, ...]:
    rng = random.Random(seed)
    return (
        _scene_match_random(rng, seed, n),
        _scene_query_random(rng, seed, n),
    )


def phase8_scenes(
    *, seed: int | None = None, fixtures: bool = False, n: int | None = None
) -> tuple[Scene, ...]:
    if fixtures:
        return PHASE8_SCENES
    return make_phase8_scenes(choose_seed(seed), n=n)
