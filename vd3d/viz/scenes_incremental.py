"""Phase 9 review scenes: incremental triple and alignment updates."""

from __future__ import annotations

import random

from matplotlib.figure import Figure

from vd3d.events import events_have_unique_z, generate_all_events
from vd3d.events.samples import planes_alignment_at_z2, planes_through_123
from vd3d.events.types import EventType
from vd3d.sweep import (
    SimultaneousEvents,
    compute_vd_around_event,
    update_2d_decomposition,
)
from vd3d.viz.plot_sweep import draw_matched_vd_pair
from vd3d.viz.random_geom import choose_seed, random_general_position_planes, resolve_n
from vd3d.viz.scene import Scene


def _planes_unique_z(rng: random.Random, n: int):
    for _ in range(40):
        planes = random_general_position_planes(rng, n)
        if events_have_unique_z(generate_all_events(planes)):
            return planes
    return random_general_position_planes(rng, n)


def _around(planes, *, isolate_alignment: bool = False):
    events = generate_all_events(planes)
    if not events:
        raise RuntimeError("expected at least one event")
    if isolate_alignment:
        event = next(e for e in events if e.type is EventType.VERTICAL_ALIGNMENT)
        around = (event,)
    else:
        event = events[0]
        around = events
    vd_before, vd_ref, z_minus, z_plus = compute_vd_around_event(
        planes, event, around
    )
    vd_after = update_2d_decomposition(vd_before, event, planes, z_plus)
    return event, vd_before, vd_after, vd_ref, z_minus, z_plus


def _draw_before_after(fig: Figure, planes, title: str, *, isolate_alignment: bool = False) -> None:
    event, vd_before, vd_after, _ref, z_minus, z_plus = _around(
        planes, isolate_alignment=isolate_alignment
    )
    draw_matched_vd_pair(
        fig,
        vd_before,
        vd_after,
        title=f"{title}  {event.type.name} z={event.z}",
        left_title=f"z− = {z_minus}  (recompute)",
        right_title=f"z+ = {z_plus}  (incremental)",
    )


def _draw_vs_ref(fig: Figure, planes, title: str, *, isolate_alignment: bool = False) -> None:
    event, _before, vd_after, vd_ref, _z_minus, z_plus = _around(
        planes, isolate_alignment=isolate_alignment
    )
    draw_matched_vd_pair(
        fig,
        vd_after,
        vd_ref,
        title=f"{title}  {event.type.name} z+={z_plus}",
        left_title="incremental",
        right_title="recompute (oracle)",
    )


def _draw_triple_update(fig: Figure) -> None:
    _draw_before_after(fig, planes_through_123(), "Step 9.1: triple update")


def _draw_triple_vs_ref(fig: Figure) -> None:
    _draw_vs_ref(fig, planes_through_123(), "Step 9.1: triple vs oracle")


def _draw_alignment_update(fig: Figure) -> None:
    _draw_before_after(
        fig, planes_alignment_at_z2(), "Step 9.2: alignment update", isolate_alignment=True
    )


def _draw_alignment_vs_ref(fig: Figure) -> None:
    _draw_vs_ref(
        fig, planes_alignment_at_z2(), "Step 9.2: alignment vs oracle", isolate_alignment=True
    )


PHASE9_SCENES: tuple[Scene, ...] = (
    Scene(
        name="triple_update",
        title="Step 9.1: incremental triple",
        caption=(
            "What you must see: side-by-side 2D VDs of the Test 21 triple, "
            "just below (recomputed) and just above (incremental update). "
            "Matched cells share a colour. Unmatched cells are black only "
            "near the small triangle; far colours agree."
        ),
        figsize=(11, 5.6),
        draw=_draw_triple_update,
    ),
    Scene(
        name="triple_vs_oracle",
        title="Step 9.1: triple incremental vs recompute",
        caption=(
            "What you must see: the incremental z+ VD next to a from-scratch "
            "compute_vd_at_z at the same height. Every cell matches (same "
            "colours, no black). If any cell is black, the triple handler "
            "disagrees with the oracle."
        ),
        figsize=(11, 5.6),
        draw=_draw_triple_vs_ref,
    ),
    Scene(
        name="alignment_update",
        title="Step 9.2: incremental alignment",
        caption=(
            "What you must see: side-by-side VDs around the z=2 alignment. "
            "The right panel is the incremental update. Matched cells share "
            "a colour; any local change is black. Far from the shared "
            "vertical wall the colours agree."
        ),
        figsize=(11, 5.6),
        draw=_draw_alignment_update,
    ),
    Scene(
        name="alignment_vs_oracle",
        title="Step 9.2: alignment incremental vs recompute",
        caption=(
            "What you must see: incremental z+ next to the recomputed oracle "
            "at the same height. Every cell matches (same colours, no black). "
            "A black cell means the alignment handler is wrong."
        ),
        figsize=(11, 5.6),
        draw=_draw_alignment_vs_ref,
    ),
)


def _scene_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (3, 4))
    planes = _planes_unique_z(rng, n)
    events = generate_all_events(planes)

    def draw(fig: Figure) -> None:
        if not events:
            ax = fig.add_subplot(111)
            ax.set_title(f"Step 9: n={n}  —  no events")
            ax.set_axis_off()
            return
        try:
            _draw_before_after(fig, planes, f"Step 9: n={n}")
        except (SimultaneousEvents, RuntimeError, StopIteration, ValueError):
            ax = fig.add_subplot(111)
            ax.set_title(f"Step 9: n={n}  —  not drawn")
            ax.set_axis_off()

    if events:
        caption = (
            f"seed={seed} n={n}. Before/after the first event at z={events[0].z}. "
            "Right panel is the incremental update. Matched cells share a "
            "colour; unmatched (local) cells are black."
        )
    else:
        caption = f"seed={seed} n={n}. No events."

    return Scene(
        name="update_random",
        title="Step 9: incremental update",
        caption=caption,
        figsize=(11, 5.6),
        draw=draw,
    )


def make_phase9_scenes(seed: int, n: int | None = None) -> tuple[Scene, ...]:
    rng = random.Random(seed)
    return (_scene_random(rng, seed, n),)


def phase9_scenes(
    *, seed: int | None = None, fixtures: bool = False, n: int | None = None
) -> tuple[Scene, ...]:
    if fixtures:
        return PHASE9_SCENES
    return make_phase9_scenes(choose_seed(seed), n=n)
