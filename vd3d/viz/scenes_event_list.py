"""Phase 7 review scenes: event timeline and VD before/after a triple."""

from __future__ import annotations

import random

from matplotlib.figure import Figure

from vd3d.events import (
    EventType,
    choose_z_below_all_events,
    compute_vd_at_z,
    generate_all_events,
)
from vd3d.events.samples import planes_alignment_at_z2, planes_through_123
from vd3d.viz.plot_event_list import draw_event_timeline, draw_vd_slice, shared_vd_limit
from vd3d.viz.random_geom import choose_seed, random_general_position_planes, resolve_n
from vd3d.viz.scene import Scene


def _draw_timeline(fig: Figure, planes, title: str) -> None:
    events = generate_all_events(planes)
    ax = fig.add_subplot(111)
    n_triple = sum(event.type is EventType.TRIPLE_INTERSECTION for event in events)
    n_align = sum(event.type is EventType.VERTICAL_ALIGNMENT for event in events)
    draw_event_timeline(
        ax,
        events,
        title=f"{title}  —  {n_triple} triple, {n_align} align",
    )


def _draw_vd_before_after(fig: Figure, planes, title: str) -> None:
    events = generate_all_events(planes)
    triples = [event for event in events if event.type is EventType.TRIPLE_INTERSECTION]
    if not triples:
        raise RuntimeError("expected a triple event for the before/after figure")
    event = triples[0]
    z_minus = event.z - 1
    z_plus = event.z + 1
    vd_before = compute_vd_at_z(planes, z_minus)
    vd_after = compute_vd_at_z(planes, z_plus)
    lim = shared_vd_limit(vd_before, vd_after)
    ax_l = fig.add_subplot(121)
    ax_r = fig.add_subplot(122)
    draw_vd_slice(ax_l, vd_before, lim=lim, title=f"z− = {z_minus}")
    draw_vd_slice(ax_r, vd_after, lim=lim, title=f"z+ = {z_plus}")
    fig.suptitle(f"{title}  —  triple at z={event.z}")


PHASE7_SCENES: tuple[Scene, ...] = (
    Scene(
        name="event_timeline",
        title="Step 7.1: event timeline",
        caption=(
            "What you must see: a 1D z-axis. Red ticks are triple events; "
            "green ticks are vertical alignments. Each tick is labeled by "
            "type and z. On this fixture the green alignment sits at z=2, "
            "between red triples. No two labels share a sort key."
        ),
        figsize=(9, 3.6),
        draw=lambda fig: _draw_timeline(
            fig, planes_alignment_at_z2(), "Step 7.1: all events"
        ),
    ),
    Scene(
        name="vd_before_after_triple",
        title="Step 7.2: VD at z− / z+",
        caption=(
            "What you must see: side-by-side 2D VDs of the Test 15 triple. "
            "Left is just below z=3, right just above. Each side has a small "
            "triangle (3 vertices) whose vertical walls are dashed red. The "
            "triangle flips through the concurrence at (1, 2); that is the "
            "only combinatorial change. Vertex counts stay 3 and 3."
        ),
        figsize=(11, 5.6),
        draw=lambda fig: _draw_vd_before_after(
            fig, planes_through_123(), "Step 7.2: triple before / after"
        ),
    ),
)


def _scene_timeline_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (3, 4))
    planes = random_general_position_planes(rng, n)
    events = generate_all_events(planes)
    zs = ", ".join(f"{event.z}" for event in events) or "none"

    def draw(fig: Figure) -> None:
        _draw_timeline(fig, planes, f"Step 7.1: n={n} timeline")

    return Scene(
        name="event_timeline_random",
        title="Step 7.1: event timeline",
        caption=(
            f"seed={seed} n={n}. {len(events)} event(s) at z = {zs}. Red = "
            "triple, green = alignment. Ticks must increase in z from left "
            "to right."
        ),
        figsize=(9, 3.6),
        draw=draw,
    )


def _scene_vd_random(rng: random.Random, seed: int, n: int | None) -> Scene:
    n = resolve_n(n, rng, (3, 4))
    planes = random_general_position_planes(rng, n)
    events = generate_all_events(planes)
    triples = [event for event in events if event.type is EventType.TRIPLE_INTERSECTION]
    z0 = choose_z_below_all_events(events)

    def draw(fig: Figure) -> None:
        if triples:
            _draw_vd_before_after(fig, planes, f"Step 7.2: n={n} triple")
            return
        vd = compute_vd_at_z(planes, z0)
        ax = fig.add_subplot(111)
        draw_vd_slice(ax, vd, lim=shared_vd_limit(vd), title=f"z0={z0} (no triples)")

    if triples:
        caption = (
            f"seed={seed} n={n}. Side-by-side VD around the first triple at "
            f"z={triples[0].z}. The small triangle (or its absence at a "
            "degeneracy) is the local change; walls stay vertical."
        )
    else:
        caption = (
            f"seed={seed} n={n}. No triples. The single VD is the initial "
            f"slice at z={z0}."
        )

    return Scene(
        name="vd_before_after_random",
        title="Step 7.2: VD at z− / z+",
        caption=caption,
        figsize=(11, 5.6),
        draw=draw,
    )


def make_phase7_scenes(seed: int, n: int | None = None) -> tuple[Scene, ...]:
    rng = random.Random(seed)
    return (
        _scene_timeline_random(rng, seed, n),
        _scene_vd_random(rng, seed, n),
    )


def phase7_scenes(
    *, seed: int | None = None, fixtures: bool = False, n: int | None = None
) -> tuple[Scene, ...]:
    if fixtures:
        return PHASE7_SCENES
    return make_phase7_scenes(choose_seed(seed), n=n)
