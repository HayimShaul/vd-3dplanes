"""Review scenes drawn from the kernel. Used by the GUI and by visual tests."""

from __future__ import annotations

import random
from matplotlib.figure import Figure

from vd3d.geometry import (
    PARALLEL,
    Plane,
    Point3D,
    intersect_planes,
    intersect_three_planes,
    slice_plane_at_z,
)
from vd3d.viz.convert import to_float
from vd3d.viz.plot_geometry import (
    draw_lifted_line2d,
    draw_line2d,
    draw_line3d,
    draw_plane,
    draw_point3d,
    set_equal_3d,
)
from vd3d.viz.random_geom import (
    choose_seed,
    display_t_range,
    format_plane,
    format_point,
    offset_along_normal,
    point_on_plane,
    random_intersecting_planes,
    random_parallel_planes,
    random_plane,
    random_triple_planes,
    view_lim,
)
from vd3d.viz.scene import Scene


def _draw_plane_eval_signs_fixed(fig: Figure) -> None:
    plane = Plane(id=1, a=1, b=1, c=1, d=-5)
    on = Point3D(1, 2, 2)
    above = Point3D(1, 2, 3)
    below = Point3D(1, 2, 1)
    ax = fig.add_subplot(111, projection="3d")
    draw_plane(ax, plane, lim=4, color="0.7", label="x+y+z-5=0")
    draw_point3d(ax, on, color="green", size=80, label=f"on  eval={plane.eval(on)}")
    draw_point3d(ax, above, color="red", size=80, label=f"above  eval={plane.eval(above)}")
    draw_point3d(ax, below, color="blue", size=80, label=f"below  eval={plane.eval(below)}")
    set_equal_3d(ax, lim=4)
    ax.set_title("Step 1.1: plane x+y+z-5=0 with signed sample points")
    ax.legend(loc="upper left")


def _draw_intersecting_planes_fixed(fig: Figure) -> None:
    p = Plane(id=1, a=1, b=0, c=0, d=0)
    q = Plane(id=2, a=0, b=1, c=0, d=0)
    line = intersect_planes(p, q)
    if line is PARALLEL:
        raise RuntimeError("x=0 and y=0 must intersect")
    ax = fig.add_subplot(111, projection="3d")
    draw_plane(ax, p, lim=3, color="steelblue", label="x=0")
    draw_plane(ax, q, lim=3, color="orange", label="y=0")
    draw_line3d(ax, line, t_min=-3, t_max=3, color="magenta", linewidth=3, label="intersection")
    set_equal_3d(ax, lim=3)
    ax.set_title("Step 1.2: x=0 ∩ y=0 is the z-axis")
    ax.legend(loc="upper left")


def _draw_parallel_planes_fixed(fig: Figure) -> None:
    p = Plane(id=1, a=1, b=0, c=0, d=0)
    q = Plane(id=2, a=1, b=0, c=0, d=-2)
    ax = fig.add_subplot(111, projection="3d")
    draw_plane(ax, p, lim=3, color="steelblue", label="x=0")
    draw_plane(ax, q, lim=3, color="orange", label="x=2")
    set_equal_3d(ax, lim=3)
    ax.set_title("Step 1.2: parallel planes x=0 and x=2 — no intersection line")
    ax.legend(loc="upper left")


def _draw_three_planes_fixed(fig: Figure) -> None:
    p = Plane(id=1, a=1, b=0, c=0, d=0)
    q = Plane(id=2, a=0, b=1, c=0, d=0)
    r = Plane(id=3, a=0, b=0, c=1, d=-1)
    point = intersect_three_planes(p, q, r)
    if point is None:
        raise RuntimeError("x=0, y=0, z=1 must meet at a point")
    ax = fig.add_subplot(111, projection="3d")
    draw_plane(ax, p, lim=2.5, color="steelblue", label="x=0")
    draw_plane(ax, q, lim=2.5, color="orange", label="y=0")
    draw_plane(ax, r, lim=2.5, color="seagreen", label="z=1")
    draw_point3d(ax, point, color="red", size=120, label="(0,0,1)")
    ax.text(
        to_float(point.x),
        to_float(point.y),
        to_float(point.z) + 0.25,
        "(0,0,1)",
        color="red",
    )
    set_equal_3d(ax, lim=2.5)
    ax.set_title("Step 1.3: x=0, y=0, z=1 meet at (0,0,1)")
    ax.legend(loc="upper left")


def _draw_slice_fixed(fig: Figure) -> None:
    plane = Plane(id=1, a=1, b=1, c=1, d=-5)
    line = slice_plane_at_z(plane, 2)
    if line is None:
        raise RuntimeError("slanted plane must have a z=2 slice")
    ax3d = fig.add_subplot(121, projection="3d")
    draw_plane(ax3d, plane, lim=4, color="0.7", label="x+y+z-5=0")
    draw_lifted_line2d(ax3d, line, 2, color="red", label="slice at z=2")
    set_equal_3d(ax3d, lim=4)
    ax3d.set_title("3D: plane and its z=2 trace (drag to rotate)")
    ax3d.legend(loc="upper left")

    ax2d = fig.add_subplot(122)
    draw_line2d(ax2d, line, lim=5, color="red", label="x+y-3=0")
    ax2d.scatter([0, 3], [3, 0], c="black", zorder=3)
    ax2d.annotate("(0,3)", (0, 3), textcoords="offset points", xytext=(6, 6))
    ax2d.annotate("(3,0)", (3, 0), textcoords="offset points", xytext=(6, 6))
    ax2d.set_xlim(-1, 5)
    ax2d.set_ylim(-1, 5)
    ax2d.set_aspect("equal")
    ax2d.grid(True, linestyle=":", alpha=0.6)
    ax2d.set_xlabel("x")
    ax2d.set_ylabel("y")
    ax2d.set_title("2D: the same line in xy")
    ax2d.legend()
    fig.suptitle("Step 1.4: slice of x+y+z-5=0 at z=2 is x+y-3=0")


PHASE1_SCENES: tuple[Scene, ...] = (
    Scene(
        name="plane_eval_signs",
        title="Step 1.1: signed residuals",
        caption=(
            "What you must see: the gray plane x+y+z-5=0, and three dots on "
            "the vertical line x=1, y=2. Green at (1,2,2) is on the plane "
            "(eval=0). Red at (1,2,3) is above (eval=1). Blue at (1,2,1) is "
            "below (eval=-1). Rotate until you can see the green dot sitting "
            "in the mesh and the other two off it."
        ),
        figsize=(8, 7),
        draw=_draw_plane_eval_signs_fixed,
    ),
    Scene(
        name="intersect_planes_z_axis",
        title="Step 1.2: intersecting planes",
        caption=(
            "What you must see: a blue plane (x=0) and an orange plane (y=0). "
            "Their intersection is the magenta line along the z-axis. Rotate "
            "around z: the magenta line should stay the crease of the two planes."
        ),
        figsize=(8, 7),
        draw=_draw_intersecting_planes_fixed,
    ),
    Scene(
        name="parallel_planes",
        title="Step 1.2: parallel planes",
        caption=(
            "What you must see: two vertical planes facing each other, blue at "
            "x=0 and orange at x=2. There is no intersection line. Rotate to "
            "confirm they never meet."
        ),
        figsize=(8, 7),
        draw=_draw_parallel_planes_fixed,
    ),
    Scene(
        name="three_planes_point",
        title="Step 1.3: three planes",
        caption=(
            "What you must see: blue x=0, orange y=0, green z=1. The fat red "
            "marker is (0,0,1) and must sit on all three meshes. Rotate until "
            "you can see it on the horizontal plane and on both vertical planes."
        ),
        figsize=(8, 7),
        draw=_draw_three_planes_fixed,
    ),
    Scene(
        name="slice_at_z",
        title="Step 1.4: horizontal slice",
        caption=(
            "What you must see: left, the gray plane with a red line lying in "
            "it at height z=2; right, that same line in xy as x+y-3=0 through "
            "(0,3) and (3,0). Rotate the left view: the red line must stay in "
            "the plane."
        ),
        figsize=(11, 6.5),
        draw=_draw_slice_fixed,
    ),
)


def _scene_eval_signs(rng: random.Random, seed: int) -> Scene:
    plane = random_plane(rng, 1)
    on = point_on_plane(plane, rng)
    above = offset_along_normal(on, plane, 1)
    below = offset_along_normal(on, plane, -1)
    lim = view_lim(on, above, below)
    eq = format_plane(plane)

    def draw(fig: Figure) -> None:
        ax = fig.add_subplot(111, projection="3d")
        draw_plane(ax, plane, lim=lim, color="0.7", label=eq)
        draw_point3d(ax, on, color="green", size=80, label=f"on  eval={plane.eval(on)}")
        draw_point3d(ax, above, color="red", size=80, label=f"above  eval={plane.eval(above)}")
        draw_point3d(ax, below, color="blue", size=80, label=f"below  eval={plane.eval(below)}")
        set_equal_3d(ax, lim=lim)
        ax.set_title(f"Step 1.1: {eq}")
        ax.legend(loc="upper left")

    return Scene(
        name="plane_eval_signs",
        title="Step 1.1: signed residuals",
        caption=(
            f"seed={seed}. Plane {eq}. Green {format_point(on)} is on the plane "
            f"(eval={plane.eval(on)}). Red {format_point(above)} is above "
            f"(eval={plane.eval(above)}). Blue {format_point(below)} is below "
            f"(eval={plane.eval(below)}). Rotate until the green dot sits in "
            "the mesh and the other two are off it."
        ),
        figsize=(8, 7),
        draw=draw,
    )


def _scene_intersecting(rng: random.Random, seed: int) -> Scene:
    p, q = random_intersecting_planes(rng)
    line = intersect_planes(p, q)
    if line is PARALLEL:
        raise RuntimeError("sampled planes were parallel")
    t_min, t_max = display_t_range(line.direction)
    lim = view_lim(line.point_at(0), line.point_at(1), line.point_at(-1))
    peq, qeq = format_plane(p), format_plane(q)

    def draw(fig: Figure) -> None:
        ax = fig.add_subplot(111, projection="3d")
        draw_plane(ax, p, lim=lim, color="steelblue", label=peq)
        draw_plane(ax, q, lim=lim, color="orange", label=qeq)
        draw_line3d(
            ax,
            line,
            t_min=t_min,
            t_max=t_max,
            color="magenta",
            linewidth=3,
            label="intersection",
        )
        set_equal_3d(ax, lim=lim)
        ax.set_title(f"Step 1.2: {peq} ∩ {qeq}")
        ax.legend(loc="upper left")

    return Scene(
        name="intersect_planes",
        title="Step 1.2: intersecting planes",
        caption=(
            f"seed={seed}. Blue {peq} and orange {qeq}. The magenta line is "
            "their intersection. Rotate: it must stay the crease of the two meshes."
        ),
        figsize=(8, 7),
        draw=draw,
    )


def _scene_parallel(rng: random.Random, seed: int) -> Scene:
    p, q = random_parallel_planes(rng)
    peq, qeq = format_plane(p), format_plane(q)
    on_p = point_on_plane(p, rng)
    on_q = point_on_plane(q, rng)
    lim = view_lim(on_p, on_q)

    def draw(fig: Figure) -> None:
        ax = fig.add_subplot(111, projection="3d")
        draw_plane(ax, p, lim=lim, color="steelblue", label=peq)
        draw_plane(ax, q, lim=lim, color="orange", label=qeq)
        set_equal_3d(ax, lim=lim)
        ax.set_title(f"Step 1.2: parallel {peq} and {qeq}")
        ax.legend(loc="upper left")

    return Scene(
        name="parallel_planes",
        title="Step 1.2: parallel planes",
        caption=(
            f"seed={seed}. Blue {peq} and orange {qeq} have parallel normals. "
            "There is no intersection line. Rotate to confirm they never meet."
        ),
        figsize=(8, 7),
        draw=draw,
    )


def _scene_triple(rng: random.Random, seed: int) -> Scene:
    point, p, q, r = random_triple_planes(rng)
    computed = intersect_three_planes(p, q, r)
    if computed is None or computed != point:
        raise RuntimeError("constructed triple did not reconstruct the point")
    lim = view_lim(point, minimum=2.5)
    peq, qeq, req = format_plane(p), format_plane(q), format_plane(r)
    label = format_point(point)

    def draw(fig: Figure) -> None:
        ax = fig.add_subplot(111, projection="3d")
        draw_plane(ax, p, lim=lim, color="steelblue", label=peq)
        draw_plane(ax, q, lim=lim, color="orange", label=qeq)
        draw_plane(ax, r, lim=lim, color="seagreen", label=req)
        draw_point3d(ax, point, color="red", size=120, label=label)
        ax.text(
            to_float(point.x),
            to_float(point.y),
            to_float(point.z) + 0.2 * lim,
            label,
            color="red",
        )
        set_equal_3d(ax, lim=lim)
        ax.set_title(f"Step 1.3: three planes meet at {label}")
        ax.legend(loc="upper left")

    return Scene(
        name="three_planes_point",
        title="Step 1.3: three planes",
        caption=(
            f"seed={seed}. Blue {peq}, orange {qeq}, green {req}. The red "
            f"marker is {label} and must sit on all three meshes. Rotate to check."
        ),
        figsize=(8, 7),
        draw=draw,
    )


def _scene_slice(rng: random.Random, seed: int) -> Scene:
    plane = random_plane(rng, 1, nonzero_xy=True)
    z = rng.randint(-3, 3)
    line = slice_plane_at_z(plane, z)
    if line is None:
        raise RuntimeError("slice of a non-horizontal plane must exist")
    p0, p1 = line.sample_points()
    eq = format_plane(plane)
    line_eq = f"{line.a}x+{line.b}y+{line.c}=0"
    on = point_on_plane(plane, rng)
    lim = view_lim(on, Point3D(p0.x, p0.y, z), Point3D(p1.x, p1.y, z))

    def draw(fig: Figure) -> None:
        ax3d = fig.add_subplot(121, projection="3d")
        draw_plane(ax3d, plane, lim=lim, color="0.7", label=eq)
        draw_lifted_line2d(ax3d, line, z, extent=lim, color="red", label=f"slice at z={z}")
        set_equal_3d(ax3d, lim=lim)
        ax3d.set_title(f"3D: {eq} at z={z} (drag to rotate)")
        ax3d.legend(loc="upper left")

        ax2d = fig.add_subplot(122)
        draw_lim = max(5.0, abs(to_float(p0.x)), abs(to_float(p0.y)), abs(to_float(p1.x)), abs(to_float(p1.y))) + 1
        draw_line2d(ax2d, line, lim=draw_lim, color="red", label=line_eq)
        ax2d.scatter([to_float(p0.x), to_float(p1.x)], [to_float(p0.y), to_float(p1.y)], c="black", zorder=3)
        ax2d.annotate(f"({p0.x},{p0.y})", (to_float(p0.x), to_float(p0.y)), textcoords="offset points", xytext=(6, 6))
        ax2d.annotate(f"({p1.x},{p1.y})", (to_float(p1.x), to_float(p1.y)), textcoords="offset points", xytext=(6, 6))
        ax2d.set_aspect("equal")
        ax2d.grid(True, linestyle=":", alpha=0.6)
        ax2d.set_xlabel("x")
        ax2d.set_ylabel("y")
        ax2d.set_title(f"2D: {line_eq}")
        ax2d.legend()
        fig.suptitle(f"Step 1.4: slice of {eq} at z={z}")

    return Scene(
        name="slice_at_z",
        title="Step 1.4: horizontal slice",
        caption=(
            f"seed={seed}. Left: {eq} with its red trace at z={z}. Right: the "
            f"same line {line_eq} in xy. Rotate the left view: the red line "
            "must stay in the plane."
        ),
        figsize=(11, 6.5),
        draw=draw,
    )


def make_phase1_scenes(seed: int, n: int | None = None) -> tuple[Scene, ...]:
    """Five Phase 1 scenes from ``seed``. Same seed → same geometry.

    ``n`` is accepted for the viewer API. Phase 1 scenes have a fixed number
    of planes (1, 2, or 3), so ``n`` does not change them.
    """
    rng = random.Random(seed)
    _ = n
    return (
        _scene_eval_signs(rng, seed),
        _scene_intersecting(rng, seed),
        _scene_parallel(rng, seed),
        _scene_triple(rng, seed),
        _scene_slice(rng, seed),
    )


def scenes_for_phase(
    phase: int,
    *,
    seed: int | None = None,
    fixtures: bool = False,
    n: int | None = None,
) -> tuple[Scene, ...]:
    if phase == 1:
        if fixtures:
            return PHASE1_SCENES
        return make_phase1_scenes(choose_seed(seed), n=n)
    if phase == 2:
        from vd3d.viz.scenes_arrangement import PHASE2_SCENES, make_phase2_scenes

        if fixtures:
            return PHASE2_SCENES
        return make_phase2_scenes(choose_seed(seed), n=n)
    if phase == 3:
        from vd3d.viz.scenes_vd import PHASE3_SCENES, make_phase3_scenes

        if fixtures:
            return PHASE3_SCENES
        return make_phase3_scenes(choose_seed(seed), n=n)
    if phase == 4:
        from vd3d.viz.scenes_zone import PHASE4_SCENES, make_phase4_scenes

        if fixtures:
            return PHASE4_SCENES
        return make_phase4_scenes(choose_seed(seed), n=n)
    if phase == 5:
        from vd3d.viz.scenes_events import PHASE5_SCENES, make_phase5_scenes

        if fixtures:
            return PHASE5_SCENES
        return make_phase5_scenes(choose_seed(seed), n=n)
    if phase == 6:
        from vd3d.viz.scenes_alignment import PHASE6_SCENES, make_phase6_scenes

        if fixtures:
            return PHASE6_SCENES
        return make_phase6_scenes(choose_seed(seed), n=n)
    if phase == 7:
        from vd3d.viz.scenes_event_list import PHASE7_SCENES, make_phase7_scenes

        if fixtures:
            return PHASE7_SCENES
        return make_phase7_scenes(choose_seed(seed), n=n)
    if phase == 8:
        from vd3d.viz.scenes_sweep import PHASE8_SCENES, make_phase8_scenes

        if fixtures:
            return PHASE8_SCENES
        return make_phase8_scenes(choose_seed(seed), n=n)
    raise ValueError(f"no interactive scenes registered for phase {phase}")
