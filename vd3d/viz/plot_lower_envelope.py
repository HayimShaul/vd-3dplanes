"""Drawing helpers for the lower-envelope GUI."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from vd3d.geometry.plane import Plane
from vd3d.lower_envelope.envelope import LowerEnvelope
from vd3d.lower_envelope.envelope_float import FloatLowerEnvelope
from vd3d.viz.convert import to_float
from vd3d.viz.plot_cells3d import ViewBox, draw_bbox_wire, plane_polygon_in_box

GEOM_EPS = 1e-8
ENVELOPE_BLUE = (0.15, 0.35, 0.85)
ENVELOPE_BLUE_EDGE = (0.05, 0.15, 0.55)
FLOAT_ORANGE = (0.95, 0.45, 0.05)


def _dual_as_float_rows(env: FloatLowerEnvelope | LowerEnvelope) -> np.ndarray:
    if isinstance(env, FloatLowerEnvelope):
        return np.asarray(env.dual_points, dtype=float)
    return np.asarray(
        [[float(a), float(b), float(c)] for a, b, c in env.dual_points], dtype=float
    )


def _planes_of(env: FloatLowerEnvelope | LowerEnvelope) -> Sequence[Plane]:
    return env.planes


def envelope_view_half(
    env: FloatLowerEnvelope | LowerEnvelope, *, minimum: float = 3.0
) -> float:
    """Cube half-width covering envelope vertices with margin."""
    values = [minimum]
    if isinstance(env, FloatLowerEnvelope):
        for vertex in env.vertices:
            values.append(float(np.max(np.abs(vertex))) + 1.0)
    else:
        for vertex in env.vertices:
            values.append(
                max(abs(to_float(vertex.x)), abs(to_float(vertex.y)), abs(to_float(vertex.z)))
                + 1.0
            )
    return float(max(values))


def draw_input_planes_grey(
    ax, env: FloatLowerEnvelope | LowerEnvelope, box: ViewBox
) -> None:
    """Faint grey input planes for context (all of them)."""
    for plane in _planes_of(env):
        poly = plane_polygon_in_box(plane, box)
        if poly is None:
            continue
        ax.add_collection3d(
            Poly3DCollection(
                [poly],
                facecolors=(0.55, 0.55, 0.55, 0.10),
                edgecolors=(0.4, 0.4, 0.4, 0.22),
                linewidths=0.35,
            )
        )


def envelope_z_on_grid(
    env: FloatLowerEnvelope | LowerEnvelope, box: ViewBox, *, n: int = 48
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sample ``z = min_i (A_i x + B_i y + C_i)`` over the xy footprint of ``box``."""
    h = box.half
    xs = np.linspace(-h, h, n)
    ys = np.linspace(-h, h, n)
    x, y = np.meshgrid(xs, ys)
    z = np.full(x.shape, np.inf, dtype=float)
    for a, b, c in _dual_as_float_rows(env):
        z = np.minimum(z, a * x + b * y + c)
    return x, y, z


def draw_envelope_surface(
    ax, env: FloatLowerEnvelope | LowerEnvelope, box: ViewBox, *, n: int = 48
) -> None:
    """Draw the lower envelope as a transparent blue piecewise-planar surface."""
    x, y, z = envelope_z_on_grid(env, box, n=n)
    z_plot = np.where(np.abs(z) <= box.half + GEOM_EPS, z, np.nan)
    ax.plot_surface(
        x,
        y,
        z_plot,
        color=(*ENVELOPE_BLUE, 0.35),
        edgecolor=(*ENVELOPE_BLUE_EDGE, 0.15),
        linewidth=0.1,
        antialiased=True,
        shade=True,
        rcount=n,
        ccount=n,
    )


def draw_envelope_vertices(ax, env: FloatLowerEnvelope | LowerEnvelope) -> None:
    if isinstance(env, FloatLowerEnvelope):
        if len(env.vertices) == 0:
            return
        xs, ys, zs = env.vertices[:, 0], env.vertices[:, 1], env.vertices[:, 2]
    else:
        if not env.vertices:
            return
        xs = [to_float(v.x) for v in env.vertices]
        ys = [to_float(v.y) for v in env.vertices]
        zs = [to_float(v.z) for v in env.vertices]
    ax.scatter(
        xs, ys, zs, color=ENVELOPE_BLUE_EDGE, s=40, depthshade=True, zorder=5, marker="o"
    )


def draw_float_vertices_overlay(ax, float_env: FloatLowerEnvelope) -> None:
    """Orange markers for float envelope vertices (compare mode)."""
    if len(float_env.vertices) == 0:
        return
    ax.scatter(
        float_env.vertices[:, 0],
        float_env.vertices[:, 1],
        float_env.vertices[:, 2],
        color=FLOAT_ORANGE,
        s=70,
        depthshade=True,
        zorder=6,
        marker="x",
        linewidths=1.5,
    )


def draw_lower_envelope_scene(
    ax,
    env: FloatLowerEnvelope | LowerEnvelope,
    box: ViewBox,
    *,
    float_overlay: FloatLowerEnvelope | None = None,
) -> None:
    draw_input_planes_grey(ax, env, box)
    draw_envelope_surface(ax, env, box)
    draw_envelope_vertices(ax, env)
    if float_overlay is not None:
        draw_float_vertices_overlay(ax, float_overlay)
    draw_bbox_wire(ax, box)
