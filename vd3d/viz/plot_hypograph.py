"""Draw hypograph VD cells (floor unbounded, ceiling on the envelope)."""

from __future__ import annotations

from fractions import Fraction

import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from vd3d.cells3d.types import Cell3D
from vd3d.lower_envelope.vd import HypographVD, Trapezoid
from vd3d.lower_envelope.vd import trapezoid_representative as exact_trap_rep
from vd3d.lower_envelope.vd_float import FloatHypographVD, FloatTrapezoid
from vd3d.lower_envelope.vd_float import trapezoid_representative as float_trap_rep
from vd3d.viz.plot_cells3d import (
    GEOM_EPS,
    ViewBox,
    cell_halfspaces,
    draw_bbox_wire,
    draw_planes_grey,
    faces_from_halfspaces,
    halfspace_vertices,
)
from vd3d.viz.plot_lower_envelope import draw_envelope_surface, envelope_view_half

AnyVD = FloatHypographVD | HypographVD
AnyTrap = FloatTrapezoid | Trapezoid


def hypograph_view_half(vd: AnyVD, *, minimum: float = 3.0) -> float:
    return envelope_view_half(vd.envelope, minimum=minimum)


def _dual_row(vd: AnyVD, plane_index: int) -> tuple[float, float, float]:
    row = vd.envelope.dual_points[plane_index]
    if isinstance(vd, FloatHypographVD):
        return float(row[0]), float(row[1]), float(row[2])
    return float(row[0]), float(row[1]), float(row[2])


def _xy_of(trap: AnyTrap, prefer_half: float) -> tuple[float, float]:
    if isinstance(trap, FloatTrapezoid):
        return float_trap_rep(trap, prefer_half=prefer_half)
    p = exact_trap_rep(
        trap, prefer_half=Fraction(prefer_half).limit_denominator(1000)
    )
    return float(p.x), float(p.y)


def interior_point_hypograph(
    vd: AnyVD, trap: AnyTrap, *, prefer_half: float = 3.0
) -> np.ndarray:
    x, y = _xy_of(trap, prefer_half)
    a, b, c = _dual_row(vd, trap.plane_index)
    z_ceil = a * x + b * y + c
    z = min(z_ceil - 0.25, prefer_half - 0.25)
    return np.array([x, y, z], dtype=float)


def clipped_hypograph_faces(
    vd: AnyVD, cell: Cell3D, trap: AnyTrap, box: ViewBox
) -> list[np.ndarray]:
    interior = interior_point_hypograph(vd, trap, prefer_half=box.half)
    hs_cell = cell_halfspaces(cell, interior)
    if not all(float(np.dot(n, interior)) <= d + 10 * GEOM_EPS for n, d in box.halfspaces()):
        interior = interior_point_hypograph(vd, trap, prefer_half=min(box.half, 3.0))
        hs_cell = cell_halfspaces(cell, interior)
    hs = hs_cell + box.halfspaces()
    verts = halfspace_vertices(hs)
    if len(verts) < 3:
        return []
    return faces_from_halfspaces(hs, verts)


def draw_hypograph_cell(
    ax,
    vd: AnyVD,
    cell: Cell3D,
    trap: AnyTrap,
    box: ViewBox,
    *,
    facecolor: tuple[float, float, float, float],
    edgecolor: tuple[float, float, float, float],
    linewidth: float = 0.6,
    mark: bool = False,
) -> None:
    faces = clipped_hypograph_faces(vd, cell, trap, box)
    if not faces:
        return
    ax.add_collection3d(
        Poly3DCollection(
            faces,
            facecolors=facecolor,
            edgecolors=edgecolor,
            linewidths=linewidth,
        )
    )
    if mark:
        p = interior_point_hypograph(vd, trap, prefer_half=box.half)
        if np.all(np.abs(p) <= box.half + GEOM_EPS):
            ax.scatter([p[0]], [p[1]], [p[2]], color=(0.7, 0.0, 0.0), s=28)


def draw_hypograph_scene(ax, vd: AnyVD, index: int, box: ViewBox) -> None:
    """Blue envelope, all hypograph cells faintly, selected cell in strong red."""
    draw_planes_grey(ax, vd.envelope.planes, box)
    draw_envelope_surface(ax, vd.envelope, box)
    for i, (cell, trap) in enumerate(zip(vd.cells, vd.trapezoids, strict=True)):
        if i == index:
            continue
        draw_hypograph_cell(
            ax,
            vd,
            cell,
            trap,
            box,
            facecolor=(0.85, 0.35, 0.3, 0.12),
            edgecolor=(0.55, 0.2, 0.15, 0.25),
            linewidth=0.4,
        )
    cell = vd.cells[index]
    trap = vd.trapezoids[index]
    draw_hypograph_cell(
        ax,
        vd,
        cell,
        trap,
        box,
        facecolor=(0.85, 0.15, 0.12, 0.50),
        edgecolor=(0.55, 0.05, 0.05, 0.95),
        linewidth=1.1,
        mark=True,
    )
    draw_bbox_wire(ax, box)
