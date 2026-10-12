"""Float vertical decomposition of the hypograph under the lower envelope.

Builds the xy minimization diagram from a :class:`FloatLowerEnvelope`, computes
its y-parallel trapezoidal map, and lifts each trapezoid to a ``Cell3D`` with
``floor=None``, ``ceiling`` = the envelope plane of that face, and vertical
walls ``x = const``.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from math import inf

import numpy as np

from vd3d.cells3d.types import Cell3D
from vd3d.geometry.plane import Plane
from vd3d.lower_envelope.envelope_float import FloatLowerEnvelope, lower_envelope_float

EPS = 1e-9
FIN = 1e6  # stand-in for ±∞ when sampling midpoints


@dataclass(frozen=True, slots=True)
class FloatDiagramEdge:
    """An edge of the minimization diagram in xy (segment or ray/line)."""

    i: int
    j: int
    # Endpoints in xy; None means unbounded in that direction along ``direction``.
    start: np.ndarray | None  # shape (2,)
    end: np.ndarray | None
    direction: np.ndarray  # unit direction from start toward end / infinity


@dataclass(frozen=True, slots=True)
class FloatTrapezoid:
    """One cell of the 2D trapezoidal map of the minimization diagram."""

    id: int
    plane_index: int
    left_x: float | None
    right_x: float | None
    # Bottom / top: y = m*x + b, or None if unbounded in that y-direction.
    bottom: tuple[float, float] | None
    top: tuple[float, float] | None


@dataclass(frozen=True, slots=True)
class FloatHypographVD:
    """Float VD of cells unbounded from below (hypograph of the lower envelope)."""

    envelope: FloatLowerEnvelope
    trapezoids: tuple[FloatTrapezoid, ...]
    cells: tuple[Cell3D, ...]
    edges: tuple[FloatDiagramEdge, ...]


def _unit(v: np.ndarray) -> np.ndarray:
    n = float(np.linalg.norm(v))
    if n < EPS:
        raise ValueError("zero direction")
    return v / n


def _line_dir_ij(dual: np.ndarray, i: int, j: int) -> np.ndarray:
    """Direction of the xy-projection of the intersection of graph-planes i, j."""
    da = dual[i, 0] - dual[j, 0]
    db = dual[i, 1] - dual[j, 1]
    # Normal to the projected line is (da, db); a direction is (-db, da).
    return _unit(np.array([-db, da], dtype=float))


def _y_on_line_ij(dual: np.ndarray, i: int, j: int, x: float) -> float | None:
    """y on the projected intersection of i,j at given x, or None if vertical."""
    da = dual[i, 0] - dual[j, 0]
    db = dual[i, 1] - dual[j, 1]
    dc = dual[i, 2] - dual[j, 2]
    # da*x + db*y + dc = 0
    if abs(db) < EPS:
        return None
    return -(da * x + dc) / db


def _mb_for_pair(dual: np.ndarray, i: int, j: int) -> tuple[float, float] | None:
    """``y = m x + b`` for the projected intersection of i,j; None if vertical."""
    da = dual[i, 0] - dual[j, 0]
    db = dual[i, 1] - dual[j, 1]
    dc = dual[i, 2] - dual[j, 2]
    if abs(db) < EPS:
        return None
    m = -da / db
    b = -dc / db
    return m, b


def _jointly_minimal(
    dual: np.ndarray, i: int, j: int, x: float, y: float, *, tol: float = 1e-7
) -> bool:
    z = dual[i, 0] * x + dual[i, 1] * y + dual[i, 2]
    zj = dual[j, 0] * x + dual[j, 1] * y + dual[j, 2]
    if abs(z - zj) > tol:
        return False
    vals = dual[:, 0] * x + dual[:, 1] * y + dual[:, 2]
    return float(np.min(vals)) >= z - tol


def build_diagram_edges(env: FloatLowerEnvelope) -> list[FloatDiagramEdge]:
    """Minimization-diagram edges from lower-hull edge incidences."""
    dual = env.dual_points
    edge_to_facets: dict[tuple[int, int], list[int]] = defaultdict(list)
    for fi, facet in enumerate(env.facets):
        a, b, c = int(facet[0]), int(facet[1]), int(facet[2])
        for u, v in ((a, b), (b, c), (c, a)):
            key = (u, v) if u <= v else (v, u)
            edge_to_facets[key].append(fi)

    edges: list[FloatDiagramEdge] = []
    for (i, j), facet_ids in sorted(edge_to_facets.items()):
        direction = _line_dir_ij(dual, i, j)
        points = [env.vertices[fi][:2].astype(float) for fi in facet_ids]
        # Sort along direction.
        origins = points[0] if points else np.zeros(2)
        keyed = sorted(
            ((float(np.dot(p - origins, direction)), p) for p in points),
            key=lambda t: t[0],
        )
        uniq: list[np.ndarray] = []
        for _, p in keyed:
            if not uniq or float(np.linalg.norm(p - uniq[-1])) > EPS:
                uniq.append(p)

        if len(facet_ids) >= 2 and len(uniq) >= 2:
            edges.append(
                FloatDiagramEdge(
                    i=i,
                    j=j,
                    start=uniq[0],
                    end=uniq[-1],
                    direction=direction,
                )
            )
        elif len(uniq) == 1:
            # Unbounded ray: pick the direction where i,j stay jointly minimal.
            p = uniq[0]
            ray_dir = None
            for cand in (direction, -direction):
                q = p + 1e-3 * cand
                if _jointly_minimal(dual, i, j, float(q[0]), float(q[1])):
                    ray_dir = cand
                    break
            if ray_dir is None:
                # Fall back: probe farther / opposite.
                for cand in (direction, -direction):
                    q = p + 0.1 * cand
                    if _jointly_minimal(dual, i, j, float(q[0]), float(q[1])):
                        ray_dir = cand
                        break
            if ray_dir is None:
                continue
            edges.append(
                FloatDiagramEdge(
                    i=i, j=j, start=p, end=None, direction=ray_dir
                )
            )
        elif len(uniq) == 0 and len(env.facets) == 0 and len(env.edges) > 0:
            # n=2 case: full line (both directions unbounded).
            y0 = _y_on_line_ij(dual, i, j, 0.0)
            if y0 is None:
                # vertical line x = const
                da = dual[i, 0] - dual[j, 0]
                dc = dual[i, 2] - dual[j, 2]
                x0 = -dc / da if abs(da) > EPS else 0.0
                p = np.array([x0, 0.0], dtype=float)
                edges.append(
                    FloatDiagramEdge(
                        i=i, j=j, start=None, end=None, direction=np.array([0.0, 1.0])
                    )
                )
            else:
                edges.append(
                    FloatDiagramEdge(
                        i=i,
                        j=j,
                        start=None,
                        end=None,
                        direction=direction,
                    )
                )
    # n=2 fallback from env.edges if no facets produced edges
    if not edges and len(env.planes) == 2 and env.edges:
        i, j = env.edges[0]
        direction = _line_dir_ij(dual, i, j)
        edges.append(
            FloatDiagramEdge(i=i, j=j, start=None, end=None, direction=direction)
        )
    return edges


def _point_on_edge_at_x(edge: FloatDiagramEdge, x: float) -> float | None:
    """If the edge crosses the vertical line at ``x``, return the y; else None."""
    # Represent edge as parametric from an origin.
    if edge.start is not None and edge.end is not None:
        x0, y0 = float(edge.start[0]), float(edge.start[1])
        x1, y1 = float(edge.end[0]), float(edge.end[1])
        lo, hi = (x0, x1) if x0 <= x1 else (x1, x0)
        if x < lo - EPS or x > hi + EPS:
            return None
        if abs(x1 - x0) < EPS:
            # vertical segment
            if abs(x - x0) > EPS:
                return None
            return 0.5 * (y0 + y1)
        t = (x - x0) / (x1 - x0)
        if t < -EPS or t > 1.0 + EPS:
            return None
        return y0 + t * (y1 - y0)

    # Ray or full line: use direction.
    d = edge.direction
    if abs(float(d[0])) < EPS:
        # vertical ray/line
        origin = edge.start if edge.start is not None else np.array([x, 0.0])
        if abs(float(origin[0]) - x) > EPS:
            return None
        return float(origin[1])

    if edge.start is None and edge.end is None:
        # Full line through origin along direction — recover using plane pair.
        # Caller should use _y_on_line_ij; store i,j on edge.
        return None  # handled separately

    origin = edge.start if edge.start is not None else edge.end
    assert origin is not None
    t = (x - float(origin[0])) / float(d[0])
    if edge.start is not None and edge.end is None:
        # ray: t >= 0 along direction from start
        if t < -EPS:
            return None
    elif edge.start is None and edge.end is not None:
        # ray ending at end: approaching end from -direction, t <= 0 if origin=end
        if t > EPS:
            return None
    y = float(origin[1]) + t * float(d[1])
    return y


def _edge_y_at_x(
    edge: FloatDiagramEdge, dual: np.ndarray, x: float
) -> float | None:
    if edge.start is None and edge.end is None:
        return _y_on_line_ij(dual, edge.i, edge.j, x)
    y = _point_on_edge_at_x(edge, x)
    if y is None:
        return None
    # Confirm the crossing lies on the envelope.
    if not _jointly_minimal(dual, edge.i, edge.j, x, y, tol=1e-5):
        return None
    return y


def bands_at_x(dual: np.ndarray, x: float) -> list[tuple[float, float, int]]:
    """1D lower envelope bands at fixed ``x``: ``(y_lo, y_hi, plane_index)``."""
    n = dual.shape[0]
    bcoef = dual[:, 1]
    dcoef = dual[:, 0] * x + dual[:, 2]
    ys: list[float] = [-inf, inf]
    for i in range(n):
        for j in range(i + 1, n):
            if abs(bcoef[i] - bcoef[j]) < EPS:
                continue
            y = (dcoef[j] - dcoef[i]) / (bcoef[i] - bcoef[j])
            z = bcoef[i] * y + dcoef[i]
            vals = bcoef * y + dcoef
            if float(np.min(vals)) >= z - 1e-7:
                ys.append(float(y))
    ys = sorted(ys)
    # unique
    uniq = [ys[0]]
    for y in ys[1:]:
        if abs(y - uniq[-1]) > EPS:
            uniq.append(y)
    bands: list[tuple[float, float, int]] = []
    for a, b in zip(uniq, uniq[1:]):
        if b <= a + EPS:
            continue
        if a == -inf and b == inf:
            ymid = 0.0
        elif a == -inf:
            ymid = b - 1.0
        elif b == inf:
            ymid = a + 1.0
        else:
            ymid = 0.5 * (a + b)
        vals = bcoef * ymid + dcoef
        plane = int(np.argmin(vals))
        bands.append((a, b, plane))
    return bands


def _wall_spans_at_vertices(
    env: FloatLowerEnvelope, edges: Sequence[FloatDiagramEdge]
) -> list[tuple[float, float, float]]:
    """For each envelope vertex, the vertical wall ``(x, y_lo, y_hi)``."""
    dual = env.dual_points
    spans: list[tuple[float, float, float]] = []
    for vertex in env.vertices:
        vx, vy = float(vertex[0]), float(vertex[1])
        # Collect y-hits of diagram edges at this x (excluding the vertex itself).
        hits = [vy]
        for edge in edges:
            y = _edge_y_at_x(edge, dual, vx)
            if y is None:
                continue
            if abs(y - vy) > EPS:
                hits.append(y)
        # Also use 1D envelope breakpoints at vx.
        for y_lo, y_hi, _p in bands_at_x(dual, vx):
            if y_lo != -inf:
                hits.append(y_lo)
            if y_hi != inf:
                hits.append(y_hi)
        hits = sorted(hits)
        # Wall through the vertex extends to the nearest hit above and below.
        below = [y for y in hits if y < vy - EPS]
        above = [y for y in hits if y > vy + EPS]
        y_lo = below[-1] if below else -inf
        y_hi = above[0] if above else inf
        spans.append((vx, y_lo, y_hi))
    return spans


def _intervals_overlap(a0: float, a1: float, b0: float, b1: float) -> bool:
    lo = max(a0, b0)
    hi = min(a1, b1)
    return hi > lo + EPS


def _finite_mid_x(xL: float, xR: float) -> float:
    if xL == -inf and xR == inf:
        return 0.0
    if xL == -inf:
        return xR - 1.0
    if xR == inf:
        return xL + 1.0
    return 0.5 * (xL + xR)


def _boundary_mb(
    dual: np.ndarray,
    edges: Sequence[FloatDiagramEdge],
    x: float,
    y: float,
    plane: int,
) -> tuple[float, float] | None:
    """Supporting ``y = m x + b`` at a band boundary point, or None if ±∞."""
    if abs(y) >= FIN / 2:
        return None
    # Prefer an edge through (x,y) incident to ``plane``.
    for edge in edges:
        if plane not in (edge.i, edge.j):
            continue
        ey = _edge_y_at_x(edge, dual, x)
        if ey is not None and abs(ey - y) < 1e-5:
            return _mb_for_pair(dual, edge.i, edge.j)
    # Fallback: pair plane with argmin among others at a nudge.
    return None


def build_trapezoids(env: FloatLowerEnvelope) -> tuple[list[FloatTrapezoid], list[FloatDiagramEdge]]:
    """Trapezoidal map of the minimization diagram."""
    dual = env.dual_points
    edges = build_diagram_edges(env)
    xs = sorted({float(v[0]) for v in env.vertices})
    # Unique x
    uniq_x: list[float] = []
    for x in xs:
        if not uniq_x or abs(x - uniq_x[-1]) > EPS:
            uniq_x.append(x)
    bounds = [-inf, *uniq_x, inf]
    wall_spans = _wall_spans_at_vertices(env, edges)

    # Candidate pieces per strip.
    pieces: list[tuple[float, float, float, float, int]] = []  # xL,xR,yLo,yHi,plane
    for xL, xR in zip(bounds, bounds[1:]):
        if xR <= xL + EPS:
            continue
        xmid = _finite_mid_x(xL, xR)
        for y_lo, y_hi, plane in bands_at_x(dual, xmid):
            pieces.append((xL, xR, y_lo, y_hi, plane))

    # Merge across walls when the open band is not blocked.
    merged: list[tuple[float, float, float, float, int]] = []
    used = [False] * len(pieces)

    def blocked(x: float, y_lo: float, y_hi: float) -> bool:
        for wx, w_lo, w_hi in wall_spans:
            if abs(wx - x) > EPS:
                continue
            if _intervals_overlap(y_lo, y_hi, w_lo, w_hi):
                return True
        return False

    for i, piece in enumerate(pieces):
        if used[i]:
            continue
        xL, xR, y_lo, y_hi, plane = piece
        used[i] = True
        changed = True
        while changed:
            changed = False
            for j, other in enumerate(pieces):
                if used[j]:
                    continue
                oxL, oxR, oy_lo, oy_hi, oplane = other
                if oplane != plane:
                    continue
                if abs(oy_lo - y_lo) > 1e-5 or abs(oy_hi - y_hi) > 1e-5:
                    # allow ±inf match
                    if not (
                        (oy_lo == y_lo or (oy_lo == -inf and y_lo == -inf))
                        and (oy_hi == y_hi or (oy_hi == inf and y_hi == inf))
                    ):
                        continue
                # adjacent strips
                if abs(oxL - xR) <= EPS and not blocked(xR, y_lo, y_hi):
                    xR = oxR
                    used[j] = True
                    changed = True
                elif abs(oxR - xL) <= EPS and not blocked(xL, y_lo, y_hi):
                    xL = oxL
                    used[j] = True
                    changed = True
        merged.append((xL, xR, y_lo, y_hi, plane))

    traps: list[FloatTrapezoid] = []
    for tid, (xL, xR, y_lo, y_hi, plane) in enumerate(merged):
        xmid = _finite_mid_x(xL, xR)
        bottom = None if y_lo == -inf else _boundary_mb(dual, edges, xmid, y_lo, plane)
        top = None if y_hi == inf else _boundary_mb(dual, edges, xmid, y_hi, plane)
        # If boundary_mb failed, synthesize from the 1D partner plane at the break.
        if y_lo != -inf and bottom is None:
            bottom = _mb_from_break(dual, xmid, y_lo, plane)
        if y_hi != inf and top is None:
            top = _mb_from_break(dual, xmid, y_hi, plane)
        traps.append(
            FloatTrapezoid(
                id=tid,
                plane_index=plane,
                left_x=None if xL == -inf else xL,
                right_x=None if xR == inf else xR,
                bottom=bottom,
                top=top,
            )
        )
    return traps, edges


def _mb_from_break(
    dual: np.ndarray, x: float, y: float, plane: int
) -> tuple[float, float] | None:
    """Recover ``y = m x + b`` from the co-minimizer at a 1D breakpoint."""
    z = dual[plane, 0] * x + dual[plane, 1] * y + dual[plane, 2]
    vals = dual[:, 0] * x + dual[:, 1] * y + dual[:, 2]
    partners = [k for k in range(dual.shape[0]) if k != plane and abs(vals[k] - z) < 1e-6]
    if not partners:
        return None
    return _mb_for_pair(dual, plane, partners[0])


def _wall_plane_x(x: float, wall_id: int) -> Plane:
    """Vertical wall ``x = const`` (parallel to the yz-plane)."""
    return Plane(id=wall_id, a=1, b=0, c=0, d=Fraction.from_float(-float(x)))


def _wall_plane_y_line(m: float, b: float, wall_id: int) -> Plane:
    """Vertical wall through the xy-line ``y = m x + b`` (parallel to the z-axis)."""
    # -m x + y - b = 0
    return Plane(
        id=wall_id,
        a=Fraction.from_float(-float(m)),
        b=1,
        c=0,
        d=Fraction.from_float(-float(b)),
    )


def trapezoid_to_cell3d(
    trap: FloatTrapezoid, planes: Sequence[Plane], *, cell_id: int
) -> Cell3D:
    """Lift a trapezoid to a hypograph prism under its envelope plane.

    Besides the ``x = const`` sides, the top/bottom diagram edges become
    vertical walls (parallel to ``z``). Without those, the prism would extend
    in ``y`` into neighboring envelope faces, where this ceiling lies *above*
    the lower envelope.
    """
    ceiling = planes[trap.plane_index]
    walls: list[Plane] = []
    wid = -1
    if trap.left_x is not None:
        walls.append(_wall_plane_x(trap.left_x, wid))
        wid -= 1
    if trap.right_x is not None:
        walls.append(_wall_plane_x(trap.right_x, wid))
        wid -= 1
    if trap.bottom is not None:
        m, b = trap.bottom
        walls.append(_wall_plane_y_line(m, b, wid))
        wid -= 1
    if trap.top is not None:
        m, b = trap.top
        walls.append(_wall_plane_y_line(m, b, wid))
    return Cell3D(
        id=cell_id,
        floor=None,
        ceiling=ceiling,
        vertical_walls=tuple(walls),
        lower_z=None,
        upper_z=None,
    )


def hypograph_vd_float(
    planes: Sequence[Plane] | FloatLowerEnvelope,
) -> FloatHypographVD:
    """Compute the float hypograph VD (cells with ``floor=None``)."""
    env = planes if isinstance(planes, FloatLowerEnvelope) else lower_envelope_float(planes)
    traps, edges = build_trapezoids(env)
    cells = tuple(
        trapezoid_to_cell3d(trap, env.planes, cell_id=trap.id) for trap in traps
    )
    return FloatHypographVD(
        envelope=env,
        trapezoids=tuple(traps),
        cells=cells,
        edges=tuple(edges),
    )


def trapezoid_representative(
    trap: FloatTrapezoid, *, prefer_half: float = 3.0
) -> tuple[float, float]:
    """An interior ``(x, y)`` sample of the trapezoid.

    Prefers a point inside ``[-prefer_half, prefer_half]^2`` when the trapezoid
    meets that square, so drawing/halfspace checks stay numerically sane.
    """
    xL = -FIN if trap.left_x is None else trap.left_x
    xR = FIN if trap.right_x is None else trap.right_x
    # Overlap of (xL, xR) with [-prefer_half, prefer_half], else a finite mid.
    lo = max(xL, -prefer_half)
    hi = min(xR, prefer_half)
    if hi > lo + EPS:
        x = 0.5 * (lo + hi)
    else:
        x = _finite_mid_x(xL, xR)
        x = max(min(x, FIN), -FIN)

    def y_at(mb: tuple[float, float] | None, default: float) -> float:
        if mb is None:
            return default
        m, b = mb
        return m * x + b

    y_lo = y_at(trap.bottom, -FIN)
    y_hi = y_at(trap.top, FIN)
    lo_y = max(y_lo, -prefer_half)
    hi_y = min(y_hi, prefer_half)
    if hi_y > lo_y + EPS:
        y = 0.5 * (lo_y + hi_y)
    elif y_hi <= y_lo + EPS:
        y = y_lo + 1.0
    else:
        y = 0.5 * (y_lo + y_hi)
        y = max(min(y, FIN), -FIN)
    return x, y


def format_float_hypograph_vd(vd: FloatHypographVD, *, seed: int | None = None) -> str:
    lines = [
        f"planes: {len(vd.envelope.planes)}",
        f"active: {len(vd.envelope.active)}",
        f"diagram_edges: {len(vd.edges)}",
        f"trapezoids: {len(vd.trapezoids)}",
        f"cells: {len(vd.cells)}",
    ]
    if seed is not None:
        lines.append(f"seed: {seed}")
    lines.append("")
    lines.append("hypograph cells (floor=unbounded):")
    for trap, cell in zip(vd.trapezoids, vd.cells, strict=True):
        ceil = cell.ceiling
        assert ceil is not None
        lx = "-∞" if trap.left_x is None else f"{trap.left_x:.6g}"
        rx = "+∞" if trap.right_x is None else f"{trap.right_x:.6g}"
        walls = len(cell.vertical_walls)
        lines.append(
            f"  cell {cell.id}: ceiling=P{ceil.id} "
            f"plane_index={trap.plane_index} x=({lx}, {rx}) walls={walls}"
        )
    return "\n".join(lines) + "\n"
