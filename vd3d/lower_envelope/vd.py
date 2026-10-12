"""Exact vertical decomposition of the hypograph under the lower envelope."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence
from dataclasses import dataclass

from vd3d.cells3d.types import Cell3D
from vd3d.geometry.plane import Plane
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar, as_scalar
from vd3d.io import format_plane
from vd3d.lower_envelope.envelope import LowerEnvelope, lower_envelope

# Sentinel half-width used only to pick finite sample points in unbounded cells.
_SAMPLE = as_scalar(10)


@dataclass(frozen=True, slots=True)
class DiagramEdge:
    """An edge of the minimization diagram in xy (segment, ray, or line)."""

    i: int
    j: int
    start: Point2D | None
    end: Point2D | None
    direction: tuple[Scalar, Scalar]


@dataclass(frozen=True, slots=True)
class Trapezoid:
    """One cell of the 2D trapezoidal map of the minimization diagram."""

    id: int
    plane_index: int
    left_x: Scalar | None
    right_x: Scalar | None
    bottom: tuple[Scalar, Scalar] | None  # y = m x + b
    top: tuple[Scalar, Scalar] | None


@dataclass(frozen=True, slots=True)
class HypographVD:
    """Exact VD of cells unbounded from below (hypograph of the lower envelope)."""

    envelope: LowerEnvelope
    trapezoids: tuple[Trapezoid, ...]
    cells: tuple[Cell3D, ...]
    edges: tuple[DiagramEdge, ...]


@dataclass(frozen=True, slots=True)
class HypographCompare:
    """Float vs exact hypograph VD comparison."""

    cell_count_match: bool
    ceiling_multiset_match: bool
    matched_cells: int
    max_x_error: float | None
    messages: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return (
            self.cell_count_match
            and self.ceiling_multiset_match
            and (self.max_x_error is None or self.max_x_error <= 1e-7)
            and not any(m.startswith("unmatched") for m in self.messages)
        )


def _line_dir_ij(
    dual: Sequence[tuple[Scalar, Scalar, Scalar]], i: int, j: int
) -> tuple[Scalar, Scalar]:
    da = dual[i][0] - dual[j][0]
    db = dual[i][1] - dual[j][1]
    return (-db, da)


def _y_on_line_ij(
    dual: Sequence[tuple[Scalar, Scalar, Scalar]], i: int, j: int, x: Scalar
) -> Scalar | None:
    da = dual[i][0] - dual[j][0]
    db = dual[i][1] - dual[j][1]
    dc = dual[i][2] - dual[j][2]
    if db == 0:
        return None
    return -(da * x + dc) / db


def _mb_for_pair(
    dual: Sequence[tuple[Scalar, Scalar, Scalar]], i: int, j: int
) -> tuple[Scalar, Scalar] | None:
    da = dual[i][0] - dual[j][0]
    db = dual[i][1] - dual[j][1]
    dc = dual[i][2] - dual[j][2]
    if db == 0:
        return None
    return (-da / db, -dc / db)


def _eval(dual_row: tuple[Scalar, Scalar, Scalar], x: Scalar, y: Scalar) -> Scalar:
    return dual_row[0] * x + dual_row[1] * y + dual_row[2]


def _jointly_minimal(
    dual: Sequence[tuple[Scalar, Scalar, Scalar]],
    i: int,
    j: int,
    x: Scalar,
    y: Scalar,
) -> bool:
    z = _eval(dual[i], x, y)
    if _eval(dual[j], x, y) != z:
        return False
    return all(_eval(row, x, y) >= z for row in dual)


def _argmin_plane(
    dual: Sequence[tuple[Scalar, Scalar, Scalar]], x: Scalar, y: Scalar
) -> int:
    best = 0
    best_z = _eval(dual[0], x, y)
    for k in range(1, len(dual)):
        z = _eval(dual[k], x, y)
        if z < best_z:
            best = k
            best_z = z
    return best


def _dot2(
    a: tuple[Scalar, Scalar], b: tuple[Scalar, Scalar]
) -> Scalar:
    return a[0] * b[0] + a[1] * b[1]


def build_diagram_edges(env: LowerEnvelope) -> list[DiagramEdge]:
    dual = env.dual_points
    edge_to_facets: dict[tuple[int, int], list[int]] = defaultdict(list)
    for fi, facet in enumerate(env.facets):
        a, b, c = facet
        for u, v in ((a, b), (b, c), (c, a)):
            key = (u, v) if u <= v else (v, u)
            edge_to_facets[key].append(fi)

    edges: list[DiagramEdge] = []
    for (i, j), facet_ids in sorted(edge_to_facets.items()):
        direction = _line_dir_ij(dual, i, j)
        points = [
            Point2D(env.vertices[fi].x, env.vertices[fi].y) for fi in facet_ids
        ]
        if points:
            origin = (points[0].x, points[0].y)
            keyed = sorted(
                (
                    (
                        _dot2(
                            (p.x - origin[0], p.y - origin[1]),
                            direction,
                        ),
                        p,
                    )
                    for p in points
                ),
                key=lambda t: t[0],
            )
        else:
            keyed = []
        uniq: list[Point2D] = []
        for _, p in keyed:
            if not uniq or (p.x != uniq[-1].x or p.y != uniq[-1].y):
                uniq.append(p)

        if len(facet_ids) >= 2 and len(uniq) >= 2:
            edges.append(
                DiagramEdge(
                    i=i, j=j, start=uniq[0], end=uniq[-1], direction=direction
                )
            )
        elif len(uniq) == 1:
            p = uniq[0]
            ray_dir = None
            step = as_scalar("1/1000")
            for sign in (1, -1):
                dx, dy = direction[0] * sign, direction[1] * sign
                qx = p.x + dx * step
                qy = p.y + dy * step
                if _jointly_minimal(dual, i, j, qx, qy):
                    ray_dir = (dx, dy)
                    break
            if ray_dir is None:
                step = as_scalar("1/10")
                for sign in (1, -1):
                    dx, dy = direction[0] * sign, direction[1] * sign
                    qx = p.x + dx * step
                    qy = p.y + dy * step
                    if _jointly_minimal(dual, i, j, qx, qy):
                        ray_dir = (dx, dy)
                        break
            if ray_dir is None:
                continue
            edges.append(
                DiagramEdge(i=i, j=j, start=p, end=None, direction=ray_dir)
            )

    if not edges and len(env.planes) == 2 and env.edges:
        i, j = env.edges[0]
        edges.append(
            DiagramEdge(
                i=i,
                j=j,
                start=None,
                end=None,
                direction=_line_dir_ij(dual, i, j),
            )
        )
    return edges


def _point_on_edge_at_x(edge: DiagramEdge, x: Scalar) -> Scalar | None:
    if edge.start is not None and edge.end is not None:
        x0, y0 = edge.start.x, edge.start.y
        x1, y1 = edge.end.x, edge.end.y
        lo, hi = (x0, x1) if x0 <= x1 else (x1, x0)
        if x < lo or x > hi:
            return None
        if x1 == x0:
            return (y0 + y1) / 2
        t = (x - x0) / (x1 - x0)
        if t < 0 or t > 1:
            return None
        return y0 + t * (y1 - y0)

    d = edge.direction
    if d[0] == 0:
        origin = edge.start if edge.start is not None else edge.end
        if origin is None:
            return None
        if origin.x != x:
            return None
        return origin.y

    if edge.start is None and edge.end is None:
        return None

    origin = edge.start if edge.start is not None else edge.end
    assert origin is not None
    t = (x - origin.x) / d[0]
    if edge.start is not None and edge.end is None and t < 0:
        return None
    if edge.start is None and edge.end is not None and t > 0:
        return None
    return origin.y + t * d[1]


def _edge_y_at_x(
    edge: DiagramEdge,
    dual: Sequence[tuple[Scalar, Scalar, Scalar]],
    x: Scalar,
) -> Scalar | None:
    if edge.start is None and edge.end is None:
        return _y_on_line_ij(dual, edge.i, edge.j, x)
    y = _point_on_edge_at_x(edge, x)
    if y is None:
        return None
    if not _jointly_minimal(dual, edge.i, edge.j, x, y):
        return None
    return y


def bands_at_x(
    dual: Sequence[tuple[Scalar, Scalar, Scalar]], x: Scalar
) -> list[tuple[Scalar | None, Scalar | None, int]]:
    """1D lower-envelope bands at fixed ``x``: ``(y_lo|None, y_hi|None, plane)``."""
    n = len(dual)
    bcoef = [row[1] for row in dual]
    dcoef = [row[0] * x + row[2] for row in dual]
    ys: list[Scalar] = []
    for i in range(n):
        for j in range(i + 1, n):
            if bcoef[i] == bcoef[j]:
                continue
            y = (dcoef[j] - dcoef[i]) / (bcoef[i] - bcoef[j])
            z = bcoef[i] * y + dcoef[i]
            if all(bcoef[k] * y + dcoef[k] >= z for k in range(n)):
                ys.append(y)
    ys = sorted(set(ys))
    bounds: list[Scalar | None] = [None, *ys, None]
    bands: list[tuple[Scalar | None, Scalar | None, int]] = []
    for a, b in zip(bounds, bounds[1:]):
        if a is not None and b is not None and b <= a:
            continue
        if a is None and b is None:
            ymid: Scalar = as_scalar(0)
        elif a is None:
            assert b is not None
            ymid = b - 1
        elif b is None:
            ymid = a + 1
        else:
            ymid = (a + b) / 2
        plane = _argmin_plane(dual, x, ymid)
        # Verify the open band is non-empty for this plane (skip degenerates).
        bands.append((a, b, plane))
    return bands


def _wall_spans_at_vertices(
    env: LowerEnvelope, edges: Sequence[DiagramEdge]
) -> list[tuple[Scalar, Scalar | None, Scalar | None]]:
    dual = env.dual_points
    spans: list[tuple[Scalar, Scalar | None, Scalar | None]] = []
    for vertex in env.vertices:
        vx, vy = vertex.x, vertex.y
        hits: list[Scalar] = [vy]
        for edge in edges:
            y = _edge_y_at_x(edge, dual, vx)
            if y is not None and y != vy:
                hits.append(y)
        for y_lo, y_hi, _p in bands_at_x(dual, vx):
            if y_lo is not None:
                hits.append(y_lo)
            if y_hi is not None:
                hits.append(y_hi)
        hits = sorted(set(hits))
        below = [y for y in hits if y < vy]
        above = [y for y in hits if y > vy]
        spans.append((vx, below[-1] if below else None, above[0] if above else None))
    return spans


def _intervals_overlap(
    a0: Scalar | None,
    a1: Scalar | None,
    b0: Scalar | None,
    b1: Scalar | None,
) -> bool:
    # Treat None as ±∞.
    lo_candidates = [v for v in (a0, b0) if v is not None]
    hi_candidates = [v for v in (a1, b1) if v is not None]
    lo = max(lo_candidates) if lo_candidates else None
    hi = min(hi_candidates) if hi_candidates else None
    if lo is None or hi is None:
        return True
    return hi > lo


def _finite_mid_x(xL: Scalar | None, xR: Scalar | None) -> Scalar:
    if xL is None and xR is None:
        return as_scalar(0)
    if xL is None:
        assert xR is not None
        return xR - 1
    if xR is None:
        return xL + 1
    return (xL + xR) / 2


def _mb_from_break(
    dual: Sequence[tuple[Scalar, Scalar, Scalar]],
    x: Scalar,
    y: Scalar,
    plane: int,
) -> tuple[Scalar, Scalar] | None:
    z = _eval(dual[plane], x, y)
    partners = [
        k
        for k in range(len(dual))
        if k != plane and _eval(dual[k], x, y) == z
    ]
    if not partners:
        return None
    return _mb_for_pair(dual, plane, partners[0])


def _boundary_mb(
    dual: Sequence[tuple[Scalar, Scalar, Scalar]],
    edges: Sequence[DiagramEdge],
    x: Scalar,
    y: Scalar,
    plane: int,
) -> tuple[Scalar, Scalar] | None:
    for edge in edges:
        if plane not in (edge.i, edge.j):
            continue
        ey = _edge_y_at_x(edge, dual, x)
        if ey is not None and ey == y:
            return _mb_for_pair(dual, edge.i, edge.j)
    return _mb_from_break(dual, x, y, plane)


def build_trapezoids(env: LowerEnvelope) -> tuple[list[Trapezoid], list[DiagramEdge]]:
    dual = env.dual_points
    edges = build_diagram_edges(env)
    xs = sorted({v.x for v in env.vertices})
    bounds: list[Scalar | None] = [None, *xs, None]
    wall_spans = _wall_spans_at_vertices(env, edges)

    pieces: list[tuple[Scalar | None, Scalar | None, Scalar | None, Scalar | None, int]] = []
    for xL, xR in zip(bounds, bounds[1:]):
        if xL is not None and xR is not None and xR <= xL:
            continue
        xmid = _finite_mid_x(xL, xR)
        for y_lo, y_hi, plane in bands_at_x(dual, xmid):
            pieces.append((xL, xR, y_lo, y_hi, plane))

    merged: list[tuple[Scalar | None, Scalar | None, Scalar | None, Scalar | None, int]] = []
    used = [False] * len(pieces)

    def blocked(x: Scalar | None, y_lo: Scalar | None, y_hi: Scalar | None) -> bool:
        if x is None:
            return False
        for wx, w_lo, w_hi in wall_spans:
            if wx != x:
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
                if oplane != plane or oy_lo != y_lo or oy_hi != y_hi:
                    continue
                if oxL == xR and not blocked(xR, y_lo, y_hi):
                    xR = oxR
                    used[j] = True
                    changed = True
                elif oxR == xL and not blocked(xL, y_lo, y_hi):
                    xL = oxL
                    used[j] = True
                    changed = True
        merged.append((xL, xR, y_lo, y_hi, plane))

    traps: list[Trapezoid] = []
    for tid, (xL, xR, y_lo, y_hi, plane) in enumerate(merged):
        xmid = _finite_mid_x(xL, xR)
        bottom = None if y_lo is None else _boundary_mb(dual, edges, xmid, y_lo, plane)
        top = None if y_hi is None else _boundary_mb(dual, edges, xmid, y_hi, plane)
        traps.append(
            Trapezoid(
                id=tid,
                plane_index=plane,
                left_x=xL,
                right_x=xR,
                bottom=bottom,
                top=top,
            )
        )
    return traps, edges


def _wall_plane_x(x: Scalar, wall_id: int) -> Plane:
    return Plane(id=wall_id, a=1, b=0, c=0, d=-x)


def _wall_plane_y_line(m: Scalar, b: Scalar, wall_id: int) -> Plane:
    return Plane(id=wall_id, a=-m, b=1, c=0, d=-b)


def trapezoid_to_cell3d(
    trap: Trapezoid, planes: Sequence[Plane], *, cell_id: int
) -> Cell3D:
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


def hypograph_vd(planes: Sequence[Plane] | LowerEnvelope) -> HypographVD:
    """Compute the exact hypograph VD (cells with ``floor=None``)."""
    env = planes if isinstance(planes, LowerEnvelope) else lower_envelope(planes)
    traps, edges = build_trapezoids(env)
    cells = tuple(
        trapezoid_to_cell3d(trap, env.planes, cell_id=trap.id) for trap in traps
    )
    return HypographVD(
        envelope=env,
        trapezoids=tuple(traps),
        cells=cells,
        edges=tuple(edges),
    )


def trapezoid_representative(
    trap: Trapezoid, *, prefer_half: Scalar | None = None
) -> Point2D:
    half = _SAMPLE if prefer_half is None else prefer_half
    xL = -half * 100 if trap.left_x is None else trap.left_x
    xR = half * 100 if trap.right_x is None else trap.right_x
    lo = max(xL, -half)
    hi = min(xR, half)
    if hi > lo:
        x = (lo + hi) / 2
    else:
        x = _finite_mid_x(trap.left_x, trap.right_x)

    def y_at(mb: tuple[Scalar, Scalar] | None, default: Scalar) -> Scalar:
        if mb is None:
            return default
        m, b = mb
        return m * x + b

    y_lo = y_at(trap.bottom, -half * 100)
    y_hi = y_at(trap.top, half * 100)
    lo_y = max(y_lo, -half)
    hi_y = min(y_hi, half)
    if hi_y > lo_y:
        y = (lo_y + hi_y) / 2
    else:
        y = (y_lo + y_hi) / 2
    return Point2D(x, y)


def format_hypograph_vd(vd: HypographVD, *, seed: int | None = None) -> str:
    lines = [
        f"planes: {len(vd.envelope.planes)}",
        f"active: {len(vd.envelope.active)}",
        f"diagram_edges: {len(vd.edges)}",
        f"trapezoids: {len(vd.trapezoids)}",
        f"cells: {len(vd.cells)}",
        "arithmetic: exact Fraction",
    ]
    if seed is not None:
        lines.append(f"seed: {seed}")
    lines.append("")
    lines.append("hypograph cells (floor=unbounded):")
    for trap, cell in zip(vd.trapezoids, vd.cells, strict=True):
        ceil = cell.ceiling
        assert ceil is not None
        lx = "-∞" if trap.left_x is None else str(trap.left_x)
        rx = "+∞" if trap.right_x is None else str(trap.right_x)
        lines.append(
            f"  cell {cell.id}: ceiling=P{ceil.id}:{format_plane(ceil)} "
            f"plane_index={trap.plane_index} x=({lx}, {rx}) "
            f"walls={len(cell.vertical_walls)}"
        )
    return "\n".join(lines) + "\n"


def compare_hypograph_float_exact(float_vd, exact: HypographVD) -> HypographCompare:
    """Compare a :class:`FloatHypographVD` to an exact :class:`HypographVD`."""
    messages: list[str] = []
    cell_count_match = len(float_vd.cells) == len(exact.cells)
    if not cell_count_match:
        messages.append(
            f"cell count: float={len(float_vd.cells)} exact={len(exact.cells)}"
        )

    def ceiling_id(cell: Cell3D) -> int:
        assert cell.ceiling is not None
        return cell.ceiling.id

    from collections import Counter

    ceil_f = Counter(ceiling_id(c) for c in float_vd.cells)
    ceil_e = Counter(ceiling_id(c) for c in exact.cells)
    ceiling_multiset_match = ceil_f == ceil_e
    if not ceiling_multiset_match:
        messages.append(f"ceiling multiset: float={dict(ceil_f)} exact={dict(ceil_e)}")

    # Match cells by (ceiling_id, left_x≈, right_x≈).
    exact_used = [False] * len(exact.trapezoids)
    matched = 0
    errors: list[float] = []
    for ft in float_vd.trapezoids:
        fceil = float_vd.envelope.planes[ft.plane_index].id
        best_j = None
        best_err = None
        for j, et in enumerate(exact.trapezoids):
            if exact_used[j]:
                continue
            eceil = exact.envelope.planes[et.plane_index].id
            if eceil != fceil:
                continue
            err = 0.0
            for fa, ea in (
                (ft.left_x, et.left_x),
                (ft.right_x, et.right_x),
            ):
                if (fa is None) != (ea is None):
                    err = 1e9
                    break
                if fa is not None and ea is not None:
                    err = max(err, abs(float(fa) - float(ea)))
            if best_err is None or err < best_err:
                best_err = err
                best_j = j
        if best_j is not None and best_err is not None and best_err < 1e-5:
            exact_used[best_j] = True
            matched += 1
            errors.append(best_err)
        else:
            messages.append(
                f"unmatched float cell plane_index={ft.plane_index} "
                f"x=({ft.left_x}, {ft.right_x})"
            )

    max_x_error = max(errors) if errors else (0.0 if matched == 0 and cell_count_match else None)
    if cell_count_match and matched == len(exact.cells) and not messages:
        messages.append("float and exact hypograph VDs match")
    return HypographCompare(
        cell_count_match=cell_count_match,
        ceiling_multiset_match=ceiling_multiset_match,
        matched_cells=matched,
        max_x_error=max_x_error,
        messages=tuple(messages),
    )


def format_hypograph_compare(cmp: HypographCompare) -> str:
    status = "OK" if cmp.ok else "MISMATCH"
    lines = [f"float vs exact hypograph VD: {status}"]
    lines.append(f"matched cells: {cmp.matched_cells}")
    if cmp.max_x_error is not None:
        lines.append(f"max x-boundary error: {cmp.max_x_error:.3e}")
    for msg in cmp.messages:
        lines.append(f"  {msg}")
    return "\n".join(lines) + "\n"
