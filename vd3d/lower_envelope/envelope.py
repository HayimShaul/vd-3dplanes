"""Exact lower envelope of planes via dual lower 3D convex hull."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from vd3d.geometry.intersections import intersect_three_planes
from vd3d.geometry.plane import Plane
from vd3d.geometry.points import Point3D
from vd3d.geometry.scalar import Scalar
from vd3d.io import format_plane
from vd3d.lower_envelope.dual import dual_points
from vd3d.lower_envelope.hull import edges_from_facets, lower_hull_facets


@dataclass(frozen=True, slots=True)
class LowerEnvelope:
    """Exact lower envelope of non-vertical planes."""

    planes: tuple[Plane, ...]
    dual_points: tuple[tuple[Scalar, Scalar, Scalar], ...]
    facets: tuple[tuple[int, int, int], ...]
    vertices: tuple[Point3D, ...]
    edges: tuple[tuple[int, int], ...]
    active: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class EnvelopeCompare:
    """Float vs exact lower-envelope comparison."""

    active_match: bool
    facets_match: bool
    edge_match: bool
    vertex_count_match: bool
    max_vertex_error: float | None
    messages: tuple[str, ...]

    @property
    def ok(self) -> bool:
        return (
            self.active_match
            and self.facets_match
            and self.edge_match
            and self.vertex_count_match
            and (self.max_vertex_error is None or self.max_vertex_error <= 1e-8)
        )


def lower_envelope(planes: Sequence[Plane]) -> LowerEnvelope:
    """Compute the exact lower envelope via dual lower hull."""
    packed = tuple(planes)
    dual = dual_points(packed)
    n = len(packed)

    if n == 1:
        facets: tuple[tuple[int, int, int], ...] = ()
        vertices: tuple[Point3D, ...] = ()
        edges: tuple[tuple[int, int], ...] = ()
        active = (0,)
    elif n == 2:
        facets = ()
        vertices = ()
        edges = ((0, 1),)
        active = (0, 1)
    else:
        facets = lower_hull_facets(dual)
        verts: list[Point3D] = []
        kept: list[tuple[int, int, int]] = []
        for i, j, k in facets:
            point = intersect_three_planes(packed[i], packed[j], packed[k])
            if point is None:
                continue
            kept.append((i, j, k))
            verts.append(point)
        facets = tuple(kept)
        vertices = tuple(verts)
        if n >= 3 and not facets:
            raise ValueError("lower hull has no facets (check plane orientations / inputs)")
        edges = edges_from_facets(facets) if facets else ((0, 1), (1, 2), (0, 2))
        active_set = {i for f in facets for i in f} if facets else set(range(n))
        active = tuple(sorted(active_set))

    return LowerEnvelope(
        planes=packed,
        dual_points=dual,
        facets=facets,
        vertices=vertices,
        edges=edges,
        active=active,
    )


def _facet_key(facet: tuple[int, ...]) -> frozenset[int]:
    return frozenset(facet)


def compare_envelopes_float_exact(
    float_env,
    exact: LowerEnvelope,
    *,
    vertex_tol: float = 1e-8,
) -> EnvelopeCompare:
    """Compare a :class:`FloatLowerEnvelope` to an exact :class:`LowerEnvelope`."""
    messages: list[str] = []
    active_f = frozenset(int(i) for i in float_env.active)
    active_e = frozenset(exact.active)
    active_match = active_f == active_e
    if not active_match:
        messages.append(f"active mismatch: float={sorted(active_f)} exact={sorted(active_e)}")

    facets_f = frozenset(_facet_key(tuple(int(x) for x in f)) for f in float_env.facets)
    facets_e = frozenset(_facet_key(f) for f in exact.facets)
    facets_match = facets_f == facets_e
    if not facets_match:
        messages.append(
            f"facet mismatch: float_only={facets_f - facets_e} exact_only={facets_e - facets_f}"
        )

    edges_f = frozenset(tuple(sorted(e)) for e in float_env.edges)
    edges_e = frozenset(tuple(sorted(e)) for e in exact.edges)
    edge_match = edges_f == edges_e
    if not edge_match:
        messages.append(
            f"edge mismatch: float_only={edges_f - edges_e} exact_only={edges_e - edges_f}"
        )

    vertex_count_match = len(float_env.vertices) == len(exact.vertices)
    if not vertex_count_match:
        messages.append(
            f"vertex count: float={len(float_env.vertices)} exact={len(exact.vertices)}"
        )

    max_err: float | None = None
    if vertex_count_match and exact.vertices:
        # Match vertices by facet plane-set.
        exact_by_key = {
            _facet_key(f): v for f, v in zip(exact.facets, exact.vertices, strict=True)
        }
        errs: list[float] = []
        for facet, fvert in zip(float_env.facets, float_env.vertices, strict=True):
            key = _facet_key(tuple(int(x) for x in facet))
            evertex = exact_by_key.get(key)
            if evertex is None:
                messages.append(f"no exact vertex for float facet {key}")
                continue
            dx = abs(float(fvert[0]) - float(evertex.x))
            dy = abs(float(fvert[1]) - float(evertex.y))
            dz = abs(float(fvert[2]) - float(evertex.z))
            errs.append(max(dx, dy, dz))
        max_err = max(errs) if errs else 0.0
        if max_err > vertex_tol:
            messages.append(f"max vertex error {max_err:.3e} > tol {vertex_tol:g}")
    elif not exact.vertices and not len(float_env.vertices):
        max_err = 0.0

    if not messages:
        messages.append("float and exact envelopes match")

    return EnvelopeCompare(
        active_match=active_match,
        facets_match=facets_match,
        edge_match=edge_match,
        vertex_count_match=vertex_count_match,
        max_vertex_error=max_err,
        messages=tuple(messages),
    )


def format_envelope(env: LowerEnvelope, *, seed: int | None = None) -> str:
    lines = [
        f"planes: {len(env.planes)}",
        f"active: {len(env.active)}",
        f"facets: {len(env.facets)}",
        f"edges: {len(env.edges)}",
        "arithmetic: exact Fraction",
    ]
    if seed is not None:
        lines.append(f"seed: {seed}")
    lines.extend(["", "input planes:"])
    for plane in env.planes:
        lines.append(f"  P{plane.id}: {format_plane(plane)}")
    lines.append("")
    lines.append("active plane indices:")
    for idx in env.active:
        lines.append(f"  [{idx}] P{env.planes[idx].id}")
    lines.append("")
    lines.append("envelope vertices (triple intersections):")
    if not env.vertices:
        lines.append("  (none)")
    else:
        for facet, vertex in zip(env.facets, env.vertices, strict=True):
            ids = tuple(env.planes[i].id for i in facet)
            lines.append(
                f"  planes {ids} -> ({vertex.x}, {vertex.y}, {vertex.z})"
            )
    return "\n".join(lines) + "\n"


def format_compare(cmp: EnvelopeCompare) -> str:
    status = "OK" if cmp.ok else "MISMATCH"
    lines = [f"float vs exact: {status}"]
    if cmp.max_vertex_error is not None:
        lines.append(f"max vertex error: {cmp.max_vertex_error:.3e}")
    for msg in cmp.messages:
        lines.append(f"  {msg}")
    return "\n".join(lines) + "\n"
