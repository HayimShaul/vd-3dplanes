"""Event-local 2D VD update (design §12, §19 Step 8).

Far vertices keep their before-event ray-hit combinatorics; vertices in
the event's x-window are reshot on the after-slice arrangement. The
recomputed ``compute_vd_at_z(z+)`` remains the oracle.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import combinations

from vd3d.arrangement2d.types import Arrangement2D, Vertex
from vd3d.events.lines import compute_intersection_lines
from vd3d.events.slice import build_2d_arrangement_at_z
from vd3d.events.types import Event, EventType
from vd3d.geometry.plane import Plane
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar
from vd3d.vertical_decomposition.cells import (
    build_cells,
    finalize_decomposition_walls,
    insert_decomposition_segment,
)
from vd3d.vertical_decomposition.geom import y_at_x
from vd3d.vertical_decomposition.rays import first_hit, vertical_rays_from
from vd3d.vertical_decomposition.types import (
    UNBOUNDED,
    Hit,
    VerticalDecomposition,
    VerticalWall,
)


@dataclass(frozen=True, slots=True)
class HitRecord:
    """Combinatorial first-hit of a vertical ray, independent of ``z`` geometry."""

    unbounded: bool
    hit_plane_id: int | None
    hit_vertex_key: tuple[int, ...] | None


def vertex_plane_key(arrangement: Arrangement2D, vertex: Vertex) -> tuple[int, ...]:
    """Sorted source-plane ids of the lines through ``vertex``."""
    ids: list[int] = []
    for line_index in vertex.line_indices:
        pid = arrangement.lines[line_index].source_plane_id
        if pid is not None:
            ids.append(pid)
    return tuple(sorted(ids))


def event_vertex_plane_keys(
    event: Event, planes: Sequence[Plane]
) -> frozenset[tuple[int, ...]]:
    """Plane-pair keys of the arrangement vertices that define ``event``."""
    if event.type is EventType.TRIPLE_INTERSECTION:
        if len(event.plane_ids) < 3:
            raise ValueError("triple event needs three plane ids")
        return frozenset(
            tuple(sorted(pair)) for pair in combinations(event.plane_ids, 2)
        )
    if event.type is EventType.VERTICAL_ALIGNMENT:
        return _alignment_vertex_keys(event, planes)
    raise ValueError(f"unknown event type: {event.type!r}")


def local_vertex_plane_keys(
    vd_before: VerticalDecomposition,
    arr_after: Arrangement2D,
    event: Event,
    planes: Sequence[Plane],
) -> frozenset[tuple[int, ...]]:
    """Event vertices plus every vertex in the closed x-window they span.

    The window is the union of the event vertices' ``x`` at ``z−`` and at
    ``z+``. That is the strip whose Steiner walls collide or flip.
    """
    event_keys = event_vertex_plane_keys(event, planes)
    xs: list[Scalar] = []
    for arrangement in (vd_before.arrangement, arr_after):
        for vertex in arrangement.vertices:
            if vertex_plane_key(arrangement, vertex) in event_keys:
                xs.append(vertex.point.x)
    extra: set[tuple[int, ...]] = set(event_keys)
    if not xs:
        return frozenset(extra)
    x_lo, x_hi = min(xs), max(xs)
    for arrangement in (vd_before.arrangement, arr_after):
        for vertex in arrangement.vertices:
            x = vertex.point.x
            if x_lo <= x <= x_hi:
                extra.add(vertex_plane_key(arrangement, vertex))
    return frozenset(extra)


def update_2d_decomposition(
    vd_before: VerticalDecomposition,
    event: Event,
    planes: Sequence[Plane],
    z_after: int | Scalar | str,
) -> VerticalDecomposition:
    """``UPDATE_2D_DECOMPOSITION``: dispatch on event type (design §12)."""
    if event.type is EventType.TRIPLE_INTERSECTION:
        return update_for_triple_intersection(vd_before, event, planes, z_after)
    if event.type is EventType.VERTICAL_ALIGNMENT:
        return update_for_vertical_alignment(vd_before, event, planes, z_after)
    raise ValueError(f"unknown event type: {event.type!r}")


def update_for_triple_intersection(
    vd_before: VerticalDecomposition,
    event: Event,
    planes: Sequence[Plane],
    z_after: int | Scalar | str,
) -> VerticalDecomposition:
    """``UPDATE_FOR_TRIPLE_INTERSECTION``: reshoot the triangle's x-window."""
    if event.type is not EventType.TRIPLE_INTERSECTION:
        raise ValueError("UPDATE_FOR_TRIPLE_INTERSECTION requires a triple event")
    return _update_local(vd_before, event, planes, z_after)


def update_for_vertical_alignment(
    vd_before: VerticalDecomposition,
    event: Event,
    planes: Sequence[Plane],
    z_after: int | Scalar | str,
) -> VerticalDecomposition:
    """``UPDATE_FOR_VERTICAL_ALIGNMENT``: reshoot the colliding-wall strip."""
    if event.type is not EventType.VERTICAL_ALIGNMENT:
        raise ValueError("UPDATE_FOR_VERTICAL_ALIGNMENT requires an alignment event")
    return _update_local(vd_before, event, planes, z_after)


def _update_local(
    vd_before: VerticalDecomposition,
    event: Event,
    planes: Sequence[Plane],
    z_after: int | Scalar | str,
) -> VerticalDecomposition:
    arr_after = build_2d_arrangement_at_z(planes, z_after)
    local_keys = local_vertex_plane_keys(vd_before, arr_after, event, planes)
    walls = _incremental_walls(vd_before, arr_after, local_keys)
    cells = build_cells(arr_after, walls)
    return VerticalDecomposition(
        arrangement=arr_after, cells=tuple(cells), walls=tuple(walls)
    )


def _alignment_vertex_keys(
    event: Event, planes: Sequence[Plane]
) -> frozenset[tuple[int, ...]]:
    if len(event.line_ids) < 2:
        raise ValueError("alignment event needs two line ids")
    lookup = {line.id: line for line in compute_intersection_lines(planes)}
    keys: list[tuple[int, ...]] = []
    for line_id in event.line_ids:
        line = lookup.get(line_id)
        if line is None or line.plane_a is None or line.plane_b is None:
            raise ValueError(f"alignment line id {line_id} is missing source planes")
        keys.append(tuple(sorted((line.plane_a, line.plane_b))))
    return frozenset(keys)


def _incremental_walls(
    vd_before: VerticalDecomposition,
    arr_after: Arrangement2D,
    local_keys: frozenset[tuple[int, ...]],
) -> list[VerticalWall]:
    combo = _wall_combinatorics(vd_before)
    walls: list[VerticalWall] = []
    for vertex in arr_after.vertices:
        key = vertex_plane_key(arr_after, vertex)
        reshoot = key in local_keys or _must_reshoot(combo, key, local_keys)
        for ray in vertical_rays_from(vertex):
            hit = None
            if not reshoot:
                record = combo.get((key, ray.direction))
                if record is not None:
                    hit = _hit_from_record(arr_after, vertex, record, ray.direction)
            if hit is None:
                hit = first_hit(arr_after, ray)
            insert_decomposition_segment(walls, vertex, ray, hit)
    return finalize_decomposition_walls(arr_after, walls)


def _must_reshoot(
    combo: dict[tuple[tuple[int, ...], int], HitRecord],
    key: tuple[int, ...],
    local_keys: frozenset[tuple[int, ...]],
) -> bool:
    for direction in (1, -1):
        record = combo.get((key, direction))
        if record is not None and record.hit_vertex_key in local_keys:
            return True
    return False


def _wall_combinatorics(
    vd: VerticalDecomposition,
) -> dict[tuple[tuple[int, ...], int], HitRecord]:
    arrangement = vd.arrangement
    records: dict[tuple[tuple[int, ...], int], HitRecord] = {}
    for vertex in arrangement.vertices:
        key = vertex_plane_key(arrangement, vertex)
        for ray in vertical_rays_from(vertex):
            hit = first_hit(arrangement, ray)
            records[(key, ray.direction)] = _record_from_hit(arrangement, hit)
    return records


def _record_from_hit(arrangement: Arrangement2D, hit: Hit) -> HitRecord:
    if hit.unbounded:
        return HitRecord(unbounded=True, hit_plane_id=None, hit_vertex_key=None)
    plane_id = None
    if hit.line_index is not None:
        plane_id = arrangement.lines[hit.line_index].source_plane_id
    vertex_key = None
    if hit.vertex_id is not None:
        vertex_key = vertex_plane_key(arrangement, arrangement.vertices[hit.vertex_id])
    return HitRecord(
        unbounded=False, hit_plane_id=plane_id, hit_vertex_key=vertex_key
    )


def _hit_from_record(
    arrangement: Arrangement2D,
    vertex: Vertex,
    record: HitRecord,
    direction: int,
) -> Hit | None:
    if record.unbounded:
        return UNBOUNDED
    if record.hit_plane_id is None:
        return None
    line_index = _line_index_for_plane(arrangement, record.hit_plane_id)
    if line_index is None:
        return None
    try:
        y = y_at_x(arrangement.lines[line_index], vertex.point.x)
    except ValueError:
        return None
    dy = y - vertex.point.y
    if direction > 0 and dy <= 0:
        return None
    if direction < 0 and dy >= 0:
        return None
    point = Point2D(vertex.point.x, y)
    return Hit(
        unbounded=False,
        point=point,
        line_index=line_index,
        vertex_id=_vertex_id_at(arrangement, point),
    )


def _line_index_for_plane(arrangement: Arrangement2D, plane_id: int) -> int | None:
    for index, line in enumerate(arrangement.lines):
        if line.source_plane_id == plane_id:
            return index
    return None


def _vertex_id_at(arrangement: Arrangement2D, point: Point2D) -> int | None:
    for vertex in arrangement.vertices:
        if vertex.point == point:
            return vertex.id
    return None
