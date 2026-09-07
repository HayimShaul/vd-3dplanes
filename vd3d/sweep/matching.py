"""Match 2D VD cells across an event (design §13)."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.cells3d.walls import side_planes, supporting_plane_ids
from vd3d.events.types import AlignmentGeometry, Event, EventType
from vd3d.geometry.points import Point2D, Point3D
from vd3d.sweep.types import CellMatch, CellSignature
from vd3d.vertical_decomposition.invariants import point_in_cell
from vd3d.vertical_decomposition.types import VDCell2D, VerticalDecomposition


def cell_signature(vd: VerticalDecomposition, cell: VDCell2D) -> CellSignature:
    """Canonical key: supporting planes + boundedness + wall-generating planes."""
    lower, upper = supporting_plane_ids(vd, cell)
    return CellSignature(
        lower_plane_id=lower,
        upper_plane_id=upper,
        left_unbounded=cell.left_x is None,
        right_unbounded=cell.right_x is None,
        left_planes=side_planes(vd, cell, cell.left_x),
        right_planes=side_planes(vd, cell, cell.right_x),
    )


def signatures_unique(vd: VerticalDecomposition) -> bool:
    keys = [cell_signature(vd, cell) for cell in vd.cells]
    return len(keys) == len(set(keys))


def same_combinatorics(vd_a: VerticalDecomposition, vd_b: VerticalDecomposition) -> bool:
    """True iff the two slices have the same multiset of cell signatures."""
    return _signature_multiset(vd_a) == _signature_multiset(vd_b)


def locate_cell(
    vd: VerticalDecomposition, point: Point2D, *, closed: bool = False
) -> VDCell2D | None:
    """The unique cell containing ``point``, or ``None`` if none / several."""
    hits = [
        cell
        for cell in vd.cells
        if point_in_cell(vd, cell, point, closed=closed)
    ]
    if len(hits) != 1:
        return None
    return hits[0]


def match_cells(
    vd_before: VerticalDecomposition,
    vd_after: VerticalDecomposition,
) -> tuple[CellMatch, ...]:
    """Mutual geometric point-in-cell, restricted to the same supporting planes.

    A before-cell continues as an after-cell when each representative lies in
    the other cell and the lower/upper source planes agree. Far from an event
    this is 1-1; the event neighbourhood is the unmatched remainder.
    """
    matches: list[CellMatch] = []
    used_after: set[int] = set()
    for cell_before in vd_before.cells:
        cell_after = locate_cell(vd_after, cell_before.representative, closed=False)
        if cell_after is None or cell_after.id in used_after:
            continue
        back = locate_cell(vd_before, cell_after.representative, closed=False)
        if back is None or back.id != cell_before.id:
            continue
        if supporting_plane_ids(vd_before, cell_before) != supporting_plane_ids(
            vd_after, cell_after
        ):
            continue
        matches.append(CellMatch(before=cell_before, after=cell_after))
        used_after.add(cell_after.id)
    return tuple(matches)


def unmatched_before(
    vd_before: VerticalDecomposition, matches: Sequence[CellMatch]
) -> tuple[VDCell2D, ...]:
    matched = {item.before.id for item in matches}
    return tuple(cell for cell in vd_before.cells if cell.id not in matched)


def unmatched_after(
    vd_after: VerticalDecomposition, matches: Sequence[CellMatch]
) -> tuple[VDCell2D, ...]:
    matched = {item.after.id for item in matches}
    return tuple(cell for cell in vd_after.cells if cell.id not in matched)


def event_points_2d(event: Event) -> tuple[Point2D, ...]:
    """The xy locations that define the event neighbourhood."""
    data = event.geometric_data
    if event.type is EventType.TRIPLE_INTERSECTION:
        if not isinstance(data, Point3D):
            raise TypeError("triple event geometric_data must be a Point3D")
        return (Point2D(data.x, data.y),)
    if event.type is EventType.VERTICAL_ALIGNMENT:
        if not isinstance(data, AlignmentGeometry):
            raise TypeError("alignment event geometric_data must be AlignmentGeometry")
        return (
            Point2D(data.point_a.x, data.point_a.y),
            Point2D(data.point_b.x, data.point_b.y),
        )
    raise ValueError(f"unknown event type: {event.type!r}")


def cell_near_event(
    vd: VerticalDecomposition, cell: VDCell2D, event: Event
) -> bool:
    """True iff the closed cell meets the event's xy feature or its x-strip."""
    points = event_points_2d(event)
    if any(point_in_cell(vd, cell, point, closed=True) for point in points):
        return True
    if any(_x_range_contains(cell, point.x) for point in points):
        return True
    if event.type is EventType.VERTICAL_ALIGNMENT and len(points) == 2:
        y0, y1 = points[0].y, points[1].y
        lo, hi = (y0, y1) if y0 <= y1 else (y1, y0)
        x = points[0].x
        sample = Point2D(x, (lo + hi) / 2)
        if point_in_cell(vd, cell, sample, closed=True):
            return True
    return False


def _x_range_contains(cell: VDCell2D, x) -> bool:
    if cell.left_x is not None and x < cell.left_x:
        return False
    if cell.right_x is not None and x > cell.right_x:
        return False
    return True


def local_cell_ids(vd: VerticalDecomposition, event: Event) -> frozenset[int]:
    """Cells whose closed geometry meets the event, plus their neighbors."""
    seed = [cell.id for cell in vd.cells if cell_near_event(vd, cell, event)]
    extra: set[int] = set(seed)
    for cell_id in seed:
        extra.update(vd.cells[cell_id].neighbors)
    return frozenset(extra)


def local_cell_ids_for_events(
    vd: VerticalDecomposition, events: Sequence[Event]
) -> frozenset[int]:
    """Union of ``local_cell_ids`` over a simultaneous group."""
    extra: set[int] = set()
    for event in events:
        extra |= local_cell_ids(vd, event)
    return frozenset(extra)


def _signature_multiset(
    vd: VerticalDecomposition,
) -> tuple[CellSignature, ...]:
    return tuple(sorted((cell_signature(vd, cell) for cell in vd.cells), key=_sig_sort))


def _sig_sort(sig: CellSignature) -> tuple:
    return (
        sig.lower_plane_id is None,
        sig.lower_plane_id or 0,
        sig.upper_plane_id is None,
        sig.upper_plane_id or 0,
        sig.left_unbounded,
        sig.right_unbounded,
        sig.left_planes,
        sig.right_planes,
    )
