"""Slice planes at a fixed ``z`` and build the 2D VD (design §10)."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.arrangement2d import Arrangement2D, build_line_arrangement
from vd3d.events.all import choose_z_below_all_events, generate_all_events
from vd3d.events.ids import require_unique_plane_ids
from vd3d.events.types import Event
from vd3d.geometry.intersections import build_slice_lines
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.plane import Plane
from vd3d.geometry.scalar import Scalar
from vd3d.vertical_decomposition import VerticalDecomposition, compute_vertical_decomposition


def build_2d_arrangement_at_z(
    planes: Sequence[Plane], z: int | Scalar | str
) -> Arrangement2D:
    """``BUILD_2D_ARRANGEMENT(P, z)``: slice every plane, then the 2D DCEL."""
    require_unique_plane_ids(planes)
    raw = build_slice_lines(list(planes), z)
    lines = [
        Line2D(
            a=line.a,
            b=line.b,
            c=line.c,
            source_plane_id=line.source_plane_id,
            id=index,
        )
        for index, line in enumerate(raw)
    ]
    return build_line_arrangement(lines)


def compute_vd_at_z(
    planes: Sequence[Plane], z: int | Scalar | str
) -> VerticalDecomposition:
    """``COMPUTE_VD_AT_Z(P, z) = VD(arrangement(slices))``."""
    return compute_vertical_decomposition(build_2d_arrangement_at_z(planes, z))


def initial_slice(
    planes: Sequence[Plane], events: Sequence[Event] | None = None
) -> VerticalDecomposition:
    """``INITIAL_SLICE``: 2D VD at a ``z`` strictly below every event."""
    if events is None:
        events = generate_all_events(planes)
    return compute_vd_at_z(planes, choose_z_below_all_events(events))
