"""Brute-force checkers for incremental 2D VD updates."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.events.slice import compute_vd_at_z
from vd3d.events.types import Event
from vd3d.geometry.plane import Plane
from vd3d.geometry.scalar import Scalar
from vd3d.sweep.invariants import equivalent_vd
from vd3d.sweep.matching import cell_signature, same_combinatorics
from vd3d.sweep.update import update_2d_decomposition
from vd3d.vertical_decomposition.invariants import verify_vd_invariants
from vd3d.vertical_decomposition.types import VerticalDecomposition, VerticalWall


def wall_keys(walls: Sequence[VerticalWall]) -> tuple:
    """Canonical wall geometry, ignoring ids."""
    return tuple(
        sorted(
            (
                wall.x,
                wall.y_min is None,
                wall.y_min if wall.y_min is not None else 0,
                wall.y_max is None,
                wall.y_max if wall.y_max is not None else 0,
            )
            for wall in walls
        )
    )


def cell_bound_keys(vd: VerticalDecomposition) -> tuple:
    """Canonical trapezoid bounds (x-range + supporting line indices)."""
    return tuple(
        sorted(
            (
                cell.left_x is None,
                cell.left_x if cell.left_x is not None else 0,
                cell.right_x is None,
                cell.right_x if cell.right_x is not None else 0,
                cell.lower_line if cell.lower_line is not None else -1,
                cell.upper_line if cell.upper_line is not None else -1,
            )
            for cell in vd.cells
        )
    )


def assert_equivalent_vd(
    incremental: VerticalDecomposition,
    reference: VerticalDecomposition,
    *,
    label: str = "",
) -> None:
    """Incremental result matches a from-scratch VD: invariants, signatures, walls."""
    prefix = f"{label}: " if label else ""
    verify_vd_invariants(incremental)
    verify_vd_invariants(reference)
    if not equivalent_vd(incremental, reference):
        inc = sorted(cell_signature(incremental, cell) for cell in incremental.cells)
        ref = sorted(cell_signature(reference, cell) for cell in reference.cells)
        raise AssertionError(
            f"{prefix}incremental VD != recomputed VD\n"
            f"  incremental signatures ({len(inc)}): {inc}\n"
            f"  reference signatures ({len(ref)}): {ref}"
        )
    if wall_keys(incremental.walls) != wall_keys(reference.walls):
        raise AssertionError(f"{prefix}incremental walls != recomputed walls")
    if cell_bound_keys(incremental) != cell_bound_keys(reference):
        raise AssertionError(f"{prefix}incremental cell bounds != recomputed cell bounds")


def assert_update_matches_recompute(
    vd_before: VerticalDecomposition,
    event: Event,
    planes: Sequence[Plane],
    z_after: Scalar,
) -> VerticalDecomposition:
    """Run the incremental handler and compare to ``compute_vd_at_z(z+)``."""
    incremental = update_2d_decomposition(vd_before, event, planes, z_after)
    reference = compute_vd_at_z(planes, z_after)
    assert_equivalent_vd(
        incremental, reference, label=f"{event.type.name} {event.stable_id} z+={z_after}"
    )
    if same_combinatorics(vd_before, incremental) and not same_combinatorics(
        vd_before, reference
    ):
        raise AssertionError(
            f"{event.stable_id}: incremental did not apply a combinatorial change"
        )
    return incremental
