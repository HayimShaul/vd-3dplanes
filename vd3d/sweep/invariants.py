"""Sweep invariants from design §21."""

from __future__ import annotations

from vd3d.cells3d.types import ActiveCell
from vd3d.sweep.matching import cell_signature, same_combinatorics, signatures_unique
from vd3d.sweep.types import CellSignature
from vd3d.vertical_decomposition.types import VerticalDecomposition


def active_count_matches_2d(
    active: dict[CellSignature, ActiveCell], vd: VerticalDecomposition
) -> bool:
    return len(active) == len(vd.cells)


def active_keys_match_vd(
    active: dict[CellSignature, ActiveCell], vd: VerticalDecomposition
) -> bool:
    return set(active) == {cell_signature(vd, cell) for cell in vd.cells}


def verify_active_against_vd(
    active: dict[CellSignature, ActiveCell], vd: VerticalDecomposition
) -> None:
    checks = (
        ("signatures unique", lambda: signatures_unique(vd)),
        (
            "|active 3D cells| == |current 2D cells|",
            lambda: active_count_matches_2d(active, vd),
        ),
        ("active keys are the current 2D signatures", lambda: active_keys_match_vd(active, vd)),
    )
    for name, pred in checks:
        if not pred():
            raise AssertionError(f"sweep invariant failed: {name}")


def equivalent_vd(
    incremental: VerticalDecomposition, reference: VerticalDecomposition
) -> bool:
    """True iff the two VDs have the same cell-signature multiset (design §12)."""
    return (
        same_combinatorics(incremental, reference)
        and len(incremental.cells) == len(reference.cells)
    )


def verify_interval_matches_vd(
    vd: VerticalDecomposition, expected: VerticalDecomposition
) -> None:
    if not same_combinatorics(vd, expected):
        raise AssertionError(
            "sweep invariant failed: CURRENT_VD != independently computed VD"
        )


def verify_incremental_matches_recompute(
    incremental: VerticalDecomposition, reference: VerticalDecomposition
) -> None:
    if not equivalent_vd(incremental, reference):
        raise AssertionError(
            "sweep invariant failed: incremental VD != compute_vd_at_z(z+)"
        )
