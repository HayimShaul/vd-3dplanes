"""Exact duality: non-vertical planes ↔ dual points for the lower envelope."""

from __future__ import annotations

from collections.abc import Sequence

from vd3d.geometry.plane import Plane
from vd3d.geometry.scalar import Scalar


def plane_to_dual(plane: Plane) -> tuple[Scalar, Scalar, Scalar]:
    """Map ``a x + b y + c z + d = 0`` to ``(A, B, C)`` for ``z = A x + B y + C``.

    Requires ``c ≠ 0``: ``A = -a/c``, ``B = -b/c``, ``C = -d/c``.
    """
    if plane.c == 0:
        raise ValueError(
            f"plane P{plane.id} is vertical (c = 0); lower-envelope duality requires c != 0"
        )
    return (-plane.a / plane.c, -plane.b / plane.c, -plane.d / plane.c)


def dual_points(planes: Sequence[Plane]) -> tuple[tuple[Scalar, Scalar, Scalar], ...]:
    if not planes:
        raise ValueError("need at least one non-vertical plane")
    return tuple(plane_to_dual(plane) for plane in planes)
