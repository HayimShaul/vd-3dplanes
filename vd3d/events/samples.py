"""Hand fixtures for events (design Tests 15–22) and the Phase 7–8 slices."""

from __future__ import annotations

from vd3d.geometry.plane import Plane


def planes_one() -> list[Plane]:
    """A single non-vertical plane (design Test 19). Slice is horizontal ``y = const``."""
    return [Plane(id=1, a=0, b=1, c=1, d=0)]  # y + z = 0


def planes_two() -> list[Plane]:
    """Two intersecting non-vertical planes (design Test 20)."""
    return [
        Plane(id=1, a=1, b=0, c=1, d=-4),  # x + z = 4
        Plane(id=2, a=0, b=1, c=1, d=-5),  # y + z = 5
    ]


def planes_through_123() -> list[Plane]:
    """Three non-vertical planes meeting at ``(1, 2, 3)`` (design Test 15)."""
    return [
        Plane(id=1, a=1, b=0, c=1, d=-4),  # x + z = 4
        Plane(id=2, a=0, b=1, c=1, d=-5),  # y + z = 5
        Plane(id=3, a=1, b=1, c=1, d=-6),  # x + y + z = 6
    ]


def planes_parallel_family() -> list[Plane]:
    """Three pairwise-parallel planes (design Test 16)."""
    return [
        Plane(id=1, a=1, b=1, c=1, d=0),
        Plane(id=2, a=1, b=1, c=1, d=-1),
        Plane(id=3, a=1, b=1, c=1, d=-2),
    ]


def planes_four_through_123() -> list[Plane]:
    """Four planes through ``(1, 2, 3)``: four triples at the same ``z``."""
    return [
        *planes_through_123(),
        Plane(id=4, a=1, b=1, c=-1, d=0),  # x + y - z = 0 at (1,2,3)
    ]


def planes_vertical_and_slanted() -> list[Plane]:
    """Vertical ``x=0``, slanted ``y+z=0``, vertical ``y=0``. Triple at the origin."""
    return [
        Plane(id=1, a=1, b=0, c=0, d=0),  # x = 0
        Plane(id=2, a=0, b=1, c=1, d=0),  # y + z = 0
        Plane(id=3, a=0, b=1, c=0, d=0),  # y = 0
    ]


def planes_alignment_at_z2() -> list[Plane]:
    """``L12`` (x=z, y=0) and ``L34`` (x=4-z, y=3) align at ``z=2``."""
    return [
        Plane(id=1, a=1, b=1, c=-1, d=0),
        Plane(id=2, a=1, b=-1, c=-1, d=0),
        Plane(id=3, a=1, b=1, c=1, d=-7),
        Plane(id=4, a=1, b=-1, c=1, d=-1),
    ]


def planes_never_align() -> list[Plane]:
    """Parallel walls ``x = z`` and ``x = z + 1``."""
    return [
        Plane(id=1, a=1, b=1, c=-1, d=0),
        Plane(id=2, a=1, b=-1, c=-1, d=0),
        Plane(id=3, a=1, b=1, c=-1, d=-4),
        Plane(id=4, a=1, b=-1, c=-1, d=2),
    ]


def planes_alignment_blocked() -> list[Plane]:
    """Test 17 plus a third vertex on the open segment at ``z = 2``."""
    return [
        *planes_alignment_at_z2(),
        Plane(id=5, a=0, b=2, c=1, d=-5),
        Plane(id=6, a=1, b=1, c=-1, d="-3/2"),
    ]
