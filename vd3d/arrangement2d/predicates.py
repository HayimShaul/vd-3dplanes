"""Exact orientation predicates. No ``atan2``, no floats."""

from __future__ import annotations

from functools import cmp_to_key

from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D
from vd3d.geometry.scalar import Scalar

Direction = tuple[Scalar, Scalar]


def line_direction(line: Line2D) -> Direction:
    """A direction vector along ``line``: ``(-b, a)``.

    This is the normal ``(a, b)`` rotated 90° clockwise, so travelling
    along it puts the positive half-plane (``eval > 0``) on the right.
    """
    return (-line.b, line.a)


def parameter_on_line(line: Line2D, point: Point2D) -> Scalar:
    """1D sort key of ``point`` along ``line``.

    Equals the dot product with ``line_direction``. On a single line this
    is strictly monotone in the forward direction.
    """
    dx, dy = line_direction(line)
    return dx * point.x + dy * point.y


def cross2(u: Direction, v: Direction) -> Scalar:
    """2D cross product ``u_x v_y - u_y v_x``. Positive iff ``v`` is CCW from ``u``."""
    return u[0] * v[1] - u[1] * v[0]


def dot2(u: Direction, v: Direction) -> Scalar:
    return u[0] * v[0] + u[1] * v[1]


def same_direction(u: Direction, v: Direction) -> bool:
    """True iff ``u`` and ``v`` are positive scalar multiples (not opposites)."""
    return cross2(u, v) == 0 and dot2(u, v) > 0


def left_normal(direction: Direction) -> Direction:
    """Rotate 90° CCW: the face-on-the-left side of a half-edge."""
    dx, dy = direction
    return (-dy, dx)


def quadrant(direction: Direction) -> int:
    """Quadrant index increasing CCW from the +x axis, including the start ray.

    - 0: ``dx > 0, dy >= 0``  (includes +x)
    - 1: ``dx <= 0, dy > 0``  (includes +y)
    - 2: ``dx < 0, dy <= 0``  (includes -x)
    - 3: ``dx >= 0, dy < 0``  (includes -y)
    """
    dx, dy = direction
    if dx == 0 and dy == 0:
        raise ValueError("zero direction")
    if dx > 0 and dy >= 0:
        return 0
    if dx <= 0 and dy > 0:
        return 1
    if dx < 0 and dy <= 0:
        return 2
    return 3


def cmp_direction_ccw(u: Direction, v: Direction) -> int:
    """``-1`` if ``u`` precedes ``v`` in CCW order from +x, ``0`` if the same ray, else ``1``."""
    if same_direction(u, v):
        return 0
    qu, qv = quadrant(u), quadrant(v)
    if qu != qv:
        return -1 if qu < qv else 1
    # Same quadrant: ``v`` is CCW from ``u`` iff cross(u, v) > 0.
    cr = cross2(u, v)
    if cr > 0:
        return -1
    if cr < 0:
        return 1
    return 0


def sort_directions_ccw(directions: list[Direction]) -> list[Direction]:
    """Stable CCW sort of directions, starting from the +x axis."""
    return sorted(directions, key=cmp_to_key(cmp_direction_ccw))
