"""Exact helpers for 2D vertical (x = const) geometry."""

from __future__ import annotations

from vd3d.geometry.line2d import Line2D
from vd3d.geometry.scalar import Scalar


def line_is_vertical(line: Line2D) -> bool:
    """True iff the line is parallel to the y-axis (``b = 0``, ``x = const``)."""
    return line.b == 0


def x_of_vertical_line(line: Line2D) -> Scalar:
    """The constant ``x`` of a vertical line ``a x + c = 0``."""
    if line.b != 0:
        raise ValueError("line is not vertical")
    return -line.c / line.a


def y_at_x(line: Line2D, x: Scalar) -> Scalar:
    """``y`` on a non-vertical line at the given ``x``."""
    if line.b == 0:
        raise ValueError("vertical line has no y(x)")
    return -(line.a * x + line.c) / line.b
