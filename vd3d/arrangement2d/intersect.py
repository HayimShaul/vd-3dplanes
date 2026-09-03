"""Exact 2D line–line intersection."""

from __future__ import annotations

from dataclasses import dataclass

from vd3d.geometry.linalg import det2, solve2
from vd3d.geometry.line2d import Line2D
from vd3d.geometry.points import Point2D


@dataclass(frozen=True, slots=True)
class ParallelLines:
    """Sentinel: distinct parallel lines, no intersection."""


@dataclass(frozen=True, slots=True)
class CoincidentLines:
    """Sentinel: the two equations describe the same line."""


PARALLEL = ParallelLines()
COINCIDENT = CoincidentLines()


def lines_parallel(left: Line2D, right: Line2D) -> bool:
    """True iff the normals are parallel (includes coincident lines)."""
    return det2(left.a, left.b, right.a, right.b) == 0


def intersect_lines_2d(
    left: Line2D,
    right: Line2D,
) -> Point2D | ParallelLines | CoincidentLines:
    """Intersect two 2D lines.

    Solves ``a x + b y = -c`` for both lines. A zero 2×2 determinant means
    the normals are parallel: the result is ``COINCIDENT`` if a sample point
    of ``left`` lies on ``right``, otherwise ``PARALLEL``.
    """
    solution = solve2(((left.a, left.b), (right.a, right.b)), (-left.c, -right.c))
    if solution is None:
        sample, _ = left.sample_points()
        if right.contains(sample):
            return COINCIDENT
        return PARALLEL
    point = Point2D(solution[0], solution[1])
    return point
