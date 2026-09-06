"""Sample simple 2D lines with exact integer coefficients."""

from __future__ import annotations

import random

from vd3d.arrangement2d.intersect import COINCIDENT, PARALLEL, intersect_lines_2d
from vd3d.arrangement2d.types import Arrangement2D
from vd3d.geometry.linalg import det3
from vd3d.geometry.line2d import Line2D


def _three_concurrent(a: Line2D, b: Line2D, c: Line2D) -> bool:
    return det3(((a.a, a.b, a.c), (b.a, b.b, b.c), (c.a, c.b, c.c))) == 0


def random_simple_lines(rng: random.Random, n: int, *, coeff: int = 8) -> list[Line2D]:
    """``n`` lines, no parallels, no coincidences, no three concurrent."""
    if n < 0:
        raise ValueError("n must be non-negative")
    lines: list[Line2D] = []
    attempts = 0
    limit = 5000 + 200 * n * n
    while len(lines) < n:
        attempts += 1
        if attempts > limit:
            raise RuntimeError(f"failed to sample {n} simple lines")
        a = rng.randint(-coeff, coeff)
        b = rng.randint(-coeff, coeff)
        c = rng.randint(-coeff, coeff)
        if a == 0 and b == 0:
            continue
        candidate = Line2D(a=a, b=b, c=c, id=len(lines))
        ok = True
        for existing in lines:
            result = intersect_lines_2d(candidate, existing)
            if result is PARALLEL or result is COINCIDENT:
                ok = False
                break
        if not ok:
            continue
        if any(
            _three_concurrent(candidate, lines[i], lines[j])
            for i in range(len(lines))
            for j in range(i + 1, len(lines))
        ):
            continue
        lines.append(candidate)
    return lines


def random_query_missing_vertices(
    rng: random.Random,
    arrangement: Arrangement2D,
    *,
    coeff: int = 8,
) -> Line2D:
    """A query line that is not coincident with any arrangement line and misses vertices."""
    attempts = 0
    limit = 4000 + 200 * max(len(arrangement.lines), 1)
    while attempts < limit:
        attempts += 1
        a = rng.randint(-coeff, coeff)
        b = rng.randint(-coeff, coeff)
        c = rng.randint(-coeff, coeff)
        if a == 0 and b == 0:
            continue
        query = Line2D(a=a, b=b, c=c)
        degenerate = False
        for line in arrangement.lines:
            result = intersect_lines_2d(query, line)
            if result is COINCIDENT:
                degenerate = True
                break
        if degenerate:
            continue
        if any(query.contains(vertex.point) for vertex in arrangement.vertices):
            continue
        return query
    raise RuntimeError("failed to sample a general-position query line")
