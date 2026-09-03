"""Formula oracle for simple 2D line arrangements."""

from vd3d.arrangement2d.sampling import random_simple_lines

__all__ = ["random_simple_lines", "simple_arrangement_counts"]


def simple_arrangement_counts(n: int) -> tuple[int, int, int]:
    """``(V, E, F)`` for a simple arrangement of ``n`` lines.

    Simple means no two parallel and no three concurrent:

    - ``V = n(n-1)/2``
    - ``E = n^2``
    - ``F = 1 + n(n+1)/2``
    """
    if n < 0:
        raise ValueError("n must be non-negative")
    vertices = n * (n - 1) // 2
    edges = n * n
    faces = 1 + n * (n + 1) // 2
    return vertices, edges, faces
