"""Exact 2x2 and 3x3 linear algebra over ``Scalar``.

All inputs are coerced with ``as_scalar``. Results stay ``Fraction``.
A singular system returns ``None`` rather than raising.
"""

from __future__ import annotations

from vd3d.geometry.scalar import Scalar, as_scalar

Vec3 = tuple[Scalar, Scalar, Scalar]
Mat2 = tuple[tuple[Scalar, Scalar], tuple[Scalar, Scalar]]
Mat3 = tuple[Vec3, Vec3, Vec3]


def _s(value: int | Scalar | str) -> Scalar:
    return as_scalar(value)


def det2(
    a: int | Scalar | str,
    b: int | Scalar | str,
    c: int | Scalar | str,
    d: int | Scalar | str,
) -> Scalar:
    """Determinant of ``[[a, b], [c, d]]``."""
    return _s(a) * _s(d) - _s(b) * _s(c)


def det3(matrix: Mat3 | tuple[tuple[int | Scalar | str, ...], ...]) -> Scalar:
    """Determinant of a 3x3 matrix given as three rows."""
    (a, b, c), (d, e, f), (g, h, i) = ((_s(x) for x in row) for row in matrix)
    return (
        a * (e * i - f * h)
        - b * (d * i - f * g)
        + c * (d * h - e * g)
    )


def cross(
    u: tuple[int | Scalar | str, int | Scalar | str, int | Scalar | str],
    v: tuple[int | Scalar | str, int | Scalar | str, int | Scalar | str],
) -> Vec3:
    """3D cross product."""
    ux, uy, uz = (_s(x) for x in u)
    vx, vy, vz = (_s(x) for x in v)
    return (
        uy * vz - uz * vy,
        uz * vx - ux * vz,
        ux * vy - uy * vx,
    )


def solve2(
    matrix: tuple[tuple[int | Scalar | str, int | Scalar | str], tuple[int | Scalar | str, int | Scalar | str]],
    rhs: tuple[int | Scalar | str, int | Scalar | str],
) -> tuple[Scalar, Scalar] | None:
    """Solve a 2x2 system by Cramer's rule. ``None`` if the determinant is 0."""
    (a, b), (c, d) = matrix
    e, f = rhs
    a, b, c, d, e, f = (_s(a), _s(b), _s(c), _s(d), _s(e), _s(f))
    determinant = det2(a, b, c, d)
    if determinant == 0:
        return None
    x = det2(e, b, f, d) / determinant
    y = det2(a, e, c, f) / determinant
    return (x, y)


def solve3(
    matrix: Mat3 | tuple[tuple[int | Scalar | str, ...], ...],
    rhs: tuple[int | Scalar | str, int | Scalar | str, int | Scalar | str],
) -> Vec3 | None:
    """Solve a 3x3 system by Cramer's rule. ``None`` if the determinant is 0."""
    rows = tuple(tuple(_s(x) for x in row) for row in matrix)
    bx, by, bz = (_s(x) for x in rhs)
    determinant = det3(rows)
    if determinant == 0:
        return None

    # Cramer's rule: replace one column at a time.
    m0: Mat3 = (
        (bx, rows[0][1], rows[0][2]),
        (by, rows[1][1], rows[1][2]),
        (bz, rows[2][1], rows[2][2]),
    )
    m1: Mat3 = (
        (rows[0][0], bx, rows[0][2]),
        (rows[1][0], by, rows[1][2]),
        (rows[2][0], bz, rows[2][2]),
    )
    m2: Mat3 = (
        (rows[0][0], rows[0][1], bx),
        (rows[1][0], rows[1][1], by),
        (rows[2][0], rows[2][1], bz),
    )
    return (
        det3(m0) / determinant,
        det3(m1) / determinant,
        det3(m2) / determinant,
    )
