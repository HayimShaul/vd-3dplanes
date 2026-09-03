from fractions import Fraction

from vd3d.geometry.linalg import cross, det2, det3, solve2, solve3
from vd3d.geometry.scalar import Scalar


def test_det2_exact():
    assert det2(1, 2, 3, 4) == Fraction(-2)
    assert isinstance(det2(1, 2, 3, 4), Scalar)


def test_det3_exact():
    matrix = (
        (1, 2, 3),
        (0, 1, 4),
        (5, 6, 0),
    )
    # 1*(0-24) - 2*(0-20) + 3*(0-5) = -24 + 40 - 15 = 1
    assert det3(matrix) == Fraction(1)


def test_cross_exact_and_orthogonal():
    result = cross((1, 0, 0), (0, 1, 0))
    assert result == (Fraction(0), Fraction(0), Fraction(1))
    assert all(isinstance(component, Scalar) for component in result)


def test_solve2_exact():
    # x + 2y = 5, 3x + 4y = 11  =>  (x, y) = (1, 2)
    result = solve2(((1, 2), (3, 4)), (5, 11))
    assert result == (Fraction(1), Fraction(2))
    assert all(isinstance(component, Scalar) for component in result)


def test_solve2_fraction_rhs():
    # 2x = 1  =>  x = 1/2; use a triangular system
    result = solve2(((2, 0), (0, 3)), (1, 1))
    assert result == (Fraction(1, 2), Fraction(1, 3))


def test_solve2_singular_is_none():
    assert solve2(((1, 2), (2, 4)), (1, 2)) is None


def test_solve3_exact():
    # Identity system keeps the right-hand side.
    result = solve3(
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        (Fraction(2, 3), 4, -1),
    )
    assert result == (Fraction(2, 3), Fraction(4), Fraction(-1))
    assert all(isinstance(component, Scalar) for component in result)


def test_solve3_known_system():
    # [[2, 1, 0], [0, 2, 1], [0, 0, 2]] [x,y,z] = [1, 1, 2]
    # z = 1, 2y + z = 1 => y = 0, 2x + y = 1 => x = 1/2
    result = solve3(((2, 1, 0), (0, 2, 1), (0, 0, 2)), (1, 1, 2))
    assert result == (Fraction(1, 2), Fraction(0), Fraction(1))


def test_solve3_singular_is_none():
    assert solve3(((1, 0, 0), (0, 1, 0), (2, 0, 0)), (1, 1, 1)) is None
