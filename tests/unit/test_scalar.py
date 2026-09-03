from fractions import Fraction

import pytest

from vd3d.geometry.scalar import Scalar, as_scalar


def test_as_scalar_accepts_int():
    value = as_scalar(2)
    assert value == Fraction(2)
    assert isinstance(value, Scalar)


def test_as_scalar_accepts_fraction():
    value = as_scalar(Fraction(3, 2))
    assert value == Fraction(3, 2)
    assert value.numerator == 3
    assert value.denominator == 2


def test_as_scalar_accepts_ratio_string():
    assert as_scalar("3/2") == Fraction(3, 2)


def test_as_scalar_rejects_float():
    with pytest.raises(TypeError, match="unsupported"):
        as_scalar(0.5)


def test_as_scalar_rejects_bool():
    with pytest.raises(TypeError, match="bool"):
        as_scalar(True)
