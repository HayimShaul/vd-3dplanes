"""Exact scalar type for the geometric kernel.

Every coordinate and coefficient in the kernel is a ``fractions.Fraction``.
Floats are rejected so inexact values cannot enter the kernel by accident.
Convert to ``float`` only in ``vd3d.viz`` when drawing.
"""

from __future__ import annotations

from fractions import Fraction
from typing import TypeAlias

Scalar: TypeAlias = Fraction


def as_scalar(value: int | Fraction | str) -> Scalar:
    """Coerce ``value`` to an exact ``Fraction``.

    Accepts ``int``, ``Fraction``, or a ratio string such as ``"3/2"``.
    Rejects ``bool`` (a subclass of ``int``) and every floating type.
    """
    if isinstance(value, bool):
        raise TypeError("refusing bool; pass an explicit 0 or 1")
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str):
        return Fraction(value)
    raise TypeError(
        f"refusing inexact or unsupported type {type(value).__name__}: {value!r}"
    )
