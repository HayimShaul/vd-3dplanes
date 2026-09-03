"""The only place kernel scalars become floats."""

from __future__ import annotations

from fractions import Fraction


def to_float(value: int | Fraction) -> float:
    """Convert an exact kernel value to ``float`` for plotting.

    Never call this from ``vd3d.geometry`` or any other kernel package.
    """
    if isinstance(value, bool):
        raise TypeError("refusing bool; pass an explicit 0 or 1")
    if isinstance(value, (int, Fraction)):
        return float(value)
    raise TypeError(
        f"to_float expects int or Fraction, not {type(value).__name__}"
    )
