"""Exact geometric kernel.

2D and 3D types live here. This package must not import ``vd3d.viz``
or any floating-point plotting library.
"""

from vd3d.geometry.linalg import cross, det2, det3, solve2, solve3
from vd3d.geometry.scalar import Scalar, as_scalar

__all__ = [
    "Scalar",
    "as_scalar",
    "cross",
    "det2",
    "det3",
    "solve2",
    "solve3",
]
