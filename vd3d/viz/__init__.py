"""Matplotlib helpers and the interactive review viewer.

Kernel packages must not import this package. Convert ``Scalar`` to
``float`` only through ``vd3d.viz.convert.to_float``.

Human review of 3D scenes::

    python -m vd3d.viz.viewer --phase 1

Human review of 2D arrangements::

    python -m vd3d.viz.viewer --phase 2
    python -m vd3d.viz.gallery --step 2
"""

from vd3d.viz.convert import to_float
from vd3d.viz.record import record_figure

__all__ = ["record_figure", "to_float"]
