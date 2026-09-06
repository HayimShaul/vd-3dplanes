"""Matplotlib helpers and the interactive review viewer.

Kernel packages must not import this package. Convert ``Scalar`` to
``float`` only through ``vd3d.viz.convert.to_float``.

Human review of 3D scenes::

    python -m vd3d.viz.viewer --phase 1

Human review of 2D arrangements::

    python -m vd3d.viz.viewer --phase 2 --n 6
    python -m vd3d.viz.gallery --step 2

Human review of 2D vertical decompositions::

    python -m vd3d.viz.viewer --phase 3 --n 6
    python -m vd3d.viz.gallery --step 3

Human review of query-line zones::

    python -m vd3d.viz.viewer --phase 4 --n 8
    python -m vd3d.viz.gallery --step 4

Human review of intersection lines and triple events::

    python -m vd3d.viz.viewer --phase 5
    python -m vd3d.viz.gallery --step 5

Human review of alignment events::

    python -m vd3d.viz.viewer --phase 6
    python -m vd3d.viz.gallery --step 6

Human review of the combined event list and slice VD::

    python -m vd3d.viz.viewer --phase 7
    python -m vd3d.viz.gallery --step 7
"""

from vd3d.viz.convert import to_float
from vd3d.viz.record import record_figure

__all__ = ["record_figure", "to_float"]
