"""Matplotlib helpers and the human-review gallery.

Kernel packages must not import this package. Convert ``Scalar`` to
``float`` only through ``vd3d.viz.convert.to_float``.
"""

from vd3d.viz.convert import to_float
from vd3d.viz.record import record_figure

__all__ = ["record_figure", "to_float"]
