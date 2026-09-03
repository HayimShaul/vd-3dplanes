"""A drawable review scene."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from matplotlib.figure import Figure


@dataclass(frozen=True)
class Scene:
    name: str
    title: str
    caption: str
    figsize: tuple[float, float]
    draw: Callable[[Figure], None]
