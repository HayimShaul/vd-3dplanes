"""Interactive 3D viewer for a finished vertical decomposition.

    python -m vd3d examples/three_planes.txt --gui
    python -m vd3d --random-planes 4 --seed 42 --gui

Keys:
  n / p       next / previous cell
  + / -       zoom in / out (shrink / grow the viewing cube)
  ↑ / ↓       rotate around the x-axis (elevation)
  ← / →       rotate around the y-axis (azimuth)
  q / Esc     quit
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from vd3d.io import format_plane
from vd3d.sweep.algorithm import SweepResult
from vd3d.viz.plot_cells3d import (
    ViewBox,
    default_view_half,
    draw_bbox_wire,
    draw_cell_red,
    draw_intersection_lines,
    draw_planes_grey,
)
from vd3d.viz.plot_geometry import set_equal_3d

_HELP = (
    "n/p cell   +/- zoom   ↑↓ elev   ←→ azim   q quit"
)


class DecompositionGui:
    """Show input planes (grey), intersection lines (dark), one cell (red)."""

    def __init__(self, result: SweepResult, *, seed: int | None = None) -> None:
        if not result.cells:
            raise ValueError("no 3D cells to display")
        self.result = result
        self.seed = seed
        self.index = 0
        self.half = default_view_half(result)
        self.elev = 22.0
        self.azim = -60.0
        self.fig = plt.figure(figsize=(9.5, 7.5))
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self._draw()

    def _title(self) -> str:
        cell = self.result.cells[self.index]
        parts = [
            f"cell {cell.id} ({self.index + 1}/{len(self.result.cells)})",
            f"box=±{self.half:g}",
        ]
        if self.seed is not None:
            parts.append(f"seed={self.seed}")
        return "  ·  ".join(parts)

    def _draw(self) -> None:
        self.ax.clear()
        box = ViewBox(self.half)
        draw_planes_grey(self.ax, self.result.planes, box)
        draw_intersection_lines(self.ax, self.result.planes, box)
        draw_cell_red(self.ax, self.result, self.result.cells[self.index], box)
        draw_bbox_wire(self.ax, box)
        set_equal_3d(self.ax, lim=box.lim)
        self.ax.view_init(elev=self.elev, azim=self.azim)
        cell = self.result.cells[self.index]
        caption = (
            f"P{cell.floor.id}:{format_plane(cell.floor)}"
            if cell.floor is not None
            else "floor=unbounded"
        )
        if cell.ceiling is not None:
            caption += f"  ceil=P{cell.ceiling.id}:{format_plane(cell.ceiling)}"
        else:
            caption += "  ceil=unbounded"
        self.ax.set_title(f"{self._title()}\n{caption}", fontsize=10)
        self.fig.subplots_adjust(bottom=0.08, top=0.90)
        self.fig.text(0.5, 0.015, _HELP, ha="center", va="bottom", fontsize=9)
        manager = getattr(self.fig.canvas, "manager", None)
        if manager is not None and hasattr(manager, "set_window_title"):
            manager.set_window_title(f"VD3D — {self._title()}")
        self.fig.canvas.draw_idle()

    def _on_key(self, event) -> None:
        key = event.key
        n = len(self.result.cells)
        if key in {"n", " "}:
            self.index = (self.index + 1) % n
            self._draw()
        elif key == "p":
            self.index = (self.index - 1) % n
            self._draw()
        elif key in {"+", "=", "plus"}:
            self.half = max(0.5, self.half * 0.8)
            self._draw()
        elif key in {"-", "minus", "_"}:
            self.half = min(200.0, self.half * 1.25)
            self._draw()
        elif key == "up":
            self.elev = (self.elev + 8.0) % 360.0
            self.ax.view_init(elev=self.elev, azim=self.azim)
            self.fig.canvas.draw_idle()
        elif key == "down":
            self.elev = (self.elev - 8.0) % 360.0
            self.ax.view_init(elev=self.elev, azim=self.azim)
            self.fig.canvas.draw_idle()
        elif key == "left":
            self.azim = (self.azim - 8.0) % 360.0
            self.ax.view_init(elev=self.elev, azim=self.azim)
            self.fig.canvas.draw_idle()
        elif key == "right":
            self.azim = (self.azim + 8.0) % 360.0
            self.ax.view_init(elev=self.elev, azim=self.azim)
            self.fig.canvas.draw_idle()
        elif key in {"q", "escape"}:
            plt.close(self.fig)

    def show(self) -> None:
        plt.show()


def show_decomposition(result: SweepResult, *, seed: int | None = None) -> DecompositionGui:
    """Open the GUI. Requires a non-Agg matplotlib backend."""
    if not plt.isinteractive() and plt.get_backend().lower() == "agg":
        raise RuntimeError(
            "matplotlib backend is Agg (no window). "
            "Unset MPLBACKEND and install a GUI backend (e.g. TkAgg)."
        )
    gui = DecompositionGui(result, seed=seed)
    gui.show()
    return gui
