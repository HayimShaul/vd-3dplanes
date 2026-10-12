"""Interactive 3D viewer for the hypograph VD.

    python -m vd3d.lower_envelope examples/three_planes.txt --vd --gui
    python -m vd3d.lower_envelope examples/three_planes.txt --exact --vd --gui

Keys:
  n / p       next / previous cell
  + / -       zoom
  ↑ / ↓       elevation
  ← / →       azimuth
  q / Esc     quit
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from vd3d.io import format_plane
from vd3d.lower_envelope.vd import HypographCompare, HypographVD
from vd3d.lower_envelope.vd_float import FloatHypographVD
from vd3d.viz.plot_geometry import set_equal_3d
from vd3d.viz.plot_hypograph import draw_hypograph_scene, hypograph_view_half
from vd3d.viz.plot_cells3d import ViewBox

_HELP = "n/p cell   +/- zoom   ↑↓ elev   ←→ azim   q quit"

AnyVD = FloatHypographVD | HypographVD


class HypographGui:
    """Grey planes, blue envelope, all hypograph cells faint, selection in red."""

    def __init__(
        self,
        vd: AnyVD,
        *,
        seed: int | None = None,
        exact: bool = False,
        compare: HypographCompare | None = None,
    ) -> None:
        if not vd.cells:
            raise ValueError("no hypograph cells to display")
        self.vd = vd
        self.seed = seed
        self.exact = exact
        self.compare = compare
        self.index = 0
        self.half = hypograph_view_half(vd)
        self.elev = 22.0
        self.azim = -60.0
        self.fig = plt.figure(figsize=(9.5, 7.5))
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self._draw()

    def _title(self) -> str:
        cell = self.vd.cells[self.index]
        kind = "exact" if self.exact else "float"
        parts = [
            f"{kind} hypograph  cell {cell.id} ({self.index + 1}/{len(self.vd.cells)})",
            f"box=±{self.half:g}",
        ]
        if self.compare is not None:
            parts.append("compare=" + ("OK" if self.compare.ok else "MISMATCH"))
        if self.seed is not None:
            parts.append(f"seed={self.seed}")
        return "  ·  ".join(parts)

    def _caption(self) -> str:
        cell = self.vd.cells[self.index]
        trap = self.vd.trapezoids[self.index]
        ceil = cell.ceiling
        assert ceil is not None
        lx = "-∞" if trap.left_x is None else f"{trap.left_x}"
        rx = "+∞" if trap.right_x is None else f"{trap.right_x}"
        return (
            f"floor=unbounded  ceil=P{ceil.id}:{format_plane(ceil)}  "
            f"x=({lx},{rx})  walls={len(cell.vertical_walls)}"
        )

    def _draw(self) -> None:
        self.ax.clear()
        box = ViewBox(self.half)
        draw_hypograph_scene(self.ax, self.vd, self.index, box)
        set_equal_3d(self.ax, lim=box.lim)
        self.ax.view_init(elev=self.elev, azim=self.azim)
        self.ax.set_title(f"{self._title()}\n{self._caption()}", fontsize=10)
        self.fig.subplots_adjust(bottom=0.08, top=0.88)
        self.fig.text(0.5, 0.015, _HELP, ha="center", va="bottom", fontsize=9)
        manager = getattr(self.fig.canvas, "manager", None)
        if manager is not None and hasattr(manager, "set_window_title"):
            manager.set_window_title(f"VD3D — {self._title()}")
        self.fig.canvas.draw_idle()

    def _on_key(self, event) -> None:
        key = event.key
        n = len(self.vd.cells)
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


def show_hypograph_vd(
    vd: AnyVD,
    *,
    seed: int | None = None,
    exact: bool = False,
    compare: HypographCompare | None = None,
) -> HypographGui:
    if not plt.isinteractive() and plt.get_backend().lower() == "agg":
        raise RuntimeError(
            "matplotlib backend is Agg (no window). "
            "Unset MPLBACKEND and install a GUI backend (e.g. TkAgg)."
        )
    gui = HypographGui(vd, seed=seed, exact=exact, compare=compare)
    gui.show()
    return gui
