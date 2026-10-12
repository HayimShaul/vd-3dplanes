"""Interactive 3D viewer for a lower envelope (float and/or exact).

    python -m vd3d.lower_envelope examples/three_planes.txt --gui
    python -m vd3d.lower_envelope examples/three_planes.txt --exact --gui

Keys:
  + / -       zoom in / out
  ↑ / ↓       rotate elevation
  ← / →       rotate azimuth
  q / Esc     quit
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from vd3d.io import format_plane
from vd3d.lower_envelope.envelope import EnvelopeCompare, LowerEnvelope
from vd3d.lower_envelope.envelope_float import FloatLowerEnvelope
from vd3d.viz.plot_geometry import set_equal_3d
from vd3d.viz.plot_lower_envelope import (
    draw_lower_envelope_scene,
    envelope_view_half,
)
from vd3d.viz.plot_cells3d import ViewBox

_HELP = "+/- zoom   ↑↓ elev   ←→ azim   q quit"


class LowerEnvelopeGui:
    """Grey planes, blue envelope; optional orange float-vertex overlay."""

    def __init__(
        self,
        env: FloatLowerEnvelope | LowerEnvelope,
        *,
        seed: int | None = None,
        float_overlay: FloatLowerEnvelope | None = None,
        compare: EnvelopeCompare | None = None,
        exact: bool = False,
    ) -> None:
        self.env = env
        self.seed = seed
        self.float_overlay = float_overlay
        self.compare = compare
        self.exact = exact
        self.half = envelope_view_half(env)
        self.elev = 22.0
        self.azim = -60.0
        self.fig = plt.figure(figsize=(9.5, 7.5))
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self._draw()

    def _title(self) -> str:
        kind = "exact" if self.exact else "float"
        n_active = len(self.env.active)
        n_planes = len(self.env.planes)
        n_facets = len(self.env.facets)
        parts = [
            f"{kind} lower envelope  active={n_active}/{n_planes}",
            f"facets={n_facets}",
            f"box=±{self.half:g}",
        ]
        if self.compare is not None:
            parts.append("compare=" + ("OK" if self.compare.ok else "MISMATCH"))
            if self.compare.max_vertex_error is not None:
                parts.append(f"err={self.compare.max_vertex_error:.2e}")
        if self.seed is not None:
            parts.append(f"seed={self.seed}")
        return "  ·  ".join(parts)

    def _caption(self) -> str:
        bits = []
        if self.float_overlay is not None:
            bits.append("blue=exact verts  orange×=float verts")
        if isinstance(self.env, LowerEnvelope):
            names = [
                f"P{self.env.planes[i].id}:{format_plane(self.env.planes[i])}"
                for i in self.env.active
            ]
        else:
            names = [
                f"P{self.env.planes[int(i)].id}:{format_plane(self.env.planes[int(i)])}"
                for i in self.env.active
            ]
        active = "active: " + ("; ".join(names) if names else "(none)")
        if bits:
            return bits[0] + "\n" + active
        return active

    def _draw(self) -> None:
        self.ax.clear()
        box = ViewBox(self.half)
        draw_lower_envelope_scene(
            self.ax, self.env, box, float_overlay=self.float_overlay
        )
        set_equal_3d(self.ax, lim=box.lim)
        self.ax.view_init(elev=self.elev, azim=self.azim)
        self.ax.set_title(f"{self._title()}\n{self._caption()}", fontsize=10)
        self.fig.subplots_adjust(bottom=0.08, top=0.86)
        self.fig.text(0.5, 0.015, _HELP, ha="center", va="bottom", fontsize=9)
        manager = getattr(self.fig.canvas, "manager", None)
        if manager is not None and hasattr(manager, "set_window_title"):
            manager.set_window_title(f"VD3D — {self._title()}")
        self.fig.canvas.draw_idle()

    def _on_key(self, event) -> None:
        key = event.key
        if key in {"+", "=", "plus"}:
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


def show_lower_envelope(
    env: FloatLowerEnvelope | LowerEnvelope,
    *,
    seed: int | None = None,
    float_overlay: FloatLowerEnvelope | None = None,
    compare: EnvelopeCompare | None = None,
    exact: bool = False,
) -> LowerEnvelopeGui:
    """Open the GUI. Requires a non-Agg matplotlib backend."""
    if not plt.isinteractive() and plt.get_backend().lower() == "agg":
        raise RuntimeError(
            "matplotlib backend is Agg (no window). "
            "Unset MPLBACKEND and install a GUI backend (e.g. TkAgg)."
        )
    gui = LowerEnvelopeGui(
        env,
        seed=seed,
        float_overlay=float_overlay,
        compare=compare,
        exact=exact,
    )
    gui.show()
    return gui
