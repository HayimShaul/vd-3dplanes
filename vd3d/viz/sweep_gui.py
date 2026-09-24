"""Interactive before/after sweep-plane viewer.

    python -m vd3d examples/three_planes.txt --show-sweep
    python -m vd3d --random-planes 4 --seed 42 --show-sweep

For each event at ``z0``, shows the 2D VD at ``z0-α`` and ``z0+α`` side by
side (α small enough that no other event lies in ``[z0-α, z0+α]``).

Keys:
  n / p       next / previous event
  ← / →       move the viewing window along x
  ↑ / ↓       move the viewing window along y
  + / -       zoom in / out (shrink / grow the window)
  z / Z       increase / decrease α
  q / Esc     quit
"""

from __future__ import annotations

import matplotlib.pyplot as plt

from vd3d.events.slice import compute_vd_at_z
from vd3d.events.types import Event
from vd3d.geometry.scalar import Scalar, as_scalar
from vd3d.sweep.algorithm import SweepResult, group_events_by_z
from vd3d.sweep.around import z_before_after
from vd3d.sweep.matching import event_points_2d
from vd3d.sweep.types import EventGroup
from vd3d.viz.convert import to_float
from vd3d.viz.plot_sweep_plane import ViewWindow2D, draw_sweep_slice, event_highlights
from vd3d.viz.plot_vd import vd_viewing_limit

_HELP = "n/p event   ←→/↑↓ pan   +/- zoom   z/Z α   q quit"
_MIN_ALPHA = as_scalar("1/1000")
_DEFAULT_ALPHA = as_scalar("1/100")


def max_alpha_for_group(group: EventGroup, events: tuple[Event, ...]) -> Scalar:
    """Largest α such that ``[z0-α, z0+α]`` contains no other event ``z``."""
    z0 = group.z
    others = [event.z for event in events if event.z != z0]
    if not others:
        return as_scalar(10)
    gap = min(abs(other - z0) for other in others)
    # Strictly inside the closed interval constraint: keep a tiny margin.
    return gap * as_scalar("999/1000")


def resolve_alpha(requested: Scalar, group: EventGroup, events: tuple[Event, ...]) -> Scalar:
    cap = max_alpha_for_group(group, events)
    if requested < _MIN_ALPHA:
        return min(_MIN_ALPHA, cap)
    if requested > cap:
        return cap
    return requested


class SweepPlaneGui:
    """Step through events; each step shows z−α and z+α side by side."""

    def __init__(self, result: SweepResult, *, seed: int | None = None) -> None:
        self.result = result
        self.seed = seed
        self.groups = group_events_by_z(result.events)
        if not self.groups:
            raise ValueError("no sweep events to display")
        self.index = 0
        self.alpha = _DEFAULT_ALPHA
        z_minus, _z_plus = z_before_after(
            self.groups[0].representative, result.events, default_eps=self._alpha()
        )
        half = float(
            max(
                vd_viewing_limit(compute_vd_at_z(result.planes, z_minus), minimum=3.0),
                3.0,
            )
        )
        # Center near the first event feature if possible.
        pts = event_points_2d(self.groups[0].representative)
        cx = to_float(pts[0].x) if pts else 0.0
        cy = to_float(pts[0].y) if pts else 0.0
        self.window = ViewWindow2D(cx=cx, cy=cy, half=half)
        self.fig = plt.figure(figsize=(12.0, 6.2))
        self.ax_before = self.fig.add_subplot(121)
        self.ax_after = self.fig.add_subplot(122)
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self._draw()

    def _group(self) -> EventGroup:
        return self.groups[self.index]

    def _alpha(self) -> Scalar:
        return resolve_alpha(self.alpha, self._group(), self.result.events)

    def _pan_step(self) -> float:
        return 0.15 * self.window.half

    def _title_bar(self) -> str:
        group = self._group()
        alpha = self._alpha()
        kinds = ",".join(sorted({event.type.name for event in group.events}))
        parts = [
            f"event {self.index + 1}/{len(self.groups)}",
            f"z={group.z}",
            f"α={alpha}",
            kinds,
        ]
        if self.seed is not None:
            parts.append(f"seed={self.seed}")
        return "  ·  ".join(parts)

    def _event_points(self):
        points = []
        for event in self._group().events:
            points.extend(event_points_2d(event))
        return tuple(points)

    def _draw(self) -> None:
        group = self._group()
        alpha = self._alpha()
        self.alpha = alpha
        z0 = group.z
        z_before = z0 - alpha
        z_after = z0 + alpha
        vd_before = compute_vd_at_z(self.result.planes, z_before)
        vd_after = compute_vd_at_z(self.result.planes, z_after)
        points = self._event_points()
        planes = self.result.planes
        hot_before = event_highlights(vd_before, group.events, planes)
        hot_after = event_highlights(vd_after, group.events, planes)
        draw_sweep_slice(
            self.ax_before,
            vd_before,
            window=self.window,
            z=z_before,
            title=f"before  z={z_before}",
            result=self.result,
            event_points=points,
            highlight_plane_ids=hot_before[0],
            highlight_wall_xs=hot_before[1],
        )
        draw_sweep_slice(
            self.ax_after,
            vd_after,
            window=self.window,
            z=z_after,
            title=f"after  z={z_after}",
            result=self.result,
            event_points=points,
            highlight_plane_ids=hot_after[0],
            highlight_wall_xs=hot_after[1],
        )
        self.fig.suptitle(self._title_bar(), fontsize=11)
        self.fig.subplots_adjust(bottom=0.10, top=0.88, wspace=0.18)
        # Clear previous help text artists by removing fig texts except via clear+redraw is hard;
        # remove old help texts tagged with _sweep_help.
        for artist in list(self.fig.texts):
            if getattr(artist, "_sweep_help", False):
                artist.remove()
        help_artist = self.fig.text(
            0.5, 0.02, _HELP, ha="center", va="bottom", fontsize=9
        )
        help_artist._sweep_help = True  # type: ignore[attr-defined]
        manager = getattr(self.fig.canvas, "manager", None)
        if manager is not None and hasattr(manager, "set_window_title"):
            manager.set_window_title(f"VD3D sweep — {self._title_bar()}")
        self.fig.canvas.draw_idle()

    def _on_key(self, event) -> None:
        key = event.key
        n = len(self.groups)
        step = self._pan_step()
        if key in {"n", " ", "pagedown"}:
            self.index = (self.index + 1) % n
            self._center_on_event()
            self._draw()
        elif key in {"p", "pageup"}:
            self.index = (self.index - 1) % n
            self._center_on_event()
            self._draw()
        elif key == "left":
            self.window.cx -= step
            self._draw()
        elif key == "right":
            self.window.cx += step
            self._draw()
        elif key == "up":
            self.window.cy += step
            self._draw()
        elif key == "down":
            self.window.cy -= step
            self._draw()
        elif key in {"+", "=", "plus"}:
            self.window.half = max(0.4, self.window.half * 0.8)
            self._draw()
        elif key in {"-", "minus", "_"}:
            self.window.half = min(200.0, self.window.half * 1.25)
            self._draw()
        elif key == "z":
            self.alpha = resolve_alpha(self.alpha * as_scalar("5/4"), self._group(), self.result.events)
            self._draw()
        elif key == "Z":
            self.alpha = resolve_alpha(self.alpha * as_scalar("4/5"), self._group(), self.result.events)
            self._draw()
        elif key in {"q", "escape"}:
            plt.close(self.fig)

    def _center_on_event(self) -> None:
        pts = self._event_points()
        if not pts:
            return
        self.window.cx = sum(to_float(p.x) for p in pts) / len(pts)
        self.window.cy = sum(to_float(p.y) for p in pts) / len(pts)

    def show(self) -> None:
        plt.show()


def show_sweep(result: SweepResult, *, seed: int | None = None) -> SweepPlaneGui:
    """Open the sweep-plane GUI. Requires a non-Agg matplotlib backend."""
    if not plt.isinteractive() and plt.get_backend().lower() == "agg":
        raise RuntimeError(
            "matplotlib backend is Agg (no window). "
            "Unset MPLBACKEND and install a GUI backend (e.g. TkAgg)."
        )
    if not result.events:
        raise ValueError("no sweep events to display")
    gui = SweepPlaneGui(result, seed=seed)
    gui.show()
    return gui
