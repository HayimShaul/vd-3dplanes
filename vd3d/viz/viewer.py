"""Interactive review window for 3D geometry scenes.

Phase 1 draws random planes, lines, and points. Pass ``--seed`` to replay
a run; with no seed the current time is used and printed::

    python -m vd3d.viz.viewer --phase 1
    python -m vd3d.viz.viewer --phase 1 --seed 42

Keys: left/right or n/p = next/previous scene, g = new random seed,
r = reset camera, q = quit. Drag a 3D view to rotate it.
"""

from __future__ import annotations

import argparse
import sys

import matplotlib.pyplot as plt

from vd3d.viz.random_geom import choose_seed
from vd3d.viz.scene import Scene
from vd3d.viz.scenes import scenes_for_phase

_HELP = "drag 3D to rotate   ←/→ or n/p change scene   g new seed   r reset view   q quit"


class ReviewViewer:
    def __init__(
        self,
        scenes: tuple[Scene, ...],
        *,
        start: int = 0,
        phase: int | None = None,
        seed: int | None = None,
    ) -> None:
        if not scenes:
            raise ValueError("no scenes to show")
        self.phase = phase
        self.seed = seed
        self.scenes = scenes
        self.index = start % len(scenes)
        self.fig = plt.figure(figsize=scenes[self.index].figsize)
        self._set_window_title("VD3D review")
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self._draw_current()

    def _seed_label(self) -> str:
        if self.seed is None:
            return ""
        return f"  seed={self.seed}"

    def _draw_current(self) -> None:
        scene = self.scenes[self.index]
        self.fig.clear()
        self.fig.set_size_inches(*scene.figsize, forward=True)
        scene.draw(self.fig)
        self.fig.subplots_adjust(bottom=0.22, top=0.90)
        self.fig.text(
            0.5,
            0.02,
            f"{scene.name}  ({self.index + 1}/{len(self.scenes)}){self._seed_label()}\n"
            f"{scene.caption}\n\n{_HELP}",
            ha="center",
            va="bottom",
            fontsize=9,
            wrap=True,
        )
        self._set_window_title(
            f"VD3D review — {scene.title} ({self.index + 1}/{len(self.scenes)})"
            f"{self._seed_label()}"
        )
        self.fig.canvas.draw_idle()

    def _set_window_title(self, title: str) -> None:
        manager = getattr(self.fig.canvas, "manager", None)
        if manager is not None and hasattr(manager, "set_window_title"):
            manager.set_window_title(title)

    def _reroll(self) -> None:
        if self.phase is None:
            return
        self.seed = choose_seed(None)
        print(f"seed={self.seed}", flush=True)
        self.scenes = scenes_for_phase(self.phase, seed=self.seed)
        self.index = 0
        self._draw_current()

    def _on_key(self, event) -> None:
        key = event.key
        if key in {"right", "n", " "}:
            self.index = (self.index + 1) % len(self.scenes)
            self._draw_current()
        elif key in {"left", "p"}:
            self.index = (self.index - 1) % len(self.scenes)
            self._draw_current()
        elif key == "g":
            self._reroll()
        elif key == "r":
            self._draw_current()
        elif key in {"q", "escape"}:
            plt.close(self.fig)

    def show(self) -> None:
        plt.show()


def show_phase(phase: int, *, scene: str | None = None, seed: int | None = None) -> ReviewViewer:
    resolved = choose_seed(seed)
    print(f"seed={resolved}", flush=True)
    scenes = scenes_for_phase(phase, seed=resolved)
    start = 0
    if scene is not None:
        names = [item.name for item in scenes]
        try:
            start = names.index(scene)
        except ValueError as exc:
            raise SystemExit(f"unknown scene {scene!r}; choose from {names}") from exc
    viewer = ReviewViewer(scenes, start=start, phase=phase, seed=resolved)
    viewer.show()
    return viewer


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Review geometry and arrangement scenes")
    parser.add_argument("--phase", type=int, default=1, help="phase number (1 or 2)")
    parser.add_argument("--scene", type=str, default=None, help="start at this scene name")
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="RNG seed (default: current time in nanoseconds)",
    )
    args = parser.parse_args(argv)
    if not plt.isinteractive() and plt.get_backend().lower() == "agg":
        print(
            "matplotlib backend is Agg (no window). "
            "Unset MPLBACKEND and run outside pytest.",
            file=sys.stderr,
        )
        return 2
    show_phase(args.phase, scene=args.scene, seed=args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
