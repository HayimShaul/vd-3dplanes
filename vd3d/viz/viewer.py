"""Interactive review window for geometry, arrangement, VD, zone, and events.

    python -m vd3d.viz.viewer --phase 1
    python -m vd3d.viz.viewer --phase 2 --n 6
    python -m vd3d.viz.viewer --phase 3 --n 6
    python -m vd3d.viz.viewer --phase 4 --n 8
    python -m vd3d.viz.viewer --phase 5
    python -m vd3d.viz.viewer --phase 6
    python -m vd3d.viz.viewer --phase 6 --seed 42 --n 4
    python -m vd3d.viz.viewer --phase 7
    python -m vd3d.viz.viewer --phase 7 --seed 42 --n 4

``--n`` is the number of lines (phases 2–4) or planes (phases 1, 5, 6, and 7).
Omit it to let each random scene pick a small count.

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
        n: int | None = None,
    ) -> None:
        if not scenes:
            raise ValueError("no scenes to show")
        self.phase = phase
        self.seed = seed
        self.n = n
        self.scenes = scenes
        self.index = start % len(scenes)
        self.fig = plt.figure(figsize=scenes[self.index].figsize)
        self._set_window_title("VD3D review")
        self.fig.canvas.mpl_connect("key_press_event", self._on_key)
        self._draw_current()

    def _seed_label(self) -> str:
        parts: list[str] = []
        if self.seed is not None:
            parts.append(f"seed={self.seed}")
        if self.n is not None:
            parts.append(f"n={self.n}")
        if not parts:
            return ""
        return "  " + " ".join(parts)

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
        if self.n is not None:
            print(f"n={self.n}", flush=True)
        self.scenes = scenes_for_phase(self.phase, seed=self.seed, n=self.n)
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


def show_phase(
    phase: int,
    *,
    scene: str | None = None,
    seed: int | None = None,
    n: int | None = None,
) -> ReviewViewer:
    resolved = choose_seed(seed)
    print(f"seed={resolved}", flush=True)
    if n is not None:
        print(f"n={n}", flush=True)
    scenes = scenes_for_phase(phase, seed=resolved, n=n)
    start = 0
    if scene is not None:
        names = [item.name for item in scenes]
        try:
            start = names.index(scene)
        except ValueError as exc:
            raise SystemExit(f"unknown scene {scene!r}; choose from {names}") from exc
    viewer = ReviewViewer(scenes, start=start, phase=phase, seed=resolved, n=n)
    viewer.show()
    return viewer


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Review geometry and arrangement scenes")
    parser.add_argument("--phase", type=int, default=1, help="phase number (1–7)")
    parser.add_argument("--scene", type=str, default=None, help="start at this scene name")
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="RNG seed (default: current time in nanoseconds)",
    )
    parser.add_argument(
        "-n",
        "--n",
        type=int,
        default=None,
        metavar="N",
        dest="n",
        help=(
            "number of lines (phases 2–4) or planes (phases 1, 5, 6, and 7). "
            "Default: each random scene picks a small count."
        ),
    )
    args = parser.parse_args(argv)
    if args.n is not None and args.n < 1:
        print("n must be a positive integer", file=sys.stderr)
        return 2
    if not plt.isinteractive() and plt.get_backend().lower() == "agg":
        print(
            "matplotlib backend is Agg (no window). "
            "Unset MPLBACKEND and run outside pytest.",
            file=sys.stderr,
        )
        return 2
    show_phase(args.phase, scene=args.scene, seed=args.seed, n=args.n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
