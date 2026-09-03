"""Viewer draws and cycles scenes without opening a window."""

from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from vd3d.viz.scenes import PHASE1_SCENES, make_phase1_scenes
from vd3d.viz.viewer import ReviewViewer, main


def test_viewer_starts_on_requested_scene():
    viewer = ReviewViewer(PHASE1_SCENES, start=2)
    assert viewer.index == 2
    assert viewer.fig.axes
    plt.close(viewer.fig)


def test_viewer_cycles_with_keys():
    viewer = ReviewViewer(PHASE1_SCENES, start=0)
    viewer._on_key(SimpleNamespace(key="right"))
    assert viewer.index == 1
    viewer._on_key(SimpleNamespace(key="left"))
    assert viewer.index == 0
    viewer._on_key(SimpleNamespace(key="p"))
    assert viewer.index == len(PHASE1_SCENES) - 1
    plt.close(viewer.fig)


def test_viewer_reroll_uses_new_seed(monkeypatch):
    monkeypatch.setattr("vd3d.viz.viewer.choose_seed", lambda seed: 99 if seed is None else seed)
    scenes = make_phase1_scenes(1)
    viewer = ReviewViewer(scenes, start=2, phase=1, seed=1)
    viewer._on_key(SimpleNamespace(key="g"))
    assert viewer.seed == 99
    assert viewer.index == 0
    assert viewer.scenes[0].caption.startswith("seed=99")
    plt.close(viewer.fig)


def test_main_refuses_agg_backend():
    assert main(["--phase", "1"]) == 2
    assert main(["--phase", "1", "--seed", "1"]) == 2
