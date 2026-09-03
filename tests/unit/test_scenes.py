"""Scenes must draw without opening a GUI (Agg / headless)."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from vd3d.viz.random_geom import choose_seed
from vd3d.viz.scenes import PHASE1_SCENES, make_phase1_scenes, scenes_for_phase


def test_phase1_fixtures_registered():
    scenes = scenes_for_phase(1, fixtures=True)
    assert [scene.name for scene in scenes] == [scene.name for scene in PHASE1_SCENES]
    assert len(scenes) == 5


def test_unknown_phase_rejected():
    try:
        scenes_for_phase(99, seed=1)
    except ValueError as exc:
        assert "99" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_each_fixture_scene_draws():
    for scene in PHASE1_SCENES:
        fig = plt.figure(figsize=scene.figsize)
        scene.draw(fig)
        assert fig.axes
        plt.close(fig)


def test_random_scenes_draw():
    for scene in make_phase1_scenes(1):
        fig = plt.figure(figsize=scene.figsize)
        scene.draw(fig)
        assert fig.axes
        plt.close(fig)


def test_same_seed_is_reproducible():
    first = make_phase1_scenes(12345)
    second = make_phase1_scenes(12345)
    assert [scene.caption for scene in first] == [scene.caption for scene in second]
    assert [scene.name for scene in first] == [scene.name for scene in second]


def test_choose_seed_explicit():
    assert choose_seed(7) == 7


def test_choose_seed_none_uses_time(monkeypatch):
    monkeypatch.setattr("vd3d.viz.random_geom.time.time_ns", lambda: 424242)
    assert choose_seed(None) == 424242
