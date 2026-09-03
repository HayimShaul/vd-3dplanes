"""Scenes must draw without opening a GUI (Agg / headless)."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from vd3d.viz.random_geom import choose_seed
from vd3d.viz.scenes import PHASE1_SCENES, make_phase1_scenes, scenes_for_phase
from vd3d.viz.scenes_arrangement import PHASE2_SCENES, make_phase2_scenes
from vd3d.viz.scenes_vd import PHASE3_SCENES, make_phase3_scenes


def test_phase1_fixtures_registered():
    scenes = scenes_for_phase(1, fixtures=True)
    assert [scene.name for scene in scenes] == [scene.name for scene in PHASE1_SCENES]
    assert len(scenes) == 5


def test_phase2_fixtures_registered():
    scenes = scenes_for_phase(2, fixtures=True)
    assert [scene.name for scene in scenes] == [scene.name for scene in PHASE2_SCENES]
    assert len(scenes) == 6


def test_phase3_fixtures_registered():
    scenes = scenes_for_phase(3, fixtures=True)
    assert [scene.name for scene in scenes] == [scene.name for scene in PHASE3_SCENES]
    assert len(scenes) == 6


def test_unknown_phase_rejected():
    try:
        scenes_for_phase(99, seed=1)
    except ValueError as exc:
        assert "99" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_each_fixture_scene_draws():
    for scene in (*PHASE1_SCENES, *PHASE2_SCENES, *PHASE3_SCENES):
        fig = plt.figure(figsize=scene.figsize)
        scene.draw(fig)
        assert fig.axes
        plt.close(fig)


def test_random_scenes_draw():
    for scene in (*make_phase1_scenes(1), *make_phase2_scenes(1), *make_phase3_scenes(1)):
        fig = plt.figure(figsize=scene.figsize)
        scene.draw(fig)
        assert fig.axes
        plt.close(fig)


def test_same_seed_is_reproducible():
    first = make_phase1_scenes(12345)
    second = make_phase1_scenes(12345)
    assert [scene.caption for scene in first] == [scene.caption for scene in second]
    assert [scene.name for scene in first] == [scene.name for scene in second]
    first2 = make_phase2_scenes(12345)
    second2 = make_phase2_scenes(12345)
    assert [scene.caption for scene in first2] == [scene.caption for scene in second2]
    first3 = make_phase3_scenes(12345)
    second3 = make_phase3_scenes(12345)
    assert [scene.caption for scene in first3] == [scene.caption for scene in second3]


def test_choose_seed_explicit():
    assert choose_seed(7) == 7


def test_choose_seed_none_uses_time(monkeypatch):
    monkeypatch.setattr("vd3d.viz.random_geom.time.time_ns", lambda: 424242)
    assert choose_seed(None) == 424242
