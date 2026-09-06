"""Scenes must draw without opening a GUI (Agg / headless)."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from vd3d.viz.random_geom import choose_seed
from vd3d.viz.scenes import PHASE1_SCENES, make_phase1_scenes, scenes_for_phase
from vd3d.viz.scenes_arrangement import PHASE2_SCENES, make_phase2_scenes
from vd3d.viz.scenes_vd import PHASE3_SCENES, make_phase3_scenes
from vd3d.viz.scenes_zone import PHASE4_SCENES, make_phase4_scenes
from vd3d.viz.scenes_events import PHASE5_SCENES, make_phase5_scenes
from vd3d.viz.scenes_alignment import PHASE6_SCENES, make_phase6_scenes
from vd3d.viz.scenes_event_list import PHASE7_SCENES, make_phase7_scenes


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


def test_phase4_fixtures_registered():
    scenes = scenes_for_phase(4, fixtures=True)
    assert [scene.name for scene in scenes] == [scene.name for scene in PHASE4_SCENES]
    assert len(scenes) == 7


def test_phase5_fixtures_registered():
    scenes = scenes_for_phase(5, fixtures=True)
    assert [scene.name for scene in scenes] == [scene.name for scene in PHASE5_SCENES]
    assert len(scenes) == 3


def test_phase6_fixtures_registered():
    scenes = scenes_for_phase(6, fixtures=True)
    assert [scene.name for scene in scenes] == [scene.name for scene in PHASE6_SCENES]
    assert len(scenes) == 5


def test_phase7_fixtures_registered():
    scenes = scenes_for_phase(7, fixtures=True)
    assert [scene.name for scene in scenes] == [scene.name for scene in PHASE7_SCENES]
    assert len(scenes) == 2


def test_unknown_phase_rejected():
    try:
        scenes_for_phase(99, seed=1)
    except ValueError as exc:
        assert "99" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_each_fixture_scene_draws():
    for scene in (
        *PHASE1_SCENES,
        *PHASE2_SCENES,
        *PHASE3_SCENES,
        *PHASE4_SCENES,
        *PHASE5_SCENES,
        *PHASE6_SCENES,
        *PHASE7_SCENES,
    ):
        fig = plt.figure(figsize=scene.figsize)
        scene.draw(fig)
        assert fig.axes
        plt.close(fig)


def test_random_scenes_draw():
    for scene in (
        *make_phase1_scenes(1),
        *make_phase2_scenes(1),
        *make_phase3_scenes(1),
        *make_phase4_scenes(1),
        *make_phase5_scenes(1),
        *make_phase6_scenes(1),
        *make_phase7_scenes(1),
    ):
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
    first4 = make_phase4_scenes(12345)
    second4 = make_phase4_scenes(12345)
    assert [scene.caption for scene in first4] == [scene.caption for scene in second4]
    first5 = make_phase5_scenes(12345)
    second5 = make_phase5_scenes(12345)
    assert [scene.caption for scene in first5] == [scene.caption for scene in second5]
    first6 = make_phase6_scenes(12345)
    second6 = make_phase6_scenes(12345)
    assert [scene.caption for scene in first6] == [scene.caption for scene in second6]
    first7 = make_phase7_scenes(12345)
    second7 = make_phase7_scenes(12345)
    assert [scene.caption for scene in first7] == [scene.caption for scene in second7]
    assert [scene.name for scene in first7] == [scene.name for scene in second7]


def test_choose_seed_explicit():
    assert choose_seed(7) == 7


def test_choose_seed_none_uses_time(monkeypatch):
    monkeypatch.setattr("vd3d.viz.random_geom.time.time_ns", lambda: 424242)
    assert choose_seed(None) == 424242


def test_n_controls_line_count_phase2_to_4():
    for scene in make_phase2_scenes(0, n=6):
        if scene.name in {"pieces_random", "faces_random"}:
            assert "n=6" in scene.caption
    for scene in make_phase3_scenes(0, n=6):
        assert "n=6" in scene.caption
    for scene in make_phase4_scenes(0, n=6):
        assert "n=6" in scene.caption
    for scene in make_phase5_scenes(0, n=4):
        assert "n=4" in scene.caption
    for scene in make_phase6_scenes(0, n=4):
        assert "n=4" in scene.caption
    for scene in make_phase7_scenes(0, n=4):
        assert "n=4" in scene.caption
    six = make_phase4_scenes(0, n=6)
    seven = make_phase4_scenes(0, n=7)
    assert [s.caption for s in six] != [s.caption for s in seven]
    four = make_phase5_scenes(0, n=4)
    three = make_phase5_scenes(0, n=3)
    assert [s.caption for s in four] != [s.caption for s in three]
    four6 = make_phase6_scenes(0, n=4)
    five6 = make_phase6_scenes(0, n=5)
    assert [s.caption for s in four6] != [s.caption for s in five6]
    four7 = make_phase7_scenes(0, n=4)
    three7 = make_phase7_scenes(0, n=3)
    assert [s.caption for s in four7] != [s.caption for s in three7]


def test_n_must_be_positive():
    try:
        make_phase4_scenes(0, n=0)
    except ValueError as exc:
        assert "positive" in str(exc)
    else:
        raise AssertionError("expected ValueError")
    try:
        make_phase5_scenes(0, n=0)
    except ValueError as exc:
        assert "positive" in str(exc)
    else:
        raise AssertionError("expected ValueError")
    try:
        make_phase6_scenes(0, n=0)
    except ValueError as exc:
        assert "positive" in str(exc)
    else:
        raise AssertionError("expected ValueError")
    try:
        make_phase7_scenes(0, n=0)
    except ValueError as exc:
        assert "positive" in str(exc)
    else:
        raise AssertionError("expected ValueError")


def test_scenes_for_phase_forwards_n():
    scenes = scenes_for_phase(4, seed=1, n=5)
    assert all("n=5" in scene.caption for scene in scenes)
    scenes5 = scenes_for_phase(5, seed=1, n=4)
    assert all("n=4" in scene.caption for scene in scenes5)
    scenes6 = scenes_for_phase(6, seed=1, n=4)
    assert all("n=4" in scene.caption for scene in scenes6)
    scenes7 = scenes_for_phase(7, seed=1, n=4)
    assert all("n=4" in scene.caption for scene in scenes7)
