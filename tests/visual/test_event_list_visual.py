"""Phase 7 snapshot figures.

    python -m vd3d.viz.gallery --step 7
    python -m vd3d.viz.viewer --phase 7
"""

import pytest
import matplotlib.pyplot as plt

from vd3d.viz.gallery import build_gallery
from vd3d.viz.record import record_figure
from vd3d.viz.scene import Scene
from vd3d.viz.scenes_event_list import PHASE7_SCENES

STEP = 7


@pytest.fixture(scope="module", autouse=True)
def _build_gallery():
    yield
    build_gallery(step=STEP)


def _record(fig, scene: Scene):
    png = record_figure(fig, step=STEP, name=scene.name, caption=scene.caption)
    plt.close(fig)
    assert png.is_file()
    assert png.stat().st_size > 0
    assert "What you must see" in png.with_suffix(".md").read_text(encoding="utf-8")
    return png


@pytest.mark.visual
@pytest.mark.parametrize("scene", PHASE7_SCENES, ids=lambda scene: scene.name)
def test_phase7_scene_snapshot(scene: Scene):
    fig = plt.figure(figsize=scene.figsize)
    scene.draw(fig)
    _record(fig, scene)
