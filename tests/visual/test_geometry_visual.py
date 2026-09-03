"""Phase 1 snapshot figures. Human review is the rotatable GUI, not the HTML page.

    python -m vd3d.viz.viewer --phase 1
"""

import pytest
import matplotlib.pyplot as plt

from vd3d.viz.gallery import build_gallery
from vd3d.viz.record import record_figure
from vd3d.viz.scenes import PHASE1_SCENES, Scene

STEP = 1


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
@pytest.mark.parametrize("scene", PHASE1_SCENES, ids=lambda scene: scene.name)
def test_phase1_scene_snapshot(scene: Scene):
    fig = plt.figure(figsize=scene.figsize)
    scene.draw(fig)
    _record(fig, scene)
