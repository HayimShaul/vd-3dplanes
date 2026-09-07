"""Phase 8 snapshot figures.

    python -m vd3d.viz.gallery --step 8
    python -m vd3d.viz.viewer --phase 8
"""

import pytest
import matplotlib.pyplot as plt

from vd3d.events.samples import planes_through_123
from vd3d.sweep import vertical_decomposition_3d
from vd3d.viz.gallery import build_gallery
from vd3d.viz.paths import step_dir
from vd3d.viz.record import record_figure
from vd3d.viz.scene import Scene
from vd3d.viz.scenes_sweep import PHASE8_SCENES

STEP = 8


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
@pytest.mark.parametrize("scene", PHASE8_SCENES, ids=lambda scene: scene.name)
def test_phase8_scene_snapshot(scene: Scene):
    fig = plt.figure(figsize=scene.figsize)
    scene.draw(fig)
    _record(fig, scene)


@pytest.mark.visual
def test_phase8_snapshot_json():
    out = step_dir(STEP)
    result = vertical_decomposition_3d(planes_through_123(), snapshot_dir=out)
    assert result.snapshots
    json_files = list(out.glob("*.json"))
    assert json_files
    text = json_files[0].read_text(encoding="utf-8")
    assert "n_active" in text
    assert "TRIPLE_INTERSECTION" in text
