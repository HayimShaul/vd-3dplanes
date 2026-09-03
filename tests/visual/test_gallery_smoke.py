"""Phase 0 smoke figure: prove the gallery pipeline writes a reviewable PNG."""

from pathlib import Path

import pytest
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

from vd3d.geometry.scalar import as_scalar
from vd3d.viz.convert import to_float
from vd3d.viz.gallery import build_gallery
from vd3d.viz.record import record_figure

STEP = 0
CAPTION = """\
**What you must see:** a unit square with corners (0,0), (1,0), (1,1), (0,1),
axis-aligned, sitting in the first quadrant. Corner dots are labeled.

**If this is wrong:** the gallery pipeline is broken. Do not proceed.
"""


@pytest.mark.visual
def test_unit_square_smoke_figure():
    origin_x = to_float(as_scalar(0))
    origin_y = to_float(as_scalar(0))
    width = to_float(as_scalar(1))
    height = to_float(as_scalar(1))

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.add_patch(
        Rectangle(
            (origin_x, origin_y),
            width,
            height,
            fill=False,
            linewidth=2,
            edgecolor="black",
        )
    )
    corners = [(0, 0), (1, 0), (1, 1), (0, 1)]
    ax.scatter([c[0] for c in corners], [c[1] for c in corners], c="black", zorder=3)
    for x, y in corners:
        ax.annotate(f"({x},{y})", (x, y), textcoords="offset points", xytext=(6, 6))
    ax.set_xlim(-0.25, 1.25)
    ax.set_ylim(-0.25, 1.25)
    ax.set_aspect("equal")
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.set_title("Phase 0 smoke: unit square")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    png = record_figure(fig, step=STEP, name="unit_square", caption=CAPTION)
    plt.close(fig)

    assert png.is_file()
    assert png.stat().st_size > 0
    caption_path = png.with_suffix(".md")
    assert caption_path.is_file()
    assert "What you must see" in caption_path.read_text(encoding="utf-8")

    index = build_gallery(step=STEP)
    assert index.is_file()
    html = index.read_text(encoding="utf-8")
    assert "How to review a step" in html
    assert "unit_square" in html
    assert Path(index).exists()
