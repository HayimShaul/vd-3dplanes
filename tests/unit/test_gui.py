"""Clipping and drawing helpers for the 3D VD GUI."""

from __future__ import annotations

import numpy as np

from vd3d.events.samples import planes_through_123
from vd3d.geometry.line3d import Line3D
from vd3d.geometry.plane import Plane
from vd3d.geometry.points import Point3D
from vd3d.sweep import vertical_decomposition_3d
from vd3d.viz.gui import DecompositionGui
from vd3d.viz.plot_cells3d import (
    ViewBox,
    clip_line_to_box,
    clipped_cell_faces,
    default_view_half,
    interior_point_of_cell,
    plane_polygon_in_box,
)
import matplotlib.pyplot as plt


def test_plane_polygon_in_box_covers_unit_plane():
    box = ViewBox(2.0)
    plane = Plane(id=1, a=0, b=0, c=1, d=0)  # z = 0
    poly = plane_polygon_in_box(plane, box)
    assert poly is not None
    assert len(poly) >= 3
    assert np.allclose(poly[:, 2], 0.0, atol=1e-8)


def test_clip_line_to_box_segment():
    box = ViewBox(1.0)
    line = Line3D(Point3D(0, 0, 0), (1, 0, 0))
    clipped = clip_line_to_box(line, box)
    assert clipped is not None
    a, b = clipped
    assert np.allclose(sorted([a[0], b[0]]), [-1.0, 1.0])
    assert np.allclose(a[1:], 0.0) and np.allclose(b[1:], 0.0)


def test_clipped_cell_faces_nonempty_for_triple():
    result = vertical_decomposition_3d(planes_through_123())
    box = ViewBox(default_view_half(result))
    nonempty = 0
    for cell in result.cells:
        assert interior_point_of_cell(result, cell.id) is not None
        faces = clipped_cell_faces(result, cell, box)
        if faces:
            nonempty += 1
            assert all(len(face) >= 3 for face in faces)
    assert nonempty >= 1


def test_gui_draws_without_show():
    result = vertical_decomposition_3d(planes_through_123())
    gui = DecompositionGui(result, seed=7)
    assert gui.index == 0
    assert gui.half > 0
    gui._on_key(type("E", (), {"key": "n"})())
    assert gui.index == 1
    gui._on_key(type("E", (), {"key": "p"})())
    assert gui.index == 0
    before = gui.half
    gui._on_key(type("E", (), {"key": "+"})())
    assert gui.half < before
    shrunk = gui.half
    gui._on_key(type("E", (), {"key": "-"})())
    assert gui.half > shrunk
    elev = gui.elev
    gui._on_key(type("E", (), {"key": "up"})())
    assert gui.elev != elev
    plt.close(gui.fig)
