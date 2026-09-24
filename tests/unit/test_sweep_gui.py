"""Sweep-plane GUI helpers."""

from __future__ import annotations

import matplotlib.pyplot as plt

from vd3d.events.samples import planes_through_123
from vd3d.events.slice import compute_vd_at_z
from vd3d.geometry.scalar import as_scalar
from vd3d.sweep import vertical_decomposition_3d
from vd3d.sweep.around import z_before_after
from vd3d.viz.plot_sweep_plane import (
    ViewWindow2D,
    adjacent_cell_colors,
    cell_clip_polygon_window,
    clip_line_to_window,
    event_highlights,
)
from vd3d.viz.sweep_gui import SweepPlaneGui, max_alpha_for_group, resolve_alpha


def test_adjacent_colors_differ_for_neighbors():
    result = vertical_decomposition_3d(planes_through_123())
    event = result.events[0]
    z_minus, _ = z_before_after(event, result.events)
    vd = compute_vd_at_z(result.planes, z_minus)
    colors = adjacent_cell_colors(vd)
    for cell in vd.cells:
        for nid in cell.neighbors:
            assert colors[cell.id] != colors[nid]


def test_clip_line_and_cell_to_window():
    result = vertical_decomposition_3d(planes_through_123())
    event = result.events[0]
    z_minus, _ = z_before_after(event, result.events)
    vd = compute_vd_at_z(result.planes, z_minus)
    window = ViewWindow2D(cx=1.0, cy=2.0, half=3.0)
    assert any(len(clip_line_to_window(line, window)) == 2 for line in vd.arrangement.lines)
    polys = [cell_clip_polygon_window(vd, cell, window) for cell in vd.cells]
    assert any(len(poly) >= 3 for poly in polys)


def test_alpha_clamped_away_from_other_events():
    result = vertical_decomposition_3d(planes_through_123())
    from vd3d.sweep import group_events_by_z

    groups = group_events_by_z(result.events)
    cap = max_alpha_for_group(groups[0], result.events)
    assert resolve_alpha(as_scalar(100), groups[0], result.events) == cap
    assert resolve_alpha(as_scalar("1/100"), groups[0], result.events) == as_scalar("1/100")


def test_event_highlights_triple_planes_and_alignment_walls():
    from vd3d.events.samples import planes_alignment_at_z2
    from vd3d.events.types import EventType

    triple = vertical_decomposition_3d(planes_through_123())
    event = next(e for e in triple.events if e.type is EventType.TRIPLE_INTERSECTION)
    z_minus, _ = z_before_after(event, triple.events)
    vd = compute_vd_at_z(triple.planes, z_minus)
    planes_hot, walls_hot = event_highlights(vd, (event,), triple.planes)
    assert planes_hot == frozenset(event.plane_ids)
    assert walls_hot == frozenset()

    align_result = vertical_decomposition_3d(planes_alignment_at_z2())
    align = next(e for e in align_result.events if e.type is EventType.VERTICAL_ALIGNMENT)
    # Isolate so z± does not land on a simultaneous triple.
    z_minus, _ = z_before_after(align, (align,))
    vd = compute_vd_at_z(align_result.planes, z_minus)
    planes_hot, walls_hot = event_highlights(vd, (align,), align_result.planes)
    assert planes_hot == frozenset()
    assert len(walls_hot) == 2


def test_sweep_gui_keys():
    result = vertical_decomposition_3d(planes_through_123())
    gui = SweepPlaneGui(result, seed=3)
    assert gui.index == 0
    before = gui.window.cx
    gui._on_key(type("E", (), {"key": "right"})())
    assert gui.window.cx > before
    half = gui.window.half
    gui._on_key(type("E", (), {"key": "+"})())
    assert gui.window.half < half
    alpha = gui.alpha
    gui._on_key(type("E", (), {"key": "z"})())
    assert gui.alpha >= alpha
    gui._on_key(type("E", (), {"key": "Z"})())
    plt.close(gui.fig)
