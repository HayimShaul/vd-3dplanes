"""Phase 1 geometry kernel: design tests 1–4 and the geometry invariants."""

from fractions import Fraction

import pytest

from vd3d.geometry import (
    PARALLEL,
    Line2D,
    Plane,
    Point2D,
    Point3D,
    build_slice_lines,
    directions_parallel,
    intersect_planes,
    intersect_three_planes,
    line_lies_on_both_planes,
    normals_parallel,
    point_lies_on_planes,
    slice_lifts_to_plane,
    slice_plane_at_z,
)


# ---------------------------------------------------------------------------
# Step 1.1 — points and planes
# ---------------------------------------------------------------------------


def test_eval_on_and_off_plane():
    plane = Plane(id=1, a=1, b=1, c=1, d=-5)
    on = Point3D(1, 2, 2)
    above = Point3D(1, 2, 3)
    below = Point3D(1, 2, 1)
    assert plane.eval(on) == 0
    assert plane.contains(on)
    assert plane.eval(above) == 1
    assert plane.eval(below) == -1
    assert not plane.contains(above)
    assert not plane.contains(below)


def test_normals_parallel_and_not():
    x0 = Plane(id=1, a=1, b=0, c=0, d=0)
    x2 = Plane(id=2, a=1, b=0, c=0, d=-2)
    y0 = Plane(id=3, a=0, b=1, c=0, d=0)
    flipped = Plane(id=4, a=-1, b=0, c=0, d=0)
    assert normals_parallel(x0, x2)
    assert normals_parallel(x0, flipped)
    assert not normals_parallel(x0, y0)


def test_zero_plane_rejected():
    with pytest.raises(ValueError, match="zero plane"):
        Plane(id=1, a=0, b=0, c=0, d=1)
    with pytest.raises(ValueError, match="zero plane"):
        Plane(id=1, a=0, b=0, c=0, d=0)


def test_point_and_plane_reject_float():
    with pytest.raises(TypeError):
        Point3D(0.5, 0, 0)
    with pytest.raises(TypeError):
        Point2D(0.5, 0)
    with pytest.raises(TypeError):
        Plane(id=1, a=1.0, b=0, c=0, d=0)


# ---------------------------------------------------------------------------
# Step 1.2 — plane–plane intersection
# ---------------------------------------------------------------------------


def test_intersect_planes_z_axis():
    """Design Test 1: x=0 ∩ y=0 is the z-axis."""
    p = Plane(id=1, a=1, b=0, c=0, d=0)
    q = Plane(id=2, a=0, b=1, c=0, d=0)
    line = intersect_planes(p, q)
    assert line is not PARALLEL
    assert directions_parallel(line.direction, (0, 0, 1))
    assert line.contains(Point3D(0, 0, 0))
    assert line.contains(Point3D(0, 0, 5))
    assert line.contains(Point3D(0, 0, -2))
    assert line.plane_a == 1
    assert line.plane_b == 2
    assert line_lies_on_both_planes(line, p, q)
    for t in (-2, 0, 1, Fraction(5, 2)):
        assert p.contains(line.point_at(t))
        assert q.contains(line.point_at(t))


def test_intersect_parallel_planes():
    """Design Test 3: x=0 and x=2 do not meet."""
    p = Plane(id=1, a=1, b=0, c=0, d=0)
    q = Plane(id=2, a=1, b=0, c=0, d=-2)
    assert intersect_planes(p, q) is PARALLEL


def test_intersect_opposite_normal_same_plane():
    p = Plane(id=1, a=1, b=0, c=0, d=0)
    q = Plane(id=2, a=-1, b=0, c=0, d=0)
    assert intersect_planes(p, q) is PARALLEL


def test_intersect_planes_generic_and_invariant():
    p = Plane(id=1, a=1, b=1, c=0, d=-1)  # x + y = 1
    q = Plane(id=2, a=0, b=0, c=1, d=-2)  # z = 2
    line = intersect_planes(p, q)
    assert line is not PARALLEL
    assert line_lies_on_both_planes(line, p, q)
    assert line.contains(Point3D(1, 0, 2))
    assert line.contains(Point3D(0, 1, 2))


# ---------------------------------------------------------------------------
# Step 1.3 — three-plane intersection
# ---------------------------------------------------------------------------


def test_intersect_three_planes_axis_aligned():
    """Design Test 2: x=0, y=0, z=1 → (0,0,1)."""
    p = Plane(id=1, a=1, b=0, c=0, d=0)
    q = Plane(id=2, a=0, b=1, c=0, d=0)
    r = Plane(id=3, a=0, b=0, c=1, d=-1)
    point = intersect_three_planes(p, q, r)
    assert point == Point3D(0, 0, 1)
    assert point_lies_on_planes(point, p, q, r)


def test_intersect_three_planes_known_point():
    p = Plane(id=1, a=1, b=0, c=0, d=-1)  # x = 1
    q = Plane(id=2, a=0, b=1, c=0, d=-2)  # y = 2
    r = Plane(id=3, a=0, b=0, c=1, d=-3)  # z = 3
    point = intersect_three_planes(p, q, r)
    assert point == Point3D(1, 2, 3)
    assert point_lies_on_planes(point, p, q, r)


def test_intersect_three_planes_generic():
    p = Plane(id=1, a=1, b=1, c=1, d=-3)  # x + y + z = 3
    q = Plane(id=2, a=1, b=-1, c=0, d=0)  # x = y
    r = Plane(id=3, a=0, b=0, c=1, d=-1)  # z = 1
    point = intersect_three_planes(p, q, r)
    assert point == Point3D(1, 1, 1)
    assert point_lies_on_planes(point, p, q, r)


def test_intersect_three_planes_parallel_pencil():
    """Three planes through the z-axis: no unique point."""
    p = Plane(id=1, a=1, b=0, c=0, d=0)  # x = 0
    q = Plane(id=2, a=0, b=1, c=0, d=0)  # y = 0
    r = Plane(id=3, a=1, b=1, c=0, d=0)  # x + y = 0
    assert intersect_three_planes(p, q, r) is None


def test_intersect_three_planes_two_parallel():
    p = Plane(id=1, a=1, b=0, c=0, d=0)
    q = Plane(id=2, a=1, b=0, c=0, d=-2)
    r = Plane(id=3, a=0, b=1, c=0, d=0)
    assert intersect_three_planes(p, q, r) is None


# ---------------------------------------------------------------------------
# Step 1.4 — horizontal slice
# ---------------------------------------------------------------------------


def test_slice_plane_at_z_design_test_4():
    """Design Test 4: x+y+z-5=0 at z=2 → x+y-3=0."""
    plane = Plane(id=7, a=1, b=1, c=1, d=-5)
    line = slice_plane_at_z(plane, 2)
    assert line is not None
    assert line.a == 1
    assert line.b == 1
    assert line.c == -3
    assert line.source_plane_id == 7
    assert slice_lifts_to_plane(line, 2, plane)
    assert line.contains(Point2D(0, 3))
    assert line.contains(Point2D(3, 0))
    assert line.contains(Point2D(1, 2))


def test_slice_points_lift_to_the_plane():
    plane = Plane(id=1, a=2, b=-1, c=1, d=4)
    z = Fraction(-3, 2)
    line = slice_plane_at_z(plane, z)
    assert line is not None
    assert slice_lifts_to_plane(line, z, plane)
    for point2d in line.sample_points():
        lifted = Point3D(point2d.x, point2d.y, z)
        assert plane.contains(lifted)


def test_horizontal_plane_has_no_slice_line():
    plane = Plane(id=1, a=0, b=0, c=1, d=-5)  # z = 5
    assert slice_plane_at_z(plane, 2) is None
    assert slice_plane_at_z(plane, 5) is None


def test_build_slice_lines_drops_horizontal():
    slanted = Plane(id=1, a=1, b=0, c=1, d=0)
    horizontal = Plane(id=2, a=0, b=0, c=1, d=-1)
    lines = build_slice_lines([slanted, horizontal], 0)
    assert len(lines) == 1
    assert lines[0].source_plane_id == 1


def test_degenerate_line2d_rejected():
    with pytest.raises(ValueError, match="degenerate"):
        Line2D(a=0, b=0, c=1)
