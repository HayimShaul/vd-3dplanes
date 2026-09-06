import random

from vd3d.events import compute_intersection_lines, generate_triple_events
from vd3d.geometry import PARALLEL, intersect_planes, intersect_three_planes, normals_parallel
from vd3d.viz.random_geom import (
    random_general_position_planes,
    random_intersecting_planes,
    random_parallel_family,
    random_parallel_planes,
    random_triple_planes,
)


def test_random_intersecting_planes_meet():
    rng = random.Random(0)
    p, q = random_intersecting_planes(rng)
    assert not normals_parallel(p, q)
    assert intersect_planes(p, q) is not PARALLEL


def test_random_parallel_planes_do_not_meet():
    rng = random.Random(0)
    p, q = random_parallel_planes(rng)
    assert normals_parallel(p, q)
    assert intersect_planes(p, q) is PARALLEL


def test_random_triple_reconstructs_the_point():
    rng = random.Random(0)
    point, p, q, r = random_triple_planes(rng)
    assert intersect_three_planes(p, q, r) == point


def test_random_general_position_planes_counts():
    rng = random.Random(5)
    planes = random_general_position_planes(rng, 4)
    assert len(planes) == 4
    assert len(compute_intersection_lines(planes)) == 6
    assert len(generate_triple_events(planes)) == 4


def test_random_parallel_family_has_no_lines():
    rng = random.Random(5)
    planes = random_parallel_family(rng, 4)
    assert len(planes) == 4
    assert compute_intersection_lines(planes) == ()
    assert generate_triple_events(planes) == ()
