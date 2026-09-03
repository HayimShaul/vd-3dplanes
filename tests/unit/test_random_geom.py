import random

from vd3d.geometry import PARALLEL, intersect_planes, intersect_three_planes, normals_parallel
from vd3d.viz.random_geom import (
    random_intersecting_planes,
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
