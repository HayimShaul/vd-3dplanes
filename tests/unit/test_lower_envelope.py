"""Exact lower envelope and hypograph VD, checked against the float reference."""

from __future__ import annotations

import random

from vd3d.geometry.sampling import random_general_position_planes
from vd3d.io import load_planes
from vd3d.lower_envelope import (
    compare_envelopes_float_exact,
    compare_hypograph_float_exact,
    hypograph_vd,
    hypograph_vd_float,
    lower_envelope,
    lower_envelope_float,
)


def test_three_planes_envelope_matches_float():
    planes = load_planes("examples/three_planes.txt")
    exact = lower_envelope(planes)
    flt = lower_envelope_float(planes)
    cmp = compare_envelopes_float_exact(flt, exact)
    assert cmp.ok
    assert len(exact.facets) == 1
    assert exact.vertices[0].x == 1
    assert exact.vertices[0].y == 2
    assert exact.vertices[0].z == 3


def test_three_planes_hypograph_vd_matches_float():
    planes = load_planes("examples/three_planes.txt")
    exact = hypograph_vd(planes)
    flt = hypograph_vd_float(planes)
    cmp = compare_hypograph_float_exact(flt, exact)
    assert cmp.ok
    assert len(exact.cells) == 4
    assert all(cell.floor is None for cell in exact.cells)
    assert all(cell.ceiling is not None for cell in exact.cells)


def test_random_envelope_and_vd_match_float():
    for seed in (42, 7, 0):
        planes = random_general_position_planes(random.Random(seed), 4)
        exact_env = lower_envelope(planes)
        float_env = lower_envelope_float(planes)
        assert compare_envelopes_float_exact(float_env, exact_env).ok

        exact_vd = hypograph_vd(exact_env)
        float_vd = hypograph_vd_float(float_env)
        assert compare_hypograph_float_exact(float_vd, exact_vd).ok
        assert all(cell.floor is None for cell in exact_vd.cells)


def test_inactive_planes_omitted_from_active():
    planes = random_general_position_planes(random.Random(7), 6)
    env = lower_envelope(planes)
    assert len(env.active) <= len(planes)
    assert set(env.active) == {i for f in env.facets for i in f}
