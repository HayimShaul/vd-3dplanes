"""Experiments for the lower-envelope / hypograph VD (paper Section 8.5-8.8).

Four experiments, selected by the first argument:

  random    scaling on random planes (wide coefficient range)
  voronoi   the Voronoi-lifting family, where Theta(n) is attained
  arith     float-vs-exact cost and agreement across coefficient magnitudes
  oracle    partition oracle (cells tile the region below the envelope)

Usage:  python paper/scripts/bench_envelope_scaling.py <experiment> <out.json>

The library's own random_general_position_planes draws coefficients from
{-4..4} and cannot produce more than ~25 general-position planes, so these
experiments use local samplers.  The envelope construction needs only c != 0,
not full general position.
"""
from __future__ import annotations

import json
import random
import sys
import time
from fractions import Fraction as F

from vd3d.geometry.plane import Plane
from vd3d.lower_envelope import (
    compare_envelopes_float_exact,
    compare_hypograph_float_exact,
    hypograph_vd,
    hypograph_vd_float,
    lower_envelope,
    lower_envelope_float,
)
from vd3d.lower_envelope.dual import dual_points, plane_to_dual


# ---------------------------------------------------------------- samplers

def sample_planes(rng, n, bound):
    """n planes with c != 0 and pairwise-distinct dual points."""
    out, seen, pid = [], set(), 1
    while len(out) < n:
        c = rng.choice((-1, 1)) * rng.randint(1, bound)
        p = Plane(id=pid, a=rng.randint(-bound, bound), b=rng.randint(-bound, bound),
                  c=c, d=rng.randint(-bound, bound))
        key = plane_to_dual(p)
        if key in seen:
            continue
        seen.add(key)
        out.append(p)
        pid += 1
    return out


def voronoi_planes(rng, n, bound=200):
    """Site (a,b) -> plane 2a x + 2b y + z - (a^2+b^2) = 0.

    Then f_i(x,y) = |p - s_i|^2 - (x^2+y^2), so the minimization diagram is the
    Voronoi diagram of the sites and every plane is on the envelope.
    """
    sites, seen = [], set()
    while len(sites) < n:
        s = (rng.randint(-bound, bound), rng.randint(-bound, bound))
        if s in seen:
            continue
        seen.add(s)
        sites.append(s)
    planes = [Plane(id=i + 1, a=2 * a, b=2 * b, c=1, d=-(a * a + b * b))
              for i, (a, b) in enumerate(sites)]
    return planes, sites


# ---------------------------------------------------------------- helpers

def in_open_cell(trap, dual, x, y, z):
    if trap.left_x is not None and x <= trap.left_x:
        return False
    if trap.right_x is not None and x >= trap.right_x:
        return False
    if trap.bottom is not None:
        m, b = trap.bottom
        if y <= m * x + b:
            return False
    if trap.top is not None:
        m, b = trap.top
        if y >= m * x + b:
            return False
    A, B, C = dual[trap.plane_index]
    return z < A * x + B * y + C


def on_boundary(trap, dual, x, y, z):
    if trap.left_x is not None and x == trap.left_x:
        return True
    if trap.right_x is not None and x == trap.right_x:
        return True
    for mb in (trap.bottom, trap.top):
        if mb is not None and y == mb[0] * x + mb[1]:
            return True
    A, B, C = dual[trap.plane_index]
    return z == A * x + B * y + C


def timed(fn):
    t = time.perf_counter()
    value = fn()
    return value, time.perf_counter() - t


# ---------------------------------------------------------------- experiments

def exp_random(out):
    for n in (4, 8, 16, 24, 32, 48, 64, 96):
        for seed in range(1, 11):
            planes = sample_planes(random.Random(seed), n, 50)
            env, t_env = timed(lambda: lower_envelope(planes))
            vd, t_vd = timed(lambda: hypograph_vd(env))
            fenv, t_fenv = timed(lambda: lower_envelope_float(planes))
            fvd, t_fvd = timed(lambda: hypograph_vd_float(fenv))
            rec = dict(exp="random", n=n, seed=seed, bound=50,
                       active=len(env.active), facets=len(env.facets),
                       cells=len(vd.cells), t_env=t_env, t_vd=t_vd,
                       t_float=t_fenv + t_fvd,
                       env_agree=bool(compare_envelopes_float_exact(fenv, env).ok),
                       vd_agree=bool(compare_hypograph_float_exact(fvd, vd).ok))
            out.append(rec)
            print(json.dumps(rec), flush=True)


def exp_voronoi(out):
    for n in (4, 8, 16, 24, 32, 48, 64):
        for seed in range(1, 6):
            planes, sites = voronoi_planes(random.Random(seed), n)
            env, t_env = timed(lambda: lower_envelope(planes))
            vd, t_vd = timed(lambda: hypograph_vd(env))
            dual = dual_points(env.planes)

            # Semantic oracle: the ceiling of the cell containing (x,y,z) must be
            # the plane of the site nearest to (x,y).
            rng = random.Random(1000 + seed)
            checked = wrong = 0
            for _ in range(300):
                x = F(rng.randint(-200000, 200000), 997)
                y = F(rng.randint(-200000, 200000), 997)
                # strictly below the envelope, so the point is in the hypograph
                z = min(A * x + B * y + C for A, B, C in dual) - 1
                hits = [t for t in vd.trapezoids if in_open_cell(t, dual, x, y, z)]
                if len(hits) != 1:
                    continue
                best = min(range(len(sites)),
                           key=lambda i: (x - sites[i][0]) ** 2 + (y - sites[i][1]) ** 2)
                checked += 1
                if hits[0].plane_index != best:
                    wrong += 1
            rec = dict(exp="voronoi", n=n, seed=seed,
                       active=len(env.active), facets=len(env.facets),
                       cells=len(vd.cells), t_env=t_env, t_vd=t_vd,
                       all_active=len(env.active) == n,
                       nearest_checked=checked, nearest_wrong=wrong)
            out.append(rec)
            print(json.dumps(rec), flush=True)


def exp_arith(out):
    for bound in (4, 50, 1000, 10 ** 6, 10 ** 12):
        for n in (8, 16, 32):
            for seed in range(1, 11):
                planes = sample_planes(random.Random(seed), n, bound)
                env, t_env = timed(lambda: lower_envelope(planes))
                vd, t_vd = timed(lambda: hypograph_vd(env))
                fenv, t_fenv = timed(lambda: lower_envelope_float(planes))
                fvd, t_fvd = timed(lambda: hypograph_vd_float(fenv))
                ec = compare_envelopes_float_exact(fenv, env)
                hc = compare_hypograph_float_exact(fvd, vd)
                rec = dict(exp="arith", n=n, seed=seed, bound=bound,
                           active=len(env.active), cells=len(vd.cells),
                           t_exact=t_env + t_vd, t_float=t_fenv + t_fvd,
                           env_agree=bool(ec.ok), vd_agree=bool(hc.ok),
                           env_active_match=bool(ec.active_match),
                           env_facets_match=bool(ec.facets_match),
                           max_vertex_error=ec.max_vertex_error)
                out.append(rec)
                print(json.dumps(rec), flush=True)


def exp_oracle(out):
    for n in (4, 8, 16, 32):
        for seed in range(1, 6):
            planes = sample_planes(random.Random(seed), n, 50)
            vd = hypograph_vd(planes)
            dual = dual_points(vd.envelope.planes)
            rng = random.Random(7000 + seed)
            below = above = bdry = gap = overlap = spurious = 0
            for _ in range(1500):
                x = F(rng.randint(-50000, 50000), 997)
                y = F(rng.randint(-50000, 50000), 997)
                z = F(rng.randint(-200000, 200000), 997)
                if any(on_boundary(t, dual, x, y, z) for t in vd.trapezoids):
                    bdry += 1
                    continue
                hits = sum(1 for t in vd.trapezoids if in_open_cell(t, dual, x, y, z))
                if z < min(A * x + B * y + C for A, B, C in dual):
                    below += 1
                    if hits == 0:
                        gap += 1
                    elif hits > 1:
                        overlap += 1
                else:
                    above += 1
                    if hits:
                        spurious += 1
            rec = dict(exp="oracle", n=n, seed=seed, cells=len(vd.cells),
                       below=below, above=above, boundary=bdry,
                       gaps=gap, overlaps=overlap, spurious=spurious)
            out.append(rec)
            print(json.dumps(rec), flush=True)


EXPERIMENTS = {"random": exp_random, "voronoi": exp_voronoi,
               "arith": exp_arith, "oracle": exp_oracle}

if __name__ == "__main__":
    which, path = sys.argv[1], sys.argv[2]
    records: list[dict] = []
    EXPERIMENTS[which](records)
    json.dump(records, open(path, "w"), indent=1)
