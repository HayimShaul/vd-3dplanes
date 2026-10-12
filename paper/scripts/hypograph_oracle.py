"""Independent partition oracle for the hypograph VD.

Below the envelope and off every cell boundary -> exactly one open cell.
Above the envelope -> no cell.  Points on a boundary are reported separately.
"""
from __future__ import annotations

import random
from fractions import Fraction as F

from vd3d.geometry.sampling import random_general_position_planes
from vd3d.io import load_planes
from vd3d.lower_envelope import hypograph_vd
from vd3d.lower_envelope.dual import dual_points


def envelope_z(dual, x, y):
    return min(A * x + B * y + C for A, B, C in dual)


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


def scan(planes, label):
    vd = hypograph_vd(planes)
    dual = dual_points(vd.envelope.planes)
    rng = random.Random(12345)
    # Denominator 997 keeps samples off the small-integer features.
    pts = [(F(rng.randint(-5000, 5000), 997),
            F(rng.randint(-5000, 5000), 997),
            F(rng.randint(-9000, 9000), 997)) for _ in range(3000)]
    below = above = bdry = gap = overlap = spurious = 0
    ex = []
    for x, y, z in pts:
        if any(on_boundary(t, dual, x, y, z) for t in vd.trapezoids):
            bdry += 1
            continue
        hits = [t.id for t in vd.trapezoids if in_open_cell(t, dual, x, y, z)]
        if z < envelope_z(dual, x, y):
            below += 1
            if len(hits) == 0:
                gap += 1
                if len(ex) < 3:
                    ex.append(("GAP", str(x), str(y), str(z)))
            elif len(hits) > 1:
                overlap += 1
                if len(ex) < 3:
                    ex.append(("OVERLAP", str(x), str(y), str(z), hits))
        else:
            above += 1
            if hits:
                spurious += 1
                if len(ex) < 3:
                    ex.append(("SPURIOUS", str(x), str(y), str(z), hits))
    status = "ok" if (gap or overlap or spurious) == 0 else "FAIL"
    print(f"{label}: cells={len(vd.cells)} below={below} above={above} "
          f"on-boundary={bdry} | gaps={gap} overlaps={overlap} spurious={spurious}  [{status}]")
    for e in ex:
        print("   ", e)
    return gap == 0 and overlap == 0 and spurious == 0


ok = scan(load_planes("examples/three_planes.txt"), "three_planes")
for n in (3, 4, 5, 6):
    for seed in (0, 7, 42, 1, 2):
        planes = random_general_position_planes(random.Random(seed), n)
        ok &= scan(planes, f"n={n} seed={seed:<3}")
print("ALL OK" if ok else "FAILURES PRESENT")
