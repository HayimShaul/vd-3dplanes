"""Robustness scan reported in Section 8.

Runs the sweep on random general-position instances and classifies aborts.
Pass ``--patched`` to monkey-patch the diagnosed fix of Section 8.2 (normalise
the direction vector before taking the left normal in ``_left_offset``) without
modifying the library.

    python paper/scripts/robustness_scan.py
    python paper/scripts/robustness_scan.py --patched
"""
from __future__ import annotations

import argparse
import random


def apply_patch() -> None:
    import vd3d.arrangement2d.build as build
    from vd3d.arrangement2d.build import _point_on_supporting_line
    from vd3d.arrangement2d.predicates import left_normal

    def normalised_left_offset(half_edge, lines, vertices):
        dx, dy = half_edge.direction
        scale = max(abs(dx), abs(dy))
        lx, ly = left_normal((dx / scale, dy / scale))
        point = _point_on_supporting_line(half_edge, lines, vertices)
        return lx * point.x + ly * point.y

    build._left_offset = normalised_left_offset


def classify(exc: BaseException) -> str:
    text = str(exc)
    if "arrangement faces" in text:
        return "a"
    if "combinatorially" in text:
        return "b"
    return type(exc).__name__


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--patched", action="store_true")
    parser.add_argument("--ns", default="3,4,5,6")
    parser.add_argument("--seeds", type=int, default=40)
    args = parser.parse_args()

    if args.patched:
        apply_patch()

    from vd3d.geometry.sampling import random_general_position_planes
    from vd3d.sweep import vertical_decomposition_3d

    for incremental in (False, True):
        completed, kinds = 0, {}
        for n in (int(x) for x in args.ns.split(",")):
            for seed in range(1, args.seeds + 1):
                planes = random_general_position_planes(random.Random(seed), n)
                try:
                    vertical_decomposition_3d(planes, incremental=incremental)
                    completed += 1
                except BaseException as exc:  # noqa: BLE001
                    kind = classify(exc)
                    kinds[kind] = kinds.get(kind, 0) + 1
        label = "patched" if args.patched else "as shipped"
        print(f"{label}, incremental={incremental}: completed {completed}, aborts {kinds}", flush=True)


if __name__ == "__main__":
    main()
