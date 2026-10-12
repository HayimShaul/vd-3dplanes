"""Benchmark harness for vd3d: runs the sweep on random instances."""
from __future__ import annotations
import json, random, sys, time, traceback
from vd3d.geometry.sampling import random_general_position_planes
from vd3d.events.all import generate_all_events
from vd3d.sweep import vertical_decomposition_3d, group_events_by_z


def run(n, seed, incremental):
    rng = random.Random(seed)
    planes = random_general_position_planes(rng, n)
    events = generate_all_events(tuple(planes))
    groups = group_events_by_z(events)
    simultaneous = sum(1 for g in groups if len(g.events) > 1)
    rec = dict(n=n, seed=seed, incremental=incremental,
               events=len(events), groups=len(groups), simultaneous=simultaneous,
               triples=sum(1 for e in events if e.type.name == "TRIPLE_INTERSECTION"),
               alignments=sum(1 for e in events if e.type.name == "VERTICAL_ALIGNMENT"))
    t0 = time.perf_counter()
    try:
        res = vertical_decomposition_3d(planes, incremental=incremental)
        rec["seconds"] = time.perf_counter() - t0
        rec["cells"] = len(res.cells)
        rec["intervals"] = len(res.intervals)
        rec["initial_cells"] = len(res.initial_vd.cells)
        rec["status"] = "ok"
    except BaseException as exc:  # noqa: BLE001
        rec["seconds"] = time.perf_counter() - t0
        rec["status"] = "fail"
        rec["error"] = type(exc).__name__
        rec["message"] = str(exc).splitlines()[0][:120] if str(exc) else ""
    return rec


if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1].split(",")]
    seeds = range(1, int(sys.argv[2]) + 1)
    modes = [m == "1" for m in sys.argv[3].split(",")]
    out = []
    for n in ns:
        for seed in seeds:
            for inc in modes:
                rec = run(n, seed, inc)
                out.append(rec)
                print(json.dumps(rec), flush=True)
    with open(sys.argv[4], "w") as fh:
        json.dump(out, fh, indent=1)
