"""Benchmark + robustness for the lower-envelope / hypograph VD path."""
from __future__ import annotations

import json, random, sys, time

from vd3d.geometry.sampling import random_general_position_planes
from vd3d.lower_envelope import (
    compare_hypograph_float_exact, hypograph_vd, hypograph_vd_float,
    lower_envelope, lower_envelope_float,
)
from vd3d.sweep import vertical_decomposition_3d

out = []
for n in [int(x) for x in sys.argv[1].split(",")]:
    for seed in range(1, int(sys.argv[2]) + 1):
        planes = random_general_position_planes(random.Random(seed), n)
        rec = {"n": n, "seed": seed}

        t = time.perf_counter()
        try:
            env = lower_envelope(planes); rec["t_env_exact"] = time.perf_counter() - t
            rec["active"] = len(env.active); rec["facets"] = len(env.facets)
            rec["env_status"] = "ok"
        except BaseException as e:
            rec["t_env_exact"] = time.perf_counter() - t
            rec["env_status"] = "fail"; rec["env_err"] = type(e).__name__
            out.append(rec); print(json.dumps(rec), flush=True); continue

        t = time.perf_counter()
        try:
            vd = hypograph_vd(env); rec["t_vd_exact"] = time.perf_counter() - t
            rec["hypo_cells"] = len(vd.cells); rec["vd_status"] = "ok"
        except BaseException as e:
            rec["t_vd_exact"] = time.perf_counter() - t
            rec["vd_status"] = "fail"; rec["vd_err"] = type(e).__name__
            out.append(rec); print(json.dumps(rec), flush=True); continue

        t = time.perf_counter()
        fenv = lower_envelope_float(planes); fvd = hypograph_vd_float(fenv)
        rec["t_float"] = time.perf_counter() - t
        rec["float_exact_agree"] = bool(compare_hypograph_float_exact(fvd, vd).ok)

        # does the full sweep survive this same instance?
        t = time.perf_counter()
        try:
            res = vertical_decomposition_3d(planes)
            rec["t_sweep"] = time.perf_counter() - t
            rec["sweep_status"] = "ok"; rec["sweep_cells"] = len(res.cells)
            rec["sweep_floor_none"] = sum(1 for c in res.cells if c.floor is None)
        except BaseException as e:
            rec["t_sweep"] = time.perf_counter() - t
            rec["sweep_status"] = "fail"
            rec["sweep_err"] = ("a" if "arrangement faces" in str(e)
                                else "b" if "combinatorially" in str(e) else type(e).__name__)
        out.append(rec); print(json.dumps(rec), flush=True)

json.dump(out, open(sys.argv[3], "w"), indent=1)
