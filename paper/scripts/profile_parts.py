"""Cost breakdown: event generation vs. sweep, and the cost of one slice VD."""
import random, time, json, sys
from vd3d.geometry.sampling import random_general_position_planes
from vd3d.events.all import generate_all_events, choose_z_below_all_events
from vd3d.events.triples import generate_triple_events
from vd3d.events.alignment import generate_alignment_events
from vd3d.events.slice import compute_vd_at_z
from vd3d.sweep import vertical_decomposition_3d

out=[]
for n in [int(x) for x in sys.argv[1].split(",")]:
    for seed in range(1, int(sys.argv[2])+1):
        planes = random_general_position_planes(random.Random(seed), n)
        rec={"n":n,"seed":seed}
        t=time.perf_counter(); tri=generate_triple_events(planes); rec["t_triples"]=time.perf_counter()-t
        t=time.perf_counter(); ali=generate_alignment_events(planes); rec["t_align"]=time.perf_counter()-t
        t=time.perf_counter(); ev=generate_all_events(tuple(planes)); rec["t_events"]=time.perf_counter()-t
        rec["events"]=len(ev)
        z0=choose_z_below_all_events(ev)
        try:
            t=time.perf_counter(); vd=compute_vd_at_z(planes,z0); rec["t_slice"]=time.perf_counter()-t
            rec["slice_cells"]=len(vd.cells)
        except BaseException as e:
            rec["t_slice"]=None; rec["slice_cells"]=None
        try:
            t=time.perf_counter(); res=vertical_decomposition_3d(planes, incremental=False)
            rec["t_total_ref"]=time.perf_counter()-t; rec["cells"]=len(res.cells); rec["status"]="ok"
        except BaseException as e:
            rec["t_total_ref"]=None; rec["status"]="fail"
        out.append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open(sys.argv[3],"w"), indent=1)
