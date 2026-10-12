"""Regenerate paper/data/*.tex from paper/data/benchmarks.json and breakdown.json.

Usage (from the repository root, with the venv active):

    python paper/scripts/bench.py 3,4,5,6,7 20 1,0 paper/data/benchmarks.json
    python paper/scripts/bench.py 8 12 1,0 /tmp/b8.json     # then merge by hand
    python paper/scripts/profile_parts.py 3,4,5,6,7 8 paper/data/breakdown.json
    python paper/scripts/make_tables.py
"""
from __future__ import annotations

import json
import math
import statistics as st
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"


def cls(rec):
    msg = rec.get("message", "")
    if "arrangement faces" in msg:
        return "a"
    if "combinatorially" in msg:
        return "b"
    return "ok" if rec["status"] == "ok" else "other"


def power_fit(rows, ns, key, incremental):
    xs, ys = [], []
    for n in ns:
        ok = [r for r in rows if r["n"] == n and r["incremental"] is incremental and r["status"] == "ok"]
        if not ok:
            continue
        xs.append(math.log(n))
        ys.append(math.log(st.median(r[key] for r in ok)))
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


def main() -> None:
    rows = json.loads((DATA / "benchmarks.json").read_text())
    parts = json.loads((DATA / "breakdown.json").read_text())
    ns = sorted({r["n"] for r in rows})

    def sel(n, inc):
        return [r for r in rows if r["n"] == n and r["incremental"] is inc]

    out = [r"\begin{table}[htbp]", r"\centering", r"\begin{tabular}{rrrrrrrr}", r"\toprule",
           r"$n$ & inst. & triples & \multicolumn{2}{c}{alignments} & \multicolumn{2}{c}{3D cells} & slice cells \\",
           r"& & $\binom{n}{3}$ & med & max & med & max & med \\", r"\midrule"]
    for n in ns:
        rs, ok = sel(n, True), [r for r in sel(n, True) if r["status"] == "ok"]
        al = [r["alignments"] for r in rs]
        out.append(f"{n} & {len(rs)} & {math.comb(n, 3)} & {st.median(al):.0f} & {max(al)} & "
                   f"{st.median(r['cells'] for r in ok):.0f} & {max(r['cells'] for r in ok)} & "
                   f"{st.median(r['initial_cells'] for r in ok):.0f} \\\\")
    out += [r"\bottomrule", r"\end{tabular}",
            r"\caption{Instance sizes on random general-position planes (seeds $1..20$; $1..12$ for $n=8$). "
            r"``slice cells'' is the number of trapezoids in the initial slice decomposition. Cell statistics "
            r"are over instances that completed.}", r"\label{tab:combinatorics}", r"\end{table}"]
    (DATA / "table_combinatorics.tex").write_text("\n".join(out) + "\n")

    fit_i, fit_r = power_fit(rows, ns, "seconds", True), power_fit(rows, ns, "seconds", False)
    out = [r"\begin{table}[htbp]", r"\centering", r"\begin{tabular}{rrrrrrrr}", r"\toprule",
           r"& \multicolumn{3}{c}{reference path (s)} & \multicolumn{3}{c}{incremental path (s)} & \\",
           r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
           r"$n$ & median & max & per event & median & max & per event & ratio \\", r"\midrule"]
    for n in ns:
        okr = [r for r in sel(n, False) if r["status"] == "ok"]
        oki = [r for r in sel(n, True) if r["status"] == "ok"]
        if not okr or not oki:
            continue
        mr, mi = st.median(r["seconds"] for r in okr), st.median(r["seconds"] for r in oki)
        out.append(f"{n} & {mr:.3f} & {max(r['seconds'] for r in okr):.3f} & "
                   f"{st.median(r['seconds'] / max(r['events'], 1) for r in okr):.4f} & "
                   f"{mi:.3f} & {max(r['seconds'] for r in oki):.3f} & "
                   f"{st.median(r['seconds'] / max(r['events'], 1) for r in oki):.4f} & {mi / mr:.2f} \\\\")
    out += [r"\midrule",
            rf"\multicolumn{{8}}{{l}}{{Power-law fit of median time: reference $t \propto n^{{{fit_r:.1f}}}$, "
            rf"incremental $t \propto n^{{{fit_i:.1f}}}$.}} \\",
            r"\bottomrule", r"\end{tabular}",
            r"\caption{Wall-clock time, single-threaded. ``per event'' is median total time divided by the "
            r"number of events. ``ratio'' is incremental median over reference median: values above $1$ mean "
            r"the incremental path is slower, because it also recomputes the slice to assert agreement.}",
            r"\label{tab:timing}", r"\end{table}"]
    (DATA / "table_timing.tex").write_text("\n".join(out) + "\n")

    out = [r"\begin{table}[htbp]", r"\centering", r"\begin{tabular}{rrrrrrr}", r"\toprule",
           r"& & \multicolumn{2}{c}{incremental path} & \multicolumn{2}{c}{reference path} & \\",
           r"\cmidrule(lr){3-4}\cmidrule(lr){5-6}",
           r"$n$ & inst. & class (a) & class (b) & class (a) & class (b) & completion rate \\", r"\midrule"]
    tot = dict(inst=0, ia=0, ib=0, ra=0, rb=0, iok=0, rok=0)
    for n in ns:
        ri, rr = sel(n, True), sel(n, False)
        vals = dict(ia=sum(1 for r in ri if cls(r) == "a"), ib=sum(1 for r in ri if cls(r) == "b"),
                    ra=sum(1 for r in rr if cls(r) == "a"), rb=sum(1 for r in rr if cls(r) == "b"),
                    iok=sum(1 for r in ri if r["status"] == "ok"),
                    rok=sum(1 for r in rr if r["status"] == "ok"), inst=len(ri))
        out.append(f"{n} & {len(ri)} & {vals['ia']} & {vals['ib']} & {vals['ra']} & {vals['rb']} & "
                   f"{100 * vals['iok'] / len(ri):.0f}\\% / {100 * vals['rok'] / len(rr):.0f}\\% \\\\")
        for k in tot:
            tot[k] += vals[k]
    out += [r"\midrule",
            f"all & {tot['inst']} & {tot['ia']} & {tot['ib']} & {tot['ra']} & {tot['rb']} & "
            f"{100 * tot['iok'] / tot['inst']:.0f}\\% / {100 * tot['rok'] / tot['inst']:.0f}\\% \\\\",
            r"\bottomrule", r"\end{tabular}",
            r"\caption{Aborted runs by class. (a) face-extraction failure in the planar layer; "
            r"(b) incremental-vs-reference divergence caught by the differential assertion. "
            r"Completion rate is given as incremental / reference.}",
            r"\label{tab:failures}", r"\end{table}"]
    (DATA / "table_failures.tex").write_text("\n".join(out) + "\n")

    out = [r"\begin{table}[htbp]", r"\centering", r"\begin{tabular}{rrrrrr}", r"\toprule",
           r"$n$ & triples (s) & alignments (s) & one slice VD (s) & total sweep (s) & event gen.\ share \\",
           r"\midrule"]
    for n in sorted({p["n"] for p in parts}):
        ps = [p for p in parts if p["n"] == n]
        ok = [p for p in ps if p["status"] == "ok"]
        if not ok:
            continue
        out.append(f"{n} & {st.median(p['t_triples'] for p in ps):.5f} & "
                   f"{st.median(p['t_align'] for p in ps):.4f} & "
                   f"{st.median(p['t_slice'] for p in ps if p['t_slice'] is not None):.4f} & "
                   f"{st.median(p['t_total_ref'] for p in ok):.3f} & "
                   f"{100 * st.median(p['t_events'] / p['t_total_ref'] for p in ok):.1f}\\% \\\\")
    out += [r"\bottomrule", r"\end{tabular}",
            r"\caption{Cost breakdown (medians, 8 seeds per $n$, reference path). ``one slice VD'' is a single "
            r"from-scratch \texttt{compute\_vd\_at\_z}. The last column is event generation as a fraction of "
            r"total sweep time.}", r"\label{tab:breakdown}", r"\end{table}"]
    (DATA / "table_breakdown.tex").write_text("\n".join(out) + "\n")

    print(f"wrote 4 tables; time fits: reference n^{fit_r:.2f}, incremental n^{fit_i:.2f}")


if __name__ == "__main__":
    main()
