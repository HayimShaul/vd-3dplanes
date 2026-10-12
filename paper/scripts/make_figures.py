"""Emit pgfplots figure sources for the paper from the measured JSON.

    python paper/scripts/make_figures.py

Writes paper/data/fig_time_constructions.tex and paper/data/fig_time_envelope_stages.tex,
which the paper \\input{}s.  The plots are drawn by LaTeX (pgfplots); no image
files are produced.  Reads only paper/data/*.json - it runs no benchmarks.
"""
from __future__ import annotations

import json
import math
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"


def medians(rows, vkey, where=None):
    sel = [r for r in rows if (where is None or where(r))]
    ns = sorted({r["n"] for r in sel})
    return ns, [st.median(vkey(r) for r in sel if r["n"] == n) for n in ns]


def slope(ns, ys):
    xs = [math.log(n) for n in ns]
    ls = [math.log(y) for y in ys]
    mx, my = sum(xs) / len(xs), sum(ls) / len(ls)
    return sum((x - mx) * (l - my) for x, l in zip(xs, ls)) / sum((x - mx) ** 2 for x in xs)


def coords(ns, ys):
    return " ".join(f"({n},{y:.6g})" for n, y in zip(ns, ys))


def crossover(ns, a, b):
    """n at which series a overtakes b, interpolated in log-log."""
    prev = None
    for n, ra in zip(ns, (x / y for x, y in zip(a, b))):
        if prev and prev[1] < 1 <= ra:
            (n0, r0) = prev
            x0, x1 = math.log(n0), math.log(n)
            y0, y1 = math.log(r0), math.log(ra)
            return math.exp(x0 + (0 - y0) * (x1 - x0) / (y1 - y0))
        prev = (n, ra)
    return None


AXIS = r"""    width=\linewidth, height=%(h)s,
    xmode=log, ymode=log, log basis x=10,
    xlabel={number of planes $n$}, ylabel={%(ylab)s},
    xtick={%(xtick)s}, xticklabels={%(xtick)s},
    ytick={1e-4,1e-3,1e-2,1e-1,1,1e1,1e2},
    log ticks with fixed point,
    xmin=%(xmin)s, xmax=%(xmax)s,%(ymm)s
    grid=both,
    major grid style={line width=0.3pt, draw=gridgray},
    minor grid style={line width=0.2pt, draw=gridgray, opacity=0.55},
    axis lines=left, axis line style={draw=inktwo, line width=0.5pt},
    tick label style={font=\footnotesize, color=inktwo},
    label style={font=\small, color=ink},
    legend cell align=left,
    legend style={draw=none, fill=none, font=\footnotesize,
                  at={%(legpos)s}, anchor=%(leganc)s, row sep=1pt,
                  /tikz/every even column/.append style={column sep=8pt}},
    legend columns=%(legcols)s,
"""


def fig_constructions() -> None:
    sweep = json.loads((DATA / "benchmarks.json").read_text())
    hypo = json.loads((DATA / "hypograph.json").read_text())
    envr = json.loads((DATA / "env_random.json").read_text())

    ok = lambda r: r["status"] == "ok"
    n_ref, t_ref = medians(sweep, lambda r: r["seconds"], lambda r: ok(r) and not r["incremental"])
    n_inc, t_inc = medians(sweep, lambda r: r["seconds"], lambda r: ok(r) and r["incremental"])
    n_env, t_env = medians(hypo, lambda r: r["t_env_exact"] + r["t_vd_exact"],
                           lambda r: r["vd_status"] == "ok")
    n_wid, t_wid = medians(envr, lambda r: r["t_env"] + r["t_vd"])

    axis = AXIS % dict(h=r"6.2cm", ylab="median running time (s)",
                       xtick="3,4,6,8,12,16,24,32,48,64,96",
                       xmin="2.7", xmax="150",
                       ymm="\n    ymin=4e-4, ymax=4e2,",
                       legpos="{(0.5,-0.20)}", leganc="north", legcols="2")
    body = rf"""\begin{{figure}}[htbp]
\centering
\begin{{tikzpicture}}
\begin{{axis}}[
{axis}]
\addplot[serA, line width=0.9pt, mark=*, mark size=1.9pt]
  coordinates {{{coords(n_ref, t_ref)}}};
\addlegendentry{{VD of the arrangement, reference sweep}}
\addplot[serB, line width=0.9pt, mark=square*, mark size=1.7pt]
  coordinates {{{coords(n_inc, t_inc)}}};
\addlegendentry{{VD of the arrangement, incremental sweep}}
\addplot[serC, line width=0.9pt, mark=triangle*, mark size=2.2pt]
  coordinates {{{coords(n_env, t_env)}}};
\addlegendentry{{VD of the lower envelope (same instances)}}
\addplot[serC, line width=0.9pt, dashed, mark=triangle, mark size=2.2pt,
         mark options={{solid, fill=white}}]
  coordinates {{{coords(n_wid, t_wid)}}};
\addlegendentry{{VD of the lower envelope (wider sampler)}}
\node[serB, font=\scriptsize, anchor=west] at (axis cs:{n_inc[-1]*1.1},{t_inc[-1]*2.2:.4g})
  {{$\propto n^{{{slope(n_inc, t_inc):.1f}}}$}};
\node[serA, font=\scriptsize, anchor=west] at (axis cs:{n_ref[-1]*1.1},{t_ref[-1]*0.42:.4g})
  {{$\propto n^{{{slope(n_ref, t_ref):.1f}}}$}};
\node[serC, font=\scriptsize, anchor=west] at (axis cs:{n_wid[-1]*1.07},{t_wid[-1]*1.5:.4g})
  {{$\propto n^{{{slope(n_wid, t_wid):.1f}}}$}};
\end{{axis}}
\end{{tikzpicture}}
\caption{{Median running time of the two constructions, log--log.  Both sweep
paths and the solid envelope series are the general-position instances of
Tables~\ref{{tab:timing}} and \ref{{tab:hypograph}} (coefficients in
$\{{-4,\dots,4\}}$, $n \le 8$); the dashed envelope series is the wider sampler of
Table~\ref{{tab:env-random}}, which reaches $n = 96$.  The two distributions agree
where they overlap.  Exponents are least-squares fits over each series' own
range; the sweep's is measured over $3 \le n \le 8$ and is not an asymptotic
claim (Section~\ref{{sec:bench-env-scaling}}).}}
\label{{fig:time-constructions}}
\end{{figure}}
"""
    (DATA / "fig_time_constructions.tex").write_text(body)
    print("wrote data/fig_time_constructions.tex")


def fig_stages() -> None:
    envr = json.loads((DATA / "env_random.json").read_text())
    envv = json.loads((DATA / "env_voronoi.json").read_text())

    n, hull = medians(envr, lambda r: r["t_env"])
    _, trap = medians(envr, lambda r: r["t_vd"])
    _, flt = medians(envr, lambda r: r["t_float"])
    n2, hull2 = medians(envv, lambda r: r["t_env"])
    _, trap2 = medians(envv, lambda r: r["t_vd"])

    cross = crossover(n, hull, trap)
    ratio = hull[-1] / trap[-1]

    def panel(title, xtick, xmin, xmax, ylab, series, legpos="{(0.03,0.97)}",
              leganc="north west"):
        ax = AXIS % dict(h=r"5.1cm", ylab=ylab, xtick=xtick, xmin=xmin, xmax=xmax,
                         ymm="\n    ymin=1e-4, ymax=2e1,", legpos=legpos, leganc=leganc,
                         legcols="1")
        plots = "\n".join(series)
        return (f"\\begin{{tikzpicture}}\n\\begin{{axis}}[\n{ax}"
                f"    title={{{title}}}, title style={{font=\\small, color=inktwo}},\n]\n"
                f"{plots}\n\\end{{axis}}\n\\end{{tikzpicture}}")

    left = panel(r"(a) random planes", "4,8,16,32,64,96", "3.6", "115",
                 "median running time (s)", [
        rf"\addplot[serA, line width=0.9pt, mark=*, mark size=1.9pt] coordinates {{{coords(n, hull)}}};",
        r"\addlegendentry{exact lower hull}",
        rf"\addplot[serB, line width=0.9pt, mark=square*, mark size=1.7pt] coordinates {{{coords(n, trap)}}};",
        r"\addlegendentry{exact trapezoids}",
        rf"\addplot[serC, line width=0.9pt, mark=triangle*, mark size=2.2pt] coordinates {{{coords(n, flt)}}};",
        r"\addlegendentry{float (both stages)}",
    ])
    right = panel(r"(b) Voronoi lifting", "4,8,16,32,64", "3.6", "78", "", [
        rf"\addplot[serA, line width=0.9pt, mark=*, mark size=1.9pt] coordinates {{{coords(n2, hull2)}}};",
        r"\addlegendentry{exact lower hull}",
        rf"\addplot[serB, line width=0.9pt, mark=square*, mark size=1.7pt] coordinates {{{coords(n2, trap2)}}};",
        r"\addlegendentry{exact trapezoids}",
    ])

    body = rf"""\begin{{figure}}[htbp]
\centering
\begin{{minipage}}{{0.49\textwidth}}\centering
{left}
\end{{minipage}}\hfill
\begin{{minipage}}{{0.49\textwidth}}\centering
{right}
\end{{minipage}}
\caption{{The envelope pipeline's two exact stages, log--log.  \textbf{{(a)}} On
random planes the $O(n^4)$ hull enumeration overtakes the trapezoid
construction at $n \approx {cross:.0f}$ and is ${ratio:.1f}\times$ its cost by
$n = 96$, because the envelope --- and therefore the trapezoid work --- stops
growing while the triple enumeration does not.  \textbf{{(b)}} On the Voronoi
lifting the ordering is reversed at every size: all $n$ planes stay active, so
the trapezoid stage has $\Theta(n)$ envelope features to decompose and dominates
throughout.  Panel (a) also shows the price of exactness, as the gap between the
two exact stages and the float series.}}
\label{{fig:time-stages}}
\end{{figure}}
"""
    (DATA / "fig_time_envelope_stages.tex").write_text(body)
    print(f"wrote data/fig_time_envelope_stages.tex  (crossover n={cross:.1f}, ratio {ratio:.2f})")


if __name__ == "__main__":
    fig_constructions()
    fig_stages()
