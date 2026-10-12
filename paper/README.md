# Paper

`vd3d-paper.tex` — *Vertical Decomposition of an Arrangement of Planes in R^3:
An Exact Sweep-Based Implementation*.

Covers both constructions in the repository: the `z`-sweep for the full
decomposition (§4) and the dual-hull hypograph construction for the cells below
the lower envelope (§5).

## Build

```bash
pdflatex vd3d-paper.tex && pdflatex vd3d-paper.tex   # twice, for refs/TOC
```

No bibtex run is needed; the bibliography is a `thebibliography` block.
Figures 2 and 3 are drawn by `pgfplots` at compile time from coordinates in
`data/fig_*.tex`; the only image file in the paper is the screenshot in
Figure 1. `make_figures.py` regenerates those coordinate files from the JSON —
it does not render any images itself.

Series colours are slots 1–3 of the data-viz reference palette, checked with
its validator on both the adjacent and all-pairs lists (all checks pass; the
aqua slot's sub-3:1 contrast warning is relieved by the per-series direct
labels and by the tables carrying the same numbers).

## Layout

```text
vd3d-paper.tex            the paper
figures/                  figures included by the paper
data/*.tex                generated tables, \input by the paper
data/benchmarks.json      sweep: sizes, timing, failure taxonomy (Tables 2, 3, 9)
data/breakdown.json       sweep: cost breakdown (Table 4)
data/hypograph.json       envelope vs sweep on identical instances (Table 5)
data/env_random.json      envelope scaling to n=96 (Table 6)
data/env_voronoi.json     Voronoi-lifting family, Theta(n) attained (Table 7)
data/env_arith.json       exact vs Qhull across coefficient ranges (Table 8)
data/env_oracle.json      partition oracle at scale (inline table, §8.10)
data/fig_*.tex            generated pgfplots sources for Figures 2 and 3
scripts/                  the measurement harnesses
```

## Reproducing the measurements

From the repository root with the virtualenv active:

```bash
# Sweep — Tables 2, 3, 9: sizes, timing, failure taxonomy
python paper/scripts/bench.py 3,4,5,6,7 20 1,0 paper/data/benchmarks.json

# Sweep — Table 4: cost breakdown
python paper/scripts/profile_parts.py 3,4,5,6,7 8 paper/data/breakdown.json

# Envelope vs sweep on identical instances — Table 5 (§8.5)
python paper/scripts/bench_hypograph.py 3,4,5,6,7 20 paper/data/hypograph.json

# Envelope experiments — Tables 6, 7, 8 and §8.6-8.10
python paper/scripts/bench_envelope_scaling.py random  paper/data/env_random.json
python paper/scripts/bench_envelope_scaling.py voronoi paper/data/env_voronoi.json
python paper/scripts/bench_envelope_scaling.py arith   paper/data/env_arith.json
python paper/scripts/bench_envelope_scaling.py oracle  paper/data/env_oracle.json

# regenerate the sweep .tex tables (2, 3, 4, 9) from the JSON
python paper/scripts/make_tables.py

# regenerate the running-time figures (2 and 3) from the JSON
python paper/scripts/make_figures.py

# Sweep robustness scan (§9), as shipped and with the diagnosed fix
python paper/scripts/robustness_scan.py
python paper/scripts/robustness_scan.py --patched

# Partition oracle on small instances (§5.4)
python paper/scripts/hypograph_oracle.py
```

### Notes on the samplers

`bench.py` and `bench_hypograph.py` use the library's
`random_general_position_planes`, whose coefficients come from `{-4..4}`. That
sampler raises `RuntimeError` at about n=26 — it runs out of general-position
coefficient tuples — so `bench_envelope_scaling.py` carries its own samplers:

- `sample_planes(rng, n, bound)` — wider coefficients, only `c != 0` enforced
  (the envelope construction does not need general position).
- `voronoi_planes(rng, n)` — the lifting `(a,b) -> 2ax + 2by + z - (a²+b²)`,
  whose minimization diagram is the Voronoi diagram of the sites. Every plane
  reaches the envelope, so this is where the `Theta(n)` bound is attained.

The oracles sample rationals with denominator 997 on purpose: round coordinates
land exactly on diagram edges and wall lines, where belonging to no *open* cell
is the correct answer rather than a coverage gap.

Timings in the paper were taken on a 13th Gen Intel Core i5-13420H running
Ubuntu 26.04 and CPython 3.14.4, single-threaded. The envelope float path needs
`scipy` (declared in the project's `dev` extra).
