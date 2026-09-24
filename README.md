# Vertical decomposition of planes in 3D

Sweep-based vertical decomposition of an arrangement of planes, following
[design.md](design.md).

This is an incremental implementation. Each phase is small on purpose. Do not
start the next phase until the current one has passed its human gate.

## Status

**Phase 10 — robustness.** Simultaneous event groups, query lines through
vertices or overlapping an arrangement edge, and a few non-general-position
fixtures (four planes at a point, vertical input planes). Exact
`Fraction` arithmetic.

Phases 0–9 are in place. Do not start further work until you have checked
the Phase 10 figures against the captions (through-vertex walk still
enters the triangle; overlap paints every face incident to L;
simultaneous groups are one before/after pair).

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Requires Python 3.11+. The virtual environment is gitignored.

## Compute a vertical decomposition

```bash
# from a plane file (each line: a b c d, or id a b c d; # comments ok)
python -m vd3d examples/three_planes.txt

# debug: N random general-position planes, optional seed
python -m vd3d --random-planes 5 --seed 42

# write the sampled planes for reuse
python -m vd3d --random-planes 5 --seed 42 --write-planes /tmp/planes.txt

# interactive 3D view (n/p cells, +/- zoom, arrows rotate, q quit)
python -m vd3d examples/three_planes.txt --gui
python -m vd3d --random-planes 4 --seed 42 --gui

# sweep-plane before/after each event (n/p events, pan/zoom, z/Z α)
python -m vd3d examples/three_planes.txt --show-sweep
python -m vd3d --random-planes 4 --seed 42 --show-sweep
```

Coefficients are exact (`int` or Fraction strings like `3/2`). The plane
equation is `a*x + b*y + c*z + d = 0`. Output lists every 3D cell with its
floor, ceiling, vertical walls, and z-extent. With `--gui`, input planes are
transparent grey, intersection lines dark grey, and the selected cell's walls
semi-transparent red (unbounded faces clipped to a viewing cube). With
`--show-sweep`, each event shows the 2D VD at `z±α` side by side: solid
arrangement lines, dashed vertical walls, adjacency-coloured cells with ids,
and a red event marker.

## How to review a step

Every step ships three things: unit tests, inspectable figures, and docs.
A step is done only when all three hold:

1. That step's unit tests are green.
2. You open the review surface and the figures match the captions.
3. The docs for that step list the invariants that are actually checked.

Do not start the next step if a figure and its caption disagree.

Commands for the current phase:

```bash
# Phase 10
pytest tests/unit/test_robustness.py tests/unit/test_zone.py tests/visual/test_robustness_visual.py
python -m vd3d.viz.gallery --step 10
python -m vd3d.viz.viewer --phase 10
python -m vd3d.viz.viewer --phase 10 --seed 42 --n 4
```

Confirm: `y=x` through the origin marks the vertex and still paints the
triangle; overlapping `y=0` paints every face that touches that side;
a simultaneous group is **one** before/after pair. Unmatched cells are
**black** only near the event strip.

`--n` sets the number of random planes (phases 1, 5–10) or lines (phases 2–4).
Omit it and each random scene picks a small count. `g` resamples with a new
seed and keeps `n`.

Earlier phases:

```bash
# Phase 9
pytest tests/unit/test_incremental.py tests/visual/test_incremental_visual.py
python -m vd3d.viz.gallery --step 9
python -m vd3d.viz.viewer --phase 9

# Phase 8
pytest tests/unit/test_sweep.py tests/visual/test_sweep_visual.py
python -m vd3d.viz.gallery --step 8
python -m vd3d.viz.viewer --phase 8

# Phase 7
pytest tests/unit/test_event_list.py tests/visual/test_event_list_visual.py
python -m vd3d.viz.gallery --step 7
python -m vd3d.viz.viewer --phase 7

# Phase 6
pytest tests/unit/test_alignment.py tests/visual/test_alignment_visual.py
python -m vd3d.viz.gallery --step 6
python -m vd3d.viz.viewer --phase 6

# Phase 5
pytest tests/unit/test_events.py tests/visual/test_events_visual.py
python -m vd3d.viz.gallery --step 5
python -m vd3d.viz.viewer --phase 5

# Phase 4
pytest tests/unit/test_zone.py tests/visual/test_zone_visual.py
python -m vd3d.viz.gallery --step 4
python -m vd3d.viz.viewer --phase 4

# Phase 3
pytest tests/unit/test_vertical_decomposition.py tests/visual/test_vd2d_visual.py
python -m vd3d.viz.gallery --step 3

# Phase 2
pytest tests/unit/test_arrangement2d.py tests/visual/test_arrangement2d_visual.py
python -m vd3d.viz.gallery --step 2

# Phase 1 (rotate the 3D scenes)
pytest tests/unit/test_geometry.py tests/unit/test_scenes.py
python -m vd3d.viz.viewer --phase 1
python -m vd3d.viz.viewer --phase 1 --seed 42

# Phase 0 (gallery pipeline)
pytest tests/unit/test_scalar.py tests/unit/test_linalg.py tests/unit/test_architecture.py tests/visual/test_gallery_smoke.py
python -m vd3d.viz.gallery --step 0
```

Later phases use the same pattern:

```bash
pytest tests/unit/test_<module>.py
python -m vd3d.viz.gallery --step N
```

## Conventions and invariants

- [docs/conventions.md](docs/conventions.md) — axes, "vertical", exact arithmetic, general position
- [docs/geometry.md](docs/geometry.md) — plane equation, orientation, intersections, slices
- [docs/arrangement2d.md](docs/arrangement2d.md) — DCEL, CCW order, Euler characteristic
- [docs/vertical_decomposition.md](docs/vertical_decomposition.md) — y-parallel walls, trapezoid cells
- [docs/zone.md](docs/zone.md) — query-line zone, supporting-line opposite vertices
- [docs/events.md](docs/events.md) — intersection lines, triple events, sort key
- [docs/alignment.md](docs/alignment.md) — y-parallel walls, alignment oracle vs zone
- [docs/event_list.md](docs/event_list.md) — combined events, initial z, VD at a slice
- [docs/sweep.md](docs/sweep.md) — reference sweep, 3D cells, matching, oracles
- [docs/incremental.md](docs/incremental.md) — local 2D updates vs recompute oracle
- [docs/robustness.md](docs/robustness.md) — simultaneous groups, zone degeneracies
- [docs/invariants.md](docs/invariants.md) — checklist from the design; boxes are checked only when code enforces them

## Layout

```text
vd3d/                 library (kernel packages must not import viz)
docs/                 human-readable notes, one concern per file
tests/unit/           exact, non-graphical tests
tests/visual/         write PNG + caption under artifacts/visual/
tests/oracles/        slow brute-force checkers
artifacts/visual/     generated review gallery (gitignored)
```
