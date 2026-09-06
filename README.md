# Vertical decomposition of planes in 3D

Sweep-based vertical decomposition of an arrangement of planes, following
[design.md](design.md).

This is an incremental implementation. Each phase is small on purpose. Do not
start the next phase until the current one has passed its human gate.

## Status

**Phase 7 — event list and initial slice.** Triples and alignments are
merged, deduplicated, and sorted by `(z, type, id)`. The 2D VD at any
`z` is recomputed from the slice arrangement. Exact `Fraction`
arithmetic.

Phases 0–6 are in place. Do not start Phase 8 until you have checked the
Phase 7 figures against the captions (timeline colors and the local
triangle flip).

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Requires Python 3.11+. The virtual environment is gitignored.

## How to review a step

Every step ships three things: unit tests, inspectable figures, and docs.
A step is done only when all three hold:

1. That step's unit tests are green.
2. You open the review surface and the figures match the captions.
3. The docs for that step list the invariants that are actually checked.

Do not start the next step if a figure and its caption disagree.

Commands for the current phase:

```bash
# Phase 7
pytest tests/unit/test_event_list.py tests/visual/test_event_list_visual.py
python -m vd3d.viz.gallery --step 7
python -m vd3d.viz.viewer --phase 7
python -m vd3d.viz.viewer --phase 7 --seed 42 --n 4
```

Confirm: the timeline is a 1D **z**-axis with **red triples** and **green
alignments** in increasing `z`; the Test 15 before/after pair has a small
triangle on each side that flips through `(1, 2)` and nothing else in the
combinatorics changes.

`--n` sets the number of random planes (phases 1, 5–7) or lines (phases 2–4).
Omit it and each random scene picks a small count. `g` resamples with a new
seed and keeps `n`.

Earlier phases:

```bash
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
