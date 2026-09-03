# Vertical decomposition of planes in 3D

Sweep-based vertical decomposition of an arrangement of planes, following
[design.md](design.md).

This is an incremental implementation. Each phase is small on purpose. Do not
start the next phase until the current one has passed its human gate.

## Status

**Phase 2 — 2D line arrangement.** Unbounded lines, DCEL faces, no bounding
box in the combinatorics. Exact `Fraction` arithmetic.

Phase 0 (scaffold) and Phase 1 (geometry kernel) are in place. Do not start
Phase 3 until you have checked the Phase 2 gallery against the captions.

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
# Phase 2
pytest tests/unit/test_arrangement2d.py tests/visual/test_arrangement2d_visual.py
python -m vd3d.viz.gallery --step 2
```

Open `artifacts/visual/index.html`. Confirm: intersection marks sit on both
lines; piece counts match the titles (1 vertex / 4 rays; 3 vertices / 6 rays
+ 3 segments); outgoing half-edge numbers increase counterclockwise from +x;
the triangle arrangement is one bounded triangle plus six unbounded cells
labeled `U`; the random `n=6` picture has no leftover slivers.

Optional GUI (cycle scenes with left/right; `g` resamples; `q` quits):

```bash
python -m vd3d.viz.viewer --phase 2
python -m vd3d.viz.viewer --phase 2 --seed 42
```

Earlier phases:

```bash
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
