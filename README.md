# Vertical decomposition of planes in 3D

Sweep-based vertical decomposition of an arrangement of planes, following
[design.md](design.md).

This is an incremental implementation. Each phase is small on purpose. Do not
start the next phase until the current one has passed its human gate.

## Status

**Phase 3 — 2D vertical decomposition.** Trapezoids with y-parallel walls
(`x = const`) on an unbounded line arrangement. Exact `Fraction` arithmetic.

Phase 0 (scaffold), Phase 1 (geometry kernel), and Phase 2 (2D arrangement)
are in place. Do not start Phase 4 until you have checked the Phase 3 gallery
against the captions.

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
# Phase 3
pytest tests/unit/test_vertical_decomposition.py tests/visual/test_vd2d_visual.py
python -m vd3d.viz.gallery --step 3
```

Open `artifacts/visual/index.html`. Confirm: cyan ±y rays from vertices either
hit the first obstacle (red mark) or reach the window edge, and never cross a
line without a hit; dashed red walls are strictly vertical; the two-line
picture has 6 trapezoids; the triangle picture has 9 cells (bounded triangle
plus split outer cells) with walls at `x=0` and `x=1`; random pictures have
no leftover slivers.

Optional GUI (cycle scenes with left/right; `g` resamples; `q` quits):

```bash
python -m vd3d.viz.viewer --phase 3
python -m vd3d.viz.viewer --phase 3 --seed 42
```

Earlier phases:

```bash
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
