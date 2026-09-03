# Vertical decomposition of planes in 3D

Sweep-based vertical decomposition of an arrangement of planes, following
[design.md](design.md).

This is an incremental implementation. Each phase is small on purpose. Do not
start the next phase until the current one has passed its human gate.

## Status

**Phase 0 — scaffold.** Package layout, exact `Fraction` scalars, linear-algebra
helpers, documentation skeleton, and the visual-review gallery.

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
2. You open the gallery and the figures match their captions.
3. The docs for that step list the invariants that are actually checked.

Do not start the next step if a figure and its caption disagree.

Commands for the current phase:

```bash
# Phase 0
pytest tests/unit/test_scalar.py tests/unit/test_linalg.py tests/unit/test_architecture.py tests/visual/test_gallery_smoke.py
python -m vd3d.viz.gallery --step 0
```

Then open `artifacts/visual/index.html`.

You should see a unit square with corners `(0,0)`, `(1,0)`, `(1,1)`, `(0,1)`.
The yellow box at the top of the page is the review checklist. If the square
is missing or blank, the gallery pipeline is broken.

Later phases use the same pattern:

```bash
pytest tests/unit/test_<module>.py tests/visual/test_<module>.py
python -m vd3d.viz.gallery --step N
```

## Conventions and invariants

- [docs/conventions.md](docs/conventions.md) — axes, "vertical", exact arithmetic, general position
- [docs/invariants.md](docs/invariants.md) — checklist from the design; boxes are checked only when code enforces them

## Layout

```text
vd3d/                 library (kernel packages must not import viz)
docs/                 human-readable notes, one concern per file
tests/unit/           exact, non-graphical tests
tests/visual/         write PNG + caption under artifacts/visual/
tests/oracles/        slow brute-force checkers (later phases)
artifacts/visual/     generated review gallery (gitignored)
```
