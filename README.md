# Vertical decomposition of planes in 3D

Sweep-based vertical decomposition of an arrangement of planes, following
[design.md](design.md).

This is an incremental implementation. Each phase is small on purpose. Do not
start the next phase until the current one has passed its human gate.

## Status

**Phase 1 — geometry kernel.** Points, planes, plane–plane and three-plane
intersection, and horizontal slices. Exact `Fraction` arithmetic.

Phase 0 (scaffold) is in place. Do not start Phase 2 until you have rotated
the Phase 1 scenes in the GUI and they match the on-screen captions.

## Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Requires Python 3.11+. The virtual environment is gitignored.

## How to review a step

Every step ships three things: unit tests, inspectable 3D scenes, and docs.
A step is done only when all three hold:

1. That step's unit tests are green.
2. You open the review GUI, rotate the 3D views, and they match the captions.
3. The docs for that step list the invariants that are actually checked.

Do not start the next step if a scene and its caption disagree.

Commands for the current phase:

```bash
# Phase 1
pytest tests/unit/test_geometry.py tests/unit/test_scenes.py
python -m vd3d.viz.viewer --phase 1
python -m vd3d.viz.viewer --phase 1 --seed 42
```

A window opens on random planes, lines, and points. With no ``--seed`` the
current time is the seed (printed as `seed=...`). Pass that value to replay
the same geometry.

Drag the 3D axes to rotate. Left/right (or `n`/`p`) cycles the five scenes;
`g` draws a new random seed; `q` quits.

Confirm by rotating: a point on/above/below a plane; two planes and their
intersection line; two parallel planes with no line; three planes meeting
at one point; and a horizontal slice whose 3D trace stays in the plane.

```bash
# Phase 0 (gallery pipeline)
pytest tests/unit/test_scalar.py tests/unit/test_linalg.py tests/unit/test_architecture.py tests/visual/test_gallery_smoke.py
python -m vd3d.viz.gallery --step 0
```

Later phases use the same pattern:

```bash
pytest tests/unit/test_<module>.py
python -m vd3d.viz.viewer --phase N
```

## Conventions and invariants

- [docs/conventions.md](docs/conventions.md) — axes, "vertical", exact arithmetic, general position
- [docs/geometry.md](docs/geometry.md) — plane equation, orientation, intersections, slices
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
