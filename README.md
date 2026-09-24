# Vertical decomposition of planes in 3D

Sweep-based vertical decomposition of an arrangement of planes, following
[design.md](design.md). Exact `Fraction` arithmetic throughout.

## Status

Implemented end to end: events, incremental 2D updates, 3D cell lifecycle,
and point location. For debugging, `--show-sweep` steps through each event
and shows how the sweep-plane vertical decomposition evolves before and after.

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

## License

This software is provided as-is, without warranty of any kind. I would like
it to be bug-free, but I cannot guarantee that. There is no support beyond
what I can do in my free time.

If you use this project in your research, please cite it.

## Author

Hayim Shaul.
