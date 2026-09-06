# Geometry kernel

Phase 1. All types use exact `Fraction` coefficients. See
[conventions.md](conventions.md) for axes and the meaning of "vertical".

## Equation form

A plane is

```text
a x + b y + c z + d = 0
```

`Plane.eval(p)` is the residual `a x + b y + c z + d`. The sign is the
orientation:

- `eval > 0` — the side the normal points toward ("above" the plane)
- `eval = 0` — on the plane
- `eval < 0` — the opposite side ("below")

Coefficients are stored as given. They are not normalized. The zero plane
`a = b = c = 0` is rejected.

A 2D line is

```text
a x + b y + c = 0
```

`a = b = 0` is degenerate and is rejected by `Line2D`. `slice_plane_at_z`
returns `None` in that case instead of constructing a line.

## Intersections

`intersect_planes(P, Q)`:

- direction is `cross(normal(P), normal(Q))`
- if that is the zero vector, the result is `PARALLEL` (includes coincident
  planes and opposite-normal copies of the same plane)
- otherwise a `Line3D` through one exact point on both planes

Invariant (enforced): sample points of the line lie on both planes.

`intersect_three_planes(P, Q, R)` solves the 3×3 system with rows
`(a, b, c)` and right-hand side `(-d, -d, -d)`. A zero determinant yields
`None` (no unique point: parallel pencil, two parallels, …).

Invariant (enforced): the returned point lies on all three planes.

## Horizontal slice

`slice_plane_at_z(P, z)` substitutes `z` into the plane equation:

```text
a x + b y + (c z + d) = 0
```

A horizontal input plane (`a = b = 0`) produces no 2D line.
`build_slice_lines` drops those.

Invariant (enforced): every sample point of the 2D line lifts to a 3D point
on `P` at that `z`.

## How to review

Do not use the static HTML gallery for Phase 1. Rotate the scenes:

```bash
python -m vd3d.viz.viewer --phase 1
python -m vd3d.viz.viewer --phase 1 --seed 42
```

`--n` is accepted (number of planes) but Phase 1 scenes have a fixed plane
count, so it does not change them.

The viewer samples random planes, lines, and points. With no `--seed` it uses
`time.time_ns()` and prints `seed=...`. Re-run with that seed to get the same
figures. `g` in the window draws a new seed.

Drag a 3D view to rotate it. Left/right (or `n`/`p`) changes scene; `q` quits.
Each scene's caption is at the bottom of the window and says what you must see.

## Checked invariants

See [invariants.md](invariants.md). Phase 1 checks the three geometry
invariants in `vd3d.geometry.invariants` and in `tests/unit/test_geometry.py`.
