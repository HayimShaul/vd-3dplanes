# 2D line arrangement

Phase 2. Unbounded lines in the Euclidean plane. There is **no bounding
box** in the combinatorics; a window is only for plots. Exact `Fraction`
coefficients. See [conventions.md](conventions.md).

## API

```text
intersect_lines_2d(L, M) -> Point2D | PARALLEL | COINCIDENT
build_line_arrangement(lines) -> Arrangement2D
```

`Arrangement2D` is a DCEL: vertices, undirected edge pieces, half-edges,
faces. 2D packages must not import `vd3d.events`, `vd3d.sweep`, or
`vd3d.cells3d`.

## Line–line intersection

`intersect_lines_2d` solves the 2×2 system of `a x + b y = -c`. A zero
determinant is `PARALLEL` unless a sample point of one line lies on the
other, in which case it is `COINCIDENT`. Coincident input is rejected by
`build_line_arrangement`; parallels are allowed.

## Vertices and edge pieces

Every pair of non-parallel lines contributes an intersection. Points with
the same exact coordinates are one vertex (a concurrent triple is one
vertex, not three copies).

Along each supporting line, vertices are sorted by the exact 1D parameter
`t = (-b, a) · (x, y)`. The line splits into:

- 0 vertices → 1 whole-line piece
- `k ≥ 1` vertices → 2 rays + `k − 1` interior segments

## Orientation (locked)

Outgoing half-edges around a vertex are ordered **counterclockwise**,
starting from the +x axis. Quadrants:

- 0: `dx > 0, dy ≥ 0` (includes +x)
- 1: `dx ≤ 0, dy > 0` (includes +y)
- 2: `dx < 0, dy ≤ 0` (includes -x)
- 3: `dx ≥ 0, dy < 0` (includes -y)

The sort uses only quadrant and the 2D cross product. No `atan2`.

Each half-edge has its incident face on the **left**. Bounded-face cycles
therefore walk counterclockwise. `twin(twin(e)) = e`.

## Faces

Walk `next` until a cycle closes. A face is unbounded iff its boundary
contains a ray or a whole line (not only segments).

`representative_point` is an exact interior sample: from a point on the
boundary, shoot along the left normal and take half the first positive
hit against another cycle line (or step 1 if there is no hit). Line
arrangement faces are convex, so this point lies in the face.

The empty arrangement is one unbounded face (the whole plane).

## Euler characteristic

For lines in the Euclidean plane:

```text
V - E + F = 1
```

`E` counts undirected edge pieces (rays, segments, whole lines). This is
the plane, not the sphere: there is no extra vertex at infinity in the
DCEL.

A simple arrangement of `n` lines (no parallels, no three concurrent)
has the closed counts

```text
V = n(n-1)/2
E = n^2
F = 1 + n(n+1)/2
```

## How to review

Phase 2 is 2D. The HTML gallery is the review surface:

```bash
pytest tests/unit/test_arrangement2d.py tests/visual/test_arrangement2d_visual.py
python -m vd3d.viz.gallery --step 2
```

Open `artifacts/visual/index.html`. Optional GUI (same figures, plus `g`
to resample):

```bash
python -m vd3d.viz.viewer --phase 2
python -m vd3d.viz.viewer --phase 2 --seed 42 --n 6
```

## Checked invariants

See [invariants.md](invariants.md). Phase 2 checks the four arrangement
invariants in `vd3d.arrangement2d.invariants` and
`tests/unit/test_arrangement2d.py`.
