# Zone of a query line

Phase 4. Independent of the sweep and of the 2D vertical decomposition.
Input is an `Arrangement2D` and a query line `L`. Exact `Fraction`
coefficients. See [conventions.md](conventions.md).

2D packages must not import `vd3d.events`, `vd3d.sweep`, or `vd3d.cells3d`.

## API

```text
compute_crossings(arrangement, L) -> tuple[Crossing, ...]
face_containing_point(arrangement, point) -> Face
face_on_other_side(arrangement, face, edge_id) -> Face
compute_zone(arrangement, L) -> Zone
compute_supporting_line_zone(arrangement, line_index) -> SupportingLineZone
```

A `Crossing` is an intersection of `L` with one arrangement edge, with
1D parameter `t` along `L` and `vertex_id` if the hit is a vertex.

A `Zone` is the ordered walk: faces (`cells`), crossed edges, vertices.
Through a vertex, `vertices` lists that feature and incident-edge hits
at the same `t` collapse to one crossing. An overlapping query uses the
incident faces of the coincident line.

## Direction along `L`

Crossings are sorted by increasing

```text
t = (-b, a) · (x, y) = parameter_on_line(L, point)
```

the same 1D parameter used to split arrangement lines. That is the
forward direction of `L`. It is **not** always left→right: for `y = k`
written `y - k = 0` the direction is `(-1, 0)`. Visual fixtures write
the horizontal query as `-y + k = 0` so `t = x` and numbers increase
left→right.

The walk starts at the point of `L` with parameter `t_min - 1` (strictly
before the first crossing) and steps across each crossed edge via
`FACE_ON_OTHER_SIDE`.

## General position (this phase)

`compute_zone` allows:

- `L` coincident with an arrangement line (overlap: incident faces of
  that line, vertices on the line as crossings)
- `L` through an arrangement vertex (one feature at that `t`; next face
  by sampling just after the vertex)

Parallel to some arrangement lines is allowed: those edges contribute no
crossing.

Near-misses are exact constructed lines (for example `y = 1/1000`) and
must still report the correct face sequence.

## Full zone of a supporting line

`COMPUTE_ZONE` cannot be run on a line already in the arrangement (`L`
overlaps its own edges). `compute_supporting_line_zone` instead collects
every finite vertex of every face incident to that line, then splits:

- `vertices_on_line` — lie on `L` (triple candidates later)
- `opposite_vertices` — other corners of those faces (alignment
  candidates)

Two crossing lines: only the origin, so **0 opposite vertices**.

Three lines forming the triangle `x=0`, `y=0`, `x+y=1`:

| supporting line | on the line     | opposite |
|-----------------|-----------------|----------|
| `x=0`           | `(0,0), (0,1)`  | `(1,0)`  |
| `y=0`           | `(0,0), (1,0)`  | `(0,1)`  |
| `x+y=1`         | `(1,0), (0,1)`  | `(0,0)`  |

Green points in the gallery are exactly those opposite vertices.

## How to review

```bash
pytest tests/unit/test_zone.py tests/visual/test_zone_visual.py
python -m vd3d.viz.gallery --step 4
```

Open `artifacts/visual/index.html`. Optional GUI:

```bash
python -m vd3d.viz.viewer --phase 4
python -m vd3d.viz.viewer --phase 4 --seed 42 --n 8
```

## Checked invariants

See [invariants.md](invariants.md). Phase 4 checks: crossings sorted along
`L`; each crossing lies on `L` and on its edge; consecutive zone faces
are the two sides of the recorded edge (or, through a vertex, the open
faces of `L` just before and after); the face sequence matches the
midpoint-location oracle; supporting-line opposite vertices are incident
and off `L`. Phase 10 covers overlap and through-vertex; see
[robustness.md](robustness.md).
