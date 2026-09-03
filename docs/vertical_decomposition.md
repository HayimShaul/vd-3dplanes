# 2D vertical decomposition

Phase 3. Input is an `Arrangement2D`. Output is a trapezoidal map whose
walls are parallel to the **y-axis** (`x = const`). Exact `Fraction`
coefficients. See [conventions.md](conventions.md).

2D packages must not import `vd3d.sweep` or `vd3d.cells3d`.

## API

```text
vertical_rays_from(vertex) -> (+y ray, -y ray)
first_hit(arrangement, ray) -> Hit | UNBOUNDED
insert_decomposition_segment(walls, vertex, ray, hit)
compute_vertical_decomposition(arrangement) -> VerticalDecomposition
verify_vd_invariants(vd)
```

A `VDCell2D` is a (possibly unbounded, possibly degenerate) trapezoid:

- `left_x` / `right_x`: a vertical wall `x = const`, or unbounded
- `lower_line` / `upper_line`: a supporting line of the arrangement, or unbounded
- `vertical_walls`: ids of the positive-length walls on those sides (`<= 4`,
  and in this construction at most 2)
- `neighbors`: cells that share a positive-length side
- `source_face`: the arrangement face that contains the cell
- `representative`: an exact interior sample

## Vertical rays and first hit

From each arrangement vertex, shoot two rays with `x` fixed: `+y` and `-y`.

`FIRST_HIT` is the closest intersection with an arrangement line that does
**not** contain the origin. Lines through the vertex are incident, not
obstacles. A vertical line at a different `x` never meets the ray.

A miss is `UNBOUNDED`: the wall continues to `y = ±∞`. A hit may be an
arrangement vertex or an interior point of an edge (a Steiner point).

## Trapezoid cells

Wall `x`-coordinates are the vertex `x` values together with any vertical
input line. Those `x` values split the plane into open vertical strips.
Inside a strip the non-vertical lines do not cross, so they sort by `y`
into bands. Each band is a candidate trapezoid.

A ray wall only occupies the y-interval from the vertex to the first hit.
If a vertex `x` does not carry a wall through some band, adjacent strip
pieces in that band are merged (`MERGE_OR_CREATE_DECOMPOSITION_CELLS`).

A cell with a zero-length vertical side (two supporting lines meet at that
`x`) is a triangle / unbounded wedge. That side is not stored as a wall.

## Invariants

Enforced by `verify_vd_invariants` and `tests/oracles/vd2d.py`:

- each cell has a valid boundary (`left_x < right_x`, lower below upper)
- `number_of_vertical_walls(cell) <= 4`
- interiors are disjoint (open point-in-cell)
- the union covers the plane (grid oracle: every sample is in exactly one
  open cell, or on a wall/line)
- adjacency is symmetric
- each cell representative lies in its `source_face`

The triangle `x=0`, `y=0`, `x+y=1` is the golden fixture
(`tests/fixtures/triangle_vd.json`): 9 cells, 5 walls. The picture
`cells_triangle` is the definition of correct; the JSON was frozen from it.

## How to review

```bash
pytest tests/unit/test_vertical_decomposition.py tests/visual/test_vd2d_visual.py
python -m vd3d.viz.gallery --step 3
```

Open `artifacts/visual/index.html`. Optional GUI:

```bash
python -m vd3d.viz.viewer --phase 3
python -m vd3d.viz.viewer --phase 3 --seed 42
```

## Checked invariants

See [invariants.md](invariants.md). Phase 3 checks the five vertical-decomposition
invariants in `vd3d.vertical_decomposition.invariants` and
`tests/unit/test_vertical_decomposition.py`.
