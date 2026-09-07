# Vertical-alignment events

Phase 6. Two intersection lines `L1`, `L2` **align** at `z` when their
slice vertices share an `x` and the open vertical segment between them
hits no other arrangement feature. Exact `Fraction` coefficients. See
[conventions.md](conventions.md).

The brute-force pair enumerator (`enumerate_alignment_pairs`) is the
**definition of correctness**. The wall-zone algorithm
(`generate_alignment_events`) must match it. If they disagree, debug
here; do not touch the sweep.

## API

```text
x_of_line(L, z) / y_of_line(L, z) / point_on_line_at_z(L, z)
build_vertical_wall(L) -> Plane          # x = a z + b, no y term
slice_planes_by_wall(planes, W, L) -> wall-frame Line2Ds
compute_line_zone_on_wall(L, W, planes) -> WallLineZone
extract_alignment_events_from_zone(...)
validate_alignment_event(event, planes, lines)
enumerate_alignment_pairs(planes)        # oracle
generate_alignment_events(planes)        # zone algorithm
deduplicate_events(events)
```

## Wall frame `(y, z)`

`L` must have `dz ≠ 0` (no vertical intersection line). Then

```text
x(z) = (dx/dz) z + (x0 - z0 dx/dz)
y(z) = (dy/dz) z + (y0 - z0 dy/dz)
```

`BUILD_VERTICAL_WALL(L)` is the unique plane parallel to the **y-axis**
that contains `L`:

```text
x = a z + b
```

i.e. `x - a z - b = 0` (coefficient of `y` is 0). This is *not* a 3D
cell wall (those are parallel to `z`). It is the locus of 2D vertical
walls `x = const` swept by `L` as `z` varies.

Coordinates on the wall (so a 2D-vertical alignment is vertical on the plot):

```text
Point2D.x = 3D z     (sweep, horizontal axis)
Point2D.y = 3D y     (2D vertical, vertical axis)
lift(z, y) = (a z + b, y, z)
```

A visible alignment is an open **vertical** segment at fixed `z` from `L`
to the other vertex. A slanted wall-edge that merely connects that
vertex to a *different* point of `L` does not block the event.

`SLICE_PLANES_BY_WALL` substitutes `x = a z + b` into each input plane.
`L` is inserted as one wall-line. The two source planes of `L` are
omitted (they are coincident with `L`). Planes parallel to the wall are
skipped. Traces coincident with an already-included wall-line are
dropped (Phase 10: four planes through one point).

## Oracle (Step 6.2)

For every pair of intersection lines:

1. Solve `x1(z) = x2(z)`. Parallel walls (`same a`, different `b`) never
   align. The same wall (`same a` and `b`) is not a discrete event.
2. Drop the pair if `L1` meets `L2` (a triple, or worse).
3. Keep the pair only if it is **visible**: at that `z` the two
   slice-vertices have the same `x` and the open vertical segment
   between them misses every other vertex and every slice edge.

## Zone extraction (Step 6.4)

On the wall arrangement, take the **supporting-line zone** of `L`
(Phase 4.4). Opposite vertices are alignment candidates; vertices on
`L` are triples. Each opposite vertex lifts to a 3D point on the wall.
The other intersection line through that point is `L2`. Keep the
candidate only if `validate_alignment_event` agrees with the oracle
predicate. The same event is found from both walls; `deduplicate_events`
keeps one, keyed by `(type, z, sorted line ids)`.

## Design tests

**Test 17.** Planes

```text
x + y - z = 0
x - y - z = 0
x + y + z = 7
x - y + z = 1
```

give `L12 = (t, 0, t)` and `L34 = (4-s, 3, s)`. They align at `z = 2`,
vertices `(2, 0)` and `(2, 3)`. A parallel-wall pair never aligns. A
fifth/sixth plane that plants a vertex on the open segment makes that
pair invisible.

**Test 18.** On the hand fixtures and random `n ≤ 6`, zone events and
oracle events have the same canonical keys.

## How to review

```bash
pytest tests/unit/test_alignment.py tests/visual/test_alignment_visual.py
python -m vd3d.viz.gallery --step 6
python -m vd3d.viz.viewer --phase 6
python -m vd3d.viz.viewer --phase 6 --seed 42 --n 4
```

Confirm: `L` lies in a y-independent wall; at `z=2` two vertices share a
dashed vertical, and they do not at `z=3`; the wall picture is a 2D line
arrangement in `(z, y)` with `L` bold red; green points sit on a vertical
(same `z`) to `L` and match the oracle table; orange points are blocked.

## Checked invariants

See [invariants.md](invariants.md). Phase 6 checks: sample points of `L`
lie on its wall; the wall has no `y` term; Test 17 / visibility /
never-align; zone keys equal oracle keys.
