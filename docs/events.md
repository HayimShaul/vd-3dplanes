# Pairwise lines and triple events

Phase 5. Input is a list of 3D planes. Output is every pairwise
intersection line (skipping parallels) and every triple-intersection
sweep event. Exact `Fraction` coefficients. See
[conventions.md](conventions.md).

2D packages must not import `vd3d.events`. Kernel packages must not
import `vd3d.viz`.

## API

```text
compute_intersection_lines(planes) -> tuple[Line3D, ...]
generate_triple_events(planes) -> tuple[Event, ...]
Event(z, type, geometric_data, plane_ids, sort_key=(z, type, stable_id))
EventType.TRIPLE_INTERSECTION | VERTICAL_ALIGNMENT
verify_intersection_lines(planes, lines)
verify_triple_events(planes, events)
```

`VERTICAL_ALIGNMENT` is generated in Phase 6; see [alignment.md](alignment.md).
The combined list and the slice VD live in Phase 7; see
[event_list.md](event_list.md).

## Intersection lines

`COMPUTE_INTERSECTION_LINES` walks every pair `i < j` and calls
`intersect_planes`. Parallel (and coincident) pairs are skipped.

For `n` planes with no parallel pair the count is

```text
n(n-1)/2
```

Each returned `Line3D` stores `plane_a` / `plane_b` (source ids) and a
dense `id` in pair order. Sample points of the line lie on both planes.

## Triple events

`GENERATE_TRIPLE_EVENTS` walks every triple `i < j < k` and calls
`intersect_three_planes`. A singular 3×3 (parallel pencil, two
parallels, …) produces no event.

An `Event` for a unique meeting point `p` is

```text
z = p.z
type = TRIPLE_INTERSECTION
geometric_data = p
plane_ids = sorted source ids
stable_id = "triple:{id}:{id}:{id}"
sort_key = (z, type, stable_id)
```

Events are returned already sorted by `sort_key`.

For `n` planes in general position (no parallels, every triple
nonsingular) the count is

```text
n(n-1)(n-2)/6
```

## Design tests

**Test 15.** Planes

```text
x + z = 4
y + z = 5
x + y + z = 6
```

meet at `(1, 2, 3)`. Exactly one triple event, at `z = 3`.

**Test 16.** A parallel family (same normal, shifted offsets) yields
zero intersection lines and zero triples.

## How to review

Phase 5 is 3D. Rotate the scenes:

```bash
pytest tests/unit/test_events.py tests/visual/test_events_visual.py
python -m vd3d.viz.gallery --step 5
python -m vd3d.viz.viewer --phase 5
python -m vd3d.viz.viewer --phase 5 --seed 42 --n 4
```

`--n` is the number of planes. Omit it and each random scene picks 3 or
4. `g` in the window draws a new seed and keeps `n`.

Drag a 3D view to rotate it. Confirm: three planes have three crease
lines; the Test 15 figure has one fat marker at `(1, 2, 3)` and one
z-axis tick at height 3; the parallel family has neither lines nor a
marker.

## Checked invariants

See [invariants.md](invariants.md). Phase 5 checks: line count equals
the number of non-parallel pairs; each line lies on both source planes;
each triple point lies on its three planes with `event.z == point.z`;
events are sorted by unique `(z, type, stable_id)` keys.
