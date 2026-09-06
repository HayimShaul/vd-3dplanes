# Event list and initial slice

Phase 7. Input is a list of 3D planes. Output is the combined sweep
event list and the 2D vertical decomposition of the slice at a chosen
`z`. Exact `Fraction` coefficients. See [conventions.md](conventions.md).

This is still not the sweep. It only builds the event list and
recomputes `VD(arrangement(slices))` from scratch at any height.

2D packages must not import `vd3d.events`. Kernel packages must not
import `vd3d.viz`.

## API

```text
generate_all_events(P) -> tuple[Event, ...]
choose_z_below_all_events(events, margin=1) -> Scalar
build_2d_arrangement_at_z(P, z) -> Arrangement2D
compute_vd_at_z(P, z) -> VerticalDecomposition
initial_slice(P, events=None) -> VerticalDecomposition
verify_all_events(events)
```

## Combined events

`GENERATE_ALL_EVENTS` is

```text
triples + alignments
→ DEDUPLICATE_EVENTS
→ sort by (z, type, stable_id)
```

Triple events come from [events.md](events.md). Alignment events come
from [alignment.md](alignment.md). The same geometric event found twice
(once from each wall) collapses to one canonical key.

On a general-position instance no two events share a `z`. The sort key
is still unique even when that assumption fails (type and `stable_id`
break ties).

## A z below every event

We cannot slice at `z = -∞`. `CHOOSE_Z_BELOW_ALL_EVENTS` returns

```text
min(event.z) - margin
```

with exact `margin` (default `1`). That height is strictly below every
event. With no events the slice combinatorics never change, so the
conventional height is `0`.

`INITIAL_SLICE` is `compute_vd_at_z(P, choose_z_below_all_events(events))`.

## VD at a given z

```text
compute_vd_at_z(P, z) = VD(BUILD_LINE_ARRANGEMENT(BUILD_SLICE_LINES(P, z)))
```

Each slice line keeps `source_plane_id` and a dense `id` in plane order
(horizontal planes are dropped). The 2D VD is the Phase 3 trapezoid
decomposition: y-parallel walls, at most four vertical walls per cell.

Just below and just above a triple, three 2D lines form a small
triangle (3 vertices). At the event they are concurrent (1 vertex).
That is the local change the sweep will process; this phase only
recomputes both slices independently.

## How to review

```bash
pytest tests/unit/test_event_list.py tests/visual/test_event_list_visual.py
python -m vd3d.viz.gallery --step 7
python -m vd3d.viz.viewer --phase 7
python -m vd3d.viz.viewer --phase 7 --seed 42 --n 4
```

Confirm: the timeline is a 1D `z`-axis with **red triples** and **green
alignments** in increasing `z`; the Test 15 before/after pair shows a
small triangle on each side and you can point to the flip through
`(1, 2)` — nothing else in the combinatorics changes.

`--n` is the number of random planes. Omit it and each random scene
picks 3 or 4. `g` in the window draws a new seed and keeps `n`.

## Checked invariants

See [invariants.md](invariants.md). Phase 7 checks: the combined list
equals triples ∪ alignments after dedup; keys are unique and sorted;
general-position fixtures have distinct event `z`; the chosen initial
`z` is strictly below every event; just below / at / just above Test 15
the vertex counts are 3 / 1 / 3.
