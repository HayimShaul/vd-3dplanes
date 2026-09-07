# Incremental 2D updates

Phase 9. After Phase 8 the reference sweep is trusted: at every event
the 2D VD can be recomputed at `z±ε`. This phase replaces the `z+`
recompute with a **local** update and keeps the recompute as an oracle.

Exact `Fraction` coefficients. See [conventions.md](conventions.md).

2D packages must not import `vd3d.sweep`. Kernel packages must not
import `vd3d.viz`.

## API

```text
update_2d_decomposition(vd_before, event, planes, z+) -> VerticalDecomposition
update_for_triple_intersection(...)
update_for_vertical_alignment(...)
event_vertex_plane_keys(event, planes)
local_vertex_plane_keys(vd_before, arr_after, event, planes)
equivalent_vd(incremental, reference)
```

Dispatcher (design §12):

```text
if event is a triple:     UPDATE_FOR_TRIPLE_INTERSECTION
if event is an alignment: UPDATE_FOR_VERTICAL_ALIGNMENT
```

After every event:

```text
equivalent(incremental_vd, compute_vd_at_z(z+))
```

`PROCESS_EVENT` still recomputes `z−` and `z+`. The `z−` slice is the
input to the handler (tight local window). The incremental `z+` result
is asserted against the recomputed `z+` slice, then used for matching
and the 3D lifecycle. `vertical_decomposition_3d(..., incremental=False)`
is the Phase 8 path.

## Local region

An event vertex is identified by the sorted source-plane ids of the
lines through it (a pair for a generic vertex).

- **Triple:** the three pairwise intersections of the three planes.
- **Alignment:** the two slice-vertices of the two intersection lines.

The **x-window** is the closed interval spanned by those vertices at
`z−` and at `z+`. Every arrangement vertex in that window is local:
its ±y rays are reshot on the after-slice. Far vertices keep their
before-event hit combinatorics (which supporting plane the ray hits);
the hit is reconstructed at the after-slice geometry.

Cells are rebuilt from the resulting walls. Far combinatorics survive;
only the event strip is allowed to change.

## Tests

- Triple fixture (Test 21): incremental `z+` equals `compute_vd_at_z`
  (signatures, wall geometry, cell bounds). Unmatched cells stay in the
  event neighbourhood.
- Alignment at `z = 2` (isolated from simultaneous triples): same
  oracle check; combinatorics change locally.
- Random general-position `n ≤ 5`: after every event the incremental
  VD equals a fresh compute. The full sweep still passes the
  mid-interval oracle.
- Handlers reject the other event type.

## How to review

```bash
pytest tests/unit/test_incremental.py tests/visual/test_incremental_visual.py
python -m vd3d.viz.gallery --step 9
python -m vd3d.viz.viewer --phase 9
python -m vd3d.viz.viewer --phase 9 --seed 42 --n 4
```

Confirm: the triple and alignment before/after pairs use the same
matching colours as Phase 8 (matched share a colour, unmatched are
**black** only near the event). The incremental-vs-oracle pair has
**no black cells** — every colour agrees.

`--n` is the number of random planes. Omit it and each random scene
picks 3 or 4. `g` in the window draws a new seed and keeps `n`.

## Checked invariants

See [invariants.md](invariants.md). Phase 9 checks: after every event
the incremental 2D VD is combinatorially equal to `compute_vd_at_z(z+)`;
walls and trapezoid bounds match the recomputed slice; far cells
continue across the event.
