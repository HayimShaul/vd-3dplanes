# Robustness (simultaneous groups and degeneracies)

Phase 10. After Phase 8–9 the reference sweep and incremental 2D updates
are trusted on general-position instances. This phase relaxes the
assumptions that were deferred in [conventions.md](conventions.md).

Exact `Fraction` coefficients. 2D packages must not import `vd3d.sweep`.
Kernel packages must not import `vd3d.viz`.

Performance (design §31 Step 9) and a numerical port of
`old-vd-3d-planes` stay out of scope: that library uses a different cell
geometry, so counts need not match 1-1.

## API

```text
GROUP_EVENTS_WITH_SAME_Z(events) -> EventGroup[]
process_event_group(P, group, ...)
update_2d_for_event_group(vd_before, group, P, z+)
z_before_after(event, events)   # gap to other heights, not other events
compute_zone(A, L)              # through-vertex and overlap
coincident_line_indices(A, L)
collapse_crossings(crossings)
```

`vertical_decomposition_3d` processes each `EventGroup` as **one**
matching of `z−` / `z+`. Pass `require_general_position=True` to keep
the old rejection (`SimultaneousEvents`).

## Zone through a vertex

A query through an arrangement vertex hits several incident edges at the
same parameter `t`. `compute_crossings` still records every hit.
`compute_zone` collapses them to one feature, records the vertex, and
locates the next face by sampling `L` just after that `t`. Wedges that
only touch `L` at the point are not in the zone.

The triangle plus `y = x`: unbounded → interior → unbounded, with a
vertex mark at the origin.

## Overlapping query

If `L` coincides with an arrangement line, crossings are the vertices of
that line (not a continuum of edge points). Zone faces are the faces
incident to the line, in order of first appearance along `L`. That set
equals `compute_supporting_line_zone`. Opposite orientation (`-L`) is
the same line.

## Simultaneous event groups

`z±ε` uses the gap to the nearest **different** event height. All events
at one `z` share one before/after pair. Incremental update reshoots the
union of those events' x-windows and is still checked against
`compute_vd_at_z(z+)`.

The alignment fixture (Test 17/22) is the hand example: triples and the
alignment share heights, and the full sweep now runs.

## Other non-general-position fixtures

- **Four planes through `(1, 2, 3)`:** four triples at `z = 3`, one group.
- **Vertical input planes** (`c = 0`): the slice line does not move with
  `z`. Two vertical planes meet in a line parallel to the z-axis; that
  line still has a well-defined `x(z)` (constant). Horizontal
  intersection lines (`dz = 0`) remain skipped for alignment.

Coincident input planes are still treated as `PARALLEL` and skipped.
Coincident 2D arrangement lines are still rejected by
`build_line_arrangement`.

## How to review

```bash
pytest tests/unit/test_robustness.py tests/unit/test_zone.py tests/visual/test_robustness_visual.py
python -m vd3d.viz.gallery --step 10
python -m vd3d.viz.viewer --phase 10
python -m vd3d.viz.viewer --phase 10 --seed 42 --n 4
```

Confirm: `y=x` through the origin marks the vertex and still paints the
triangle; the overlapping `y=0` query paints every face that touches
that side; the simultaneous-group before/after pair is one picture for
several events at that `z`.

## Checked invariants

See [invariants.md](invariants.md). Phase 10 checks: through-vertex face
sequence matches midpoint location on unique `t`; overlap faces equal
the supporting-line incident set; a same-`z` group is one snapshot;
incremental group update equals `compute_vd_at_z(z+)`; four concurrent
triples and vertical input planes pass the mid-interval and partition
oracles.
