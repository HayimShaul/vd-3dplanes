# Reference sweep and 3D cells

Phase 8. Input is a list of 3D planes. Output is a 3D vertical
decomposition: prisms grown by sweeping the 2D VD along `z`. Exact
`Fraction` coefficients. See [conventions.md](conventions.md).

This is the **reference** algorithm ([design.md](../design.md) §19): at
every event the 2D VD is **recomputed** at `z±ε`. Incremental 2D updates
are Phase 9 and must not start until this phase is signed off.

2D packages must not import `vd3d.sweep` or `vd3d.cells3d`. Kernel
packages must not import `vd3d.viz`.

## API

```text
match_cells(vd_before, vd_after) -> tuple[CellMatch, ...]
cell_signature(vd, cell) -> CellSignature
compute_vd_around_event(P, event, events) -> vd_before, vd_after, z−, z+
start_3d_cell / continue_3d_cell / end_3d_cell
process_event(...) -> vd_after, active, SweepSnapshot
vertical_decomposition_3d(P) -> SweepResult
locate_cell3d(result, point) -> cell id or None
```

## Matching across an event

A 2D cell continues as the same 3D cell when:

1. each representative point lies in the other cell (open point-in-cell), and
2. the lower / upper supporting plane ids agree.

Far from the event this is 1-1. The event neighbourhood is the unmatched
remainder (ended on the `z−` side, started on the `z+` side). Unmatched
cells are drawn **black**; matched cells share a colour.

Between events the active map is keyed by `CellSignature` (supporting
planes, left/right boundedness, and the planes of vertices on that
cell's own vertical sides). Signatures are unique in a slice and stable
inside an open `z` interval.

`z±ε` uses `ε = 1/100`, or half the gap to the neighbouring event if
that is smaller. Two events at the same `z` raise `SimultaneousEvents`
(Phase 10).

## 3D cell lifecycle

```text
Cell3D:
    floor, ceiling     # lower / upper supporting input planes, or none
    vertical_walls     # Steiner y-parallel planes, length ≤ 4
    lower_z, upper_z   # None = ±∞
```

Every initial 2D cell starts a 3D cell with `lower_z = −∞`. At an event,
unmatched before-cells end, matched cells continue (walls may accumulate),
unmatched after-cells start. After the last event, remaining cells end
at `+∞`. Invariant: `|active 3D cells| == |current 2D cells|`.

## Tests 19–22

- **19.** One plane (`y + z = 0`): two 3D cells (below / above), no
  vertical walls, no events.
- **20.** Two intersecting planes: four 3D cells, at most one vertical
  wall each, no events.
- **21.** Three planes, one triple: only the local neighbourhood fails
  to match 1-1; far cells keep their colour.
- **22.** Alignment at `z = 2` (isolated from the fixture's simultaneous
  triples): combinatorics change locally; far cells still match.

## Global oracles

For each open interval between events, a sample `z` must have the same
signatures as `compute_vd_at_z`. Random 3D points: independently locate
the 2D cell at `p.z`, then the active 3D cell in that interval. An
exact grid in a box is a partition (interior points in one cell,
boundary points in none).

## How to review

```bash
pytest tests/unit/test_sweep.py tests/visual/test_sweep_visual.py
python -m vd3d.viz.gallery --step 8
python -m vd3d.viz.viewer --phase 8
python -m vd3d.viz.viewer --phase 8 --seed 42 --n 4
```

Confirm: matched cells share a colour and unmatched ones are black only
near the event; one plane yields two half-spaces with no vertical walls;
query points sit in the cell of their colour.

`--n` is the number of random planes. Omit it and each random scene
picks 3 or 4. `g` in the window draws a new seed and keeps `n`.

## Checked invariants

See [invariants.md](invariants.md). Phase 8 checks: far cells match
across an event; `|active| == |2D cells|` after every event; floor /
ceiling at most one; vertical walls ≤ 4; finite `lower_z < upper_z`;
mid-interval VD equals a fresh compute; 3D point location agrees with
the independent 2D-then-active path.
