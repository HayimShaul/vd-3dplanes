# Invariants

Copied from [design.md](../design.md) §21. A box is checked only when code
enforces that invariant (a test or an explicit `verify_*` call).

Phase 1 implements the geometry kernel. Phase 2 implements the 2D line
arrangement. Phase 3 implements the 2D vertical decomposition. Phase 4
implements the zone of a query line. Phase 5 implements pairwise
intersection lines and triple events. Phase 6 implements vertical
alignment events. Phase 7 implements the combined event list and the
2D VD of a slice at a given `z`. Phase 8 implements the reference
`z`-sweep and 3D cells. Phase 9 implements incremental 2D updates
checked against that recompute.

## Phase 0 discipline

- [x] Kernel scalars are exact `Fraction` values (`as_scalar` rejects floats)
- [x] 2x2 / 3x3 determinants and solves stay exact
- [x] Float conversion lives only in `vd3d.viz`
- [x] 2D packages do not import sweep / 3D-cell code
- [x] Kernel packages do not import `matplotlib` or `vd3d.viz`

## Geometry invariants

- [x] `INTERSECT_PLANES(P, Q)`: the returned line lies on both planes (`line_lies_on_both_planes`, `tests/unit/test_geometry.py`)
- [x] `INTERSECT_THREE_PLANES(P, Q, R)`: the returned point lies on all three planes (`point_lies_on_planes`)
- [x] `SLICE_PLANE_AT_Z(P, z)`: every returned 2D point lifts to a 3D point on `P` (`slice_lifts_to_plane`)

## Arrangement invariants

- [x] Every arrangement vertex lies on all of its incident edges (`vertex_lies_on_incident_edges`, `tests/unit/test_arrangement2d.py`)
- [x] Every edge lies on its supporting line (`edge_lies_on_supporting_line`)
- [x] `twin(twin(e)) == e` (`twin_involution`)
- [x] Every face has a valid boundary cycle (`face_cycles_valid`; Euler `V - E + F = 1` for the plane; each `representative_point` lies in exactly one face)

## Vertical decomposition invariants

- [x] Every decomposition cell has a valid boundary (`cell_has_valid_boundary`, `tests/unit/test_vertical_decomposition.py`)
- [x] Every vertical wall belongs to the underlying arrangement / decomposition (`wall_belongs_to_decomposition`)
- [x] No two decomposition cells overlap in their interiors (`interiors_disjoint_at_representatives`; grid oracle in `tests/oracles/vd2d.py`)
- [x] The union of decomposition cells equals the underlying arrangement domain (`assert_grid_partition`, `face_representatives_covered`)
- [x] `number_of_vertical_walls(cell) <= 4` (`vertical_walls_at_most_four`)

## Zone invariants

- [x] Crossings of `L` are sorted by the 1D parameter along `L` (`crossings_sorted`, `tests/unit/test_zone.py`)
- [x] Each crossing lies on `L` and on the recorded arrangement edge (`crossings_lie_on_query_and_edge`)
- [x] Consecutive zone faces are the two sides of the recorded edge (`zone_faces_match_crossings`)
- [x] The zone face sequence matches independent midpoint location on `L` (`tests/oracles/zone.py`)
- [x] Supporting-line opposite vertices are vertices of incident faces and do not lie on `L` (`supporting_split_matches_line`)
- [ ] Query through a vertex (deferred; `QueryThroughVertex`)
- [ ] Query overlapping an arrangement edge (deferred; `QueryOverlapsArrangement`)

## Event invariants

- [x] `COMPUTE_INTERSECTION_LINES`: line count equals the number of non-parallel pairs (`verify_intersection_lines`, `tests/unit/test_events.py`)
- [x] Each intersection line lies on both source planes (`each_line_lies_on_source_planes`)
- [x] `n` planes with no parallel pair yield `n(n-1)/2` lines
- [x] Test 15: planes through `(1, 2, 3)` emit exactly one triple event at `z = 3`
- [x] Test 16: a parallel family emits zero triples
- [x] Every triple event point lies on its three planes and `event.z == point.z` (`triple_event_point_on_planes`)
- [x] Triple events are sorted by unique `(z, type, stable_id)` keys
- [x] `BUILD_VERTICAL_WALL(L)` contains `L` and has no `y` term (`tests/unit/test_alignment.py`)
- [x] Test 17: one visible alignment at `z = 2`; never-align and blocked-visibility cases
- [x] Alignment events have the same `x`, distinct `y`, and `event.z` equal to both points' `z`
- [x] Test 18: `generate_alignment_events` keys equal the pair-enumeration oracle (`tests/oracles/alignment.py`)
- [x] `GENERATE_ALL_EVENTS` is triples ∪ alignments, deduped, sorted by unique `(z, type, stable_id)` (`verify_all_events`, `tests/unit/test_event_list.py`)
- [x] Combined keys equal the union of the triple enumerator and the alignment oracle (`tests/oracles/event_list.py`)
- [x] No two events share `z` on general-position fixtures (`events_have_unique_z`)
- [x] `CHOOSE_Z_BELOW_ALL_EVENTS` is strictly below every event `z` (`z_strictly_below_all_events`)
- [x] `COMPUTE_VD_AT_Z` just below / at / just above Test 15: 3 / 1 / 3 arrangement vertices (triangle vs concurrent)

## Sweep invariants

- [x] At every open `z` interval, `CURRENT_VD` equals an independently computed VD at a sample `z` (`assert_mid_interval_matches_recompute`, `tests/oracles/sweep.py`)
- [x] `|active 3D cells| == |current 2D cells|` (`verify_active_against_vd`, `tests/unit/test_sweep.py`)
- [x] After every event, the incremental 2D VD equals `compute_vd_at_z(z+)` (`equivalent_vd`, `assert_equivalent_vd`, `tests/unit/test_incremental.py`)
- [x] Incremental walls and trapezoid bounds match the recomputed slice (`wall_keys`, `cell_bound_keys`)
- [x] Unmatched cells across a triple or alignment stay in the event neighbourhood

## 3D cell invariants

- [x] `cell.floor` is none or exactly one plane (`floor_at_most_one`, `verify_cell3d`)
- [x] `cell.ceiling` is none or exactly one plane (`ceiling_at_most_one`)
- [x] `len(cell.vertical_walls) <= 4` (`vertical_walls_at_most_four`)
- [x] For bounded cells, `lower_z < upper_z` (`z_extent_ordered`; floor/ceiling are the y-supporting planes, z-extent is `lower_z`/`upper_z`)
