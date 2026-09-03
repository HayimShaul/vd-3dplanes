# Invariants

Copied from [design.md](../design.md) §21. A box is checked only when code
enforces that invariant (a test or an explicit `verify_*` call).

Phase 1 implements the geometry kernel. Phase 2 implements the 2D line
arrangement. Phase 3 implements the 2D vertical decomposition. Sweep
invariants are still unchecked.

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

## Sweep invariants

- [ ] At every open `z` interval, `CURRENT_VD` equals an independently computed VD at a sample `z`
- [ ] `|active 3D cells| == |current 2D cells|`

## 3D cell invariants

- [ ] `cell.floor` is none or exactly one plane
- [ ] `cell.ceiling` is none or exactly one plane
- [ ] `len(cell.vertical_walls) <= 4`
- [ ] For bounded cells, floor is below ceiling
