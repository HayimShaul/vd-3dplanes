# Invariants

Copied from [design.md](../design.md) §21. A box is checked only when code
enforces that invariant (a test or an explicit `verify_*` call).

Phase 1 implements the geometry kernel. Arrangement, VD, and sweep
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

- [ ] Every arrangement vertex lies on all of its incident edges
- [ ] Every edge lies on its supporting line
- [ ] `twin(twin(e)) == e`
- [ ] Every face has a valid boundary cycle

## Vertical decomposition invariants

- [ ] Every decomposition cell has a valid boundary
- [ ] Every vertical wall belongs to the underlying arrangement / decomposition
- [ ] No two decomposition cells overlap in their interiors
- [ ] The union of decomposition cells equals the underlying arrangement domain
- [ ] `number_of_vertical_walls(cell) <= 4`

## Sweep invariants

- [ ] At every open `z` interval, `CURRENT_VD` equals an independently computed VD at a sample `z`
- [ ] `|active 3D cells| == |current 2D cells|`

## 3D cell invariants

- [ ] `cell.floor` is none or exactly one plane
- [ ] `cell.ceiling` is none or exactly one plane
- [ ] `len(cell.vertical_walls) <= 4`
- [ ] For bounded cells, floor is below ceiling
