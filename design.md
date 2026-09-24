# Vertical Decomposition of Planes in 3D — Modular Pseudocode, Tests, and Implementation Plan

## 0. Model and assumptions

Input:

```text
P = {P1, P2, ..., Pn}
```

where every `Pi` is an infinite plane in 3D.

The sweep coordinate is `z`.

At any fixed `z = z0`, intersect every plane with the horizontal sweep plane

```text
H(z0) = {(x,y,z0)}
```

to obtain a set of 2D lines.

The 2D vertical decomposition of those lines is the decomposition maintained during the sweep.

The sweep events are:

```text
E1: triple intersections of 3 planes
E2: vertical alignments of two intersection lines
```

We assume general position initially:

```text
- no plane pair is parallel unless explicitly supported
- no 4 planes intersect in a common point
- no two independent event conditions occur at exactly the same z
- intersection lines are not degenerate
- no alignment event is tangent/degenerate
```

These assumptions should be relaxed later by grouping simultaneous events.

---

# 1. Core data structures

## 1.1 Basic geometry

```text
Plane:
    id
    equation: a*x + b*y + c*z + d = 0

Line3D:
    id
    plane_a
    plane_b
    point
    direction

Point3D:
    x, y, z

Point2D:
    x, y

Line2D:
    id
    equation / point+direction
    source_plane_id
```

---

## 1.2 2D arrangement

Use a DCEL/half-edge representation if the implementation will become substantial.

```text
Arrangement2D:
    vertices
    half_edges
    edges
    faces

ArrangementVertex:
    position
    incident_edges
    source_information

ArrangementEdge:
    curve = segment/ray/line portion
    left_face
    right_face

ArrangementFace:
    boundary
    unbounded?
```

---

## 1.3 2D vertical decomposition

Each decomposition cell should know its vertical walls and neighboring cells.

```text
VDCell2D:
    id

    left_boundary
    right_boundary

    lower_boundary
    upper_boundary

    vertical_walls[]

    neighboring_cells[]

    source_face
```

Depending on the exact definition of "vertical" in your 2D decomposition, rename these boundaries accordingly.

The important invariant is:

```text
number_of_vertical_walls(VDCell2D) <= 4
```

---

## 1.4 Sweep events

```text
Event:
    z
    type                // TRIPLE_INTERSECTION or VERTICAL_ALIGNMENT

    geometric_data

    involved_planes[]
    involved_lines[]
```

For deterministic processing:

```text
sort key = (z, event_type, stable_id)
```

Later, simultaneous events can be grouped:

```text
EventGroup:
    z
    events[]
```

---

## 1.5 3D cell

A 3D cell is constructed incrementally from the 2D cells encountered during the sweep.

```text
Cell3D:
    id

    floor              // plane or NONE
    ceiling            // plane or NONE

    vertical_walls[]   // <= 4

    lower_z
    upper_z

    cross_sections[]   // optional debugging representation
```

The intended invariant is:

```text
at most one floor
at most one ceiling
at most four vertical walls
```

---

## 1.6 Active cell record

```text
ActiveCell:
    cell3d_id

    current_2d_cell_id

    start_z

    current_floor
    current_ceiling

    current_vertical_walls[]
```

The active-cell table maps:

```text
2D cell -> currently growing 3D cell
```

---

# 2. Geometry module

This module should be implemented and tested before anything involving sweep logic.

## 2.1 Plane-plane intersection

```text
function INTERSECT_PLANES(P, Q):
    if NORMAL(P) parallel NORMAL(Q):
        return PARALLEL

    direction = CROSS(NORMAL(P), NORMAL(Q))

    point = solve_one_point_on_both_planes(P, Q)

    return Line3D(point, direction, P.id, Q.id)
```

Test:

```text
P: x = 0
Q: y = 0

=> line:
   (0, 0, t)
```

---

## 2.2 Plane-plane-plane intersection

```text
function INTERSECT_THREE_PLANES(P, Q, R):
    A =
        [Px.a Px.b Px.c
         Qx.a Qx.b Qx.c
         Rx.a Rx.b Rx.c]

    b =
        [-Px.d
         -Qx.d
         -Rx.d]

    if determinant(A) == 0:
        return NONE

    p = solve(A, b)

    return p
```

---

## 2.3 Plane with horizontal sweep plane

```text
function SLICE_PLANE_AT_Z(P, z):
    substitute z into:

        a*x + b*y + c*z + d = 0

    obtain:

        A*x + B*y + C = 0

    return Line2D(A, B, C, source_plane=P.id)
```

---

## 2.4 Build the complete slice arrangement

```text
function BUILD_SLICE_LINES(P, z):
    lines = []

    for plane in P:
        line = SLICE_PLANE_AT_Z(plane, z)

        if line is not degenerate:
            lines.append(line)

    return lines
```

```text
function BUILD_2D_ARRANGEMENT(P, z):
    lines = BUILD_SLICE_LINES(P, z)

    arrangement = BUILD_LINE_ARRANGEMENT(lines)

    return arrangement
```

---

# 3. 2D line arrangement module

This module is independent of the 3D sweep.

```text
function BUILD_LINE_ARRANGEMENT(lines):

    create one unbounded face

    for every pair (Li, Lj):
        if Li and Lj intersect:
            create arrangement vertex
            split Li and Lj at the intersection

    sort edge pieces along each line

    connect consecutive edge pieces

    construct half-edge structure

    construct faces

    return Arrangement2D
```

This should be testable using only 2D lines.

---

# 4. 2D vertical decomposition module

This module is also independent.

Input:

```text
Arrangement2D
```

Output:

```text
VD2D
```

A convenient modular decomposition is:

```text
function COMPUTE_VERTICAL_DECOMPOSITION(arrangement):

    vd = initialize_with_arrangement_faces(arrangement)

    vertices = arrangement.vertices

    for vertex in vertices:

        rays = VERTICAL_RAYS_FROM(vertex)

        for ray in rays:
            hit = FIRST_INTERSECTION_WITH_ARRANGEMENT(ray)

            INSERT_DECOMPOSITION_SEGMENT(
                vd,
                vertex,
                hit
            )

    MERGE_OR_CREATE_DECOMPOSITION_CELLS(vd)

    VERIFY_VD_INVARIANTS(vd)

    return vd
```

The geometry-specific operations are isolated:

```text
VERTICAL_RAYS_FROM(vertex)

FIRST_INTERSECTION_WITH_ARRANGEMENT(ray)

INSERT_DECOMPOSITION_SEGMENT(...)

MERGE_OR_CREATE_DECOMPOSITION_CELLS(...)
```

so the vertical-decomposition convention can be changed without touching the sweep code.

---

# 5. Zone computation module

This should be a completely independent module.

For a 2D arrangement and a query line `L`, compute the zone of `L`.

```text
function COMPUTE_ZONE(arrangement, L):

    crossings = []

    for every arrangement edge E:

        if INTERSECTS(L, E):
            p = INTERSECTION(L, E)
            crossings.append((p, E))

    sort crossings along L

    zone_cells = []
    zone_edges = []
    zone_vertices = []

    current_face = FACE_CONTAINING_POINT(
        point slightly before first crossing
    )

    append current_face to zone_cells

    for crossing in crossings:

        record crossing.edge
        record crossing.vertex if applicable

        next_face = FACE_ON_OTHER_SIDE(
            current_face,
            crossing.edge
        )

        append next_face to zone_cells

        current_face = next_face

    return Zone(
        cells = zone_cells,
        edges = zone_edges,
        vertices = zone_vertices
    )
```

Important invariant:

```text
The zone is the ordered sequence of arrangement features
crossed by L.
```

---

## 5.1 Zone test

Example arrangement:

```text
three lines forming a triangle
```

Query line:

```text
L crossing the triangle from left to right
```

Expected:

```text
unbounded face
    -> triangle boundary
    -> interior face
    -> triangle boundary
    -> unbounded face
```

Additional randomized test:

```text
generate random lines
generate random query line

compare COMPUTE_ZONE(...)
against brute-force sorting of all intersections
```

This test is particularly valuable because it isolates a difficult geometric operation.

---

# 6. Intersection-line preprocessing

First compute every pairwise intersection line.

```text
function COMPUTE_INTERSECTION_LINES(P):

    lines = []

    for i = 1 .. n:
        for j = i+1 .. n:

            result = INTERSECT_PLANES(P[i], P[j])

            if result != PARALLEL:
                lines.append(result)

    return lines
```

---

# 7. Triple-intersection event generation

```text
function GENERATE_TRIPLE_EVENTS(P):

    events = []

    for i = 1 .. n:
        for j = i+1 .. n:
            for k = j+1 .. n:

                p = INTERSECT_THREE_PLANES(
                    P[i], P[j], P[k]
                )

                if p exists:

                    events.append(
                        Event(
                            z = p.z,
                            type = TRIPLE_INTERSECTION,
                            geometric_data = p,
                            involved_planes = [i,j,k]
                        )
                    )

    return events
```

This is intentionally simple in the first implementation.

Optimization can come later.

---

# 8. Vertical-alignment events

This is the second major geometric module.

The implementation should follow the procedure you described.

## 8.1 Vertical wall associated with an intersection line

```text
function BUILD_VERTICAL_WALL(L):

    projection = PROJECT_TO_XY(L)

    W = VERTICAL_WALL_THROUGH(projection)

    return W
```

The exact definition depends on the chosen notion of vertical alignment.

It should be isolated behind this interface so that the rest of the algorithm does not depend on the representation.

---

## 8.2 Intersect planes with the wall

```text
function SLICE_PLANES_BY_WALL(P, W):

    wall_lines = []

    for plane in P:

        result = INTERSECT_PLANE_WITH_VERTICAL_WALL(
            plane,
            W
        )

        if result is not empty:
            wall_lines.append(
                Line2D(
                    geometry = result,
                    source_plane_id = plane.id
                )
            )

    return wall_lines
```

---

## 8.3 Compute the zone of an intersection line on the wall

```text
function COMPUTE_LINE_ZONE_ON_WALL(L, W, P):

    wall_lines = SLICE_PLANES_BY_WALL(P, W)

    arrangement = BUILD_2D_ARRANGEMENT(wall_lines)

    query_line = REPRESENT_L_INSIDE_WALL(L, W)

    zone = COMPUTE_ZONE(
        arrangement,
        query_line
    )

    return zone
```

---

## 8.4 Convert zone vertices to alignment events

```text
function EXTRACT_ALIGNMENT_EVENTS_FROM_ZONE(
    L,
    W,
    zone
):

    events = []

    for vertex in zone.vertices:

        candidate = INTERPRET_ZONE_VERTEX_AS_ALIGNMENT(
            vertex,
            L,
            W
        )

        if candidate exists:

            if candidate.represents_two_distinct_intersection_lines:

                z = COMPUTE_ALIGNMENT_Z(candidate)

                events.append(
                    Event(
                        z = z,
                        type = VERTICAL_ALIGNMENT,
                        geometric_data = candidate,
                        involved_lines =
                            [candidate.L1, candidate.L2]
                    )
                )

    return events
```

The validation step is important:

```text
function VALIDATE_ALIGNMENT_EVENT(candidate):

    compute the two intersection lines L1, L2

    verify:
        their projections satisfy the vertical-alignment predicate
        at z = candidate.z

    return TRUE / FALSE
```

---

## 8.5 Generate all vertical-alignment events

```text
function GENERATE_ALIGNMENT_EVENTS(P):

    intersection_lines =
        COMPUTE_INTERSECTION_LINES(P)

    events = []

    for L in intersection_lines:

        W = BUILD_VERTICAL_WALL(L)

        zone =
            COMPUTE_LINE_ZONE_ON_WALL(
                L,
                W,
                P
            )

        local_events =
            EXTRACT_ALIGNMENT_EVENTS_FROM_ZONE(
                L,
                W,
                zone
            )

        for event in local_events:
            events.append(event)

    return events
```

Because the same event may be discovered more than once:

```text
function DEDUPLICATE_EVENTS(events):

    canonical_events = hash_map()

    for e in events:

        key = CANONICAL_EVENT_KEY(e)

        canonical_events[key] = MERGE_EVENT(
            canonical_events[key],
            e
        )

    return values(canonical_events)
```

---

# 9. Complete event generation

```text
function GENERATE_ALL_EVENTS(P):

    triple_events =
        GENERATE_TRIPLE_EVENTS(P)

    alignment_events =
        GENERATE_ALIGNMENT_EVENTS(P)

    events =
        triple_events + alignment_events

    events =
        DEDUPLICATE_EVENTS(events)

    events =
        SORT_EVENTS_BY_Z(events)

    return events
```

---

# 10. Initial slice at z = -infinity

We cannot literally use `z = -infinity`.

Instead choose:

```text
z0 < minimum event z
```

with enough margin that no event occurs below it.

A robust implementation can compute:

```text
z0 = min_event_z - margin
```

where `margin` is determined from the problem scale.

For exact arithmetic, a symbolic `-infinity` state is even better.

```text
function INITIAL_SLICE(P, events):

    z0 = CHOOSE_Z_BELOW_ALL_EVENTS(events)

    arrangement =
        BUILD_2D_ARRANGEMENT(P, z0)

    vd =
        COMPUTE_VERTICAL_DECOMPOSITION(arrangement)

    return vd
```

---

# 11. Processing one event

The central sweep operation is:

```text
old_vd
    |
    | process event
    v
new_vd
```

We need to determine which 2D cells continue, which end, and which start.

---

## 11.1 Compute slice immediately before and after event

For event at `ze` choose:

```text
z_minus = NEXT_REPRESENTABLE_VALUE_BELOW(ze)
z_plus  = NEXT_REPRESENTABLE_VALUE_ABOVE(ze)
```

or symbolic perturbations:

```text
ze - epsilon
ze + epsilon
```

Then:

```text
function COMPUTE_VD_AROUND_EVENT(P, event):

    z_minus = event.z - epsilon
    z_plus  = event.z + epsilon

    vd_before =
        COMPUTE_VERTICAL_DECOMPOSITION(
            BUILD_2D_ARRANGEMENT(P, z_minus)
        )

    vd_after =
        COMPUTE_VERTICAL_DECOMPOSITION(
            BUILD_2D_ARRANGEMENT(P, z_plus)
        )

    return vd_before, vd_after
```

For efficiency, the production implementation should update locally instead of recomputing globally.

But recomputation is extremely useful for testing.

---

# 12. Event-local 2D update

A useful abstraction:

```text
function UPDATE_VD_AT_EVENT(
    vd_before,
    event
):

    affected_region =
        FIND_LOCAL_REGION_AFFECTED_BY_EVENT(
            vd_before,
            event
        )

    remove:
        edges/cells whose combinatorics disappear

    add:
        edges/cells whose combinatorics appear

    repair adjacency

    return vd_after
```

There should be separate handlers:

```text
UPDATE_FOR_TRIPLE_INTERSECTION(...)

UPDATE_FOR_VERTICAL_ALIGNMENT(...)
```

Dispatcher:

```text
function UPDATE_2D_DECOMPOSITION(vd, event):

    if event.type == TRIPLE_INTERSECTION:
        return UPDATE_FOR_TRIPLE_INTERSECTION(vd, event)

    if event.type == VERTICAL_ALIGNMENT:
        return UPDATE_FOR_VERTICAL_ALIGNMENT(vd, event)

    error("unknown event type")
```

This separation makes debugging much easier.

---

# 13. Matching 2D cells across an event

We need to decide whether a cell before the event is the same 3D cell after the event.

A robust conceptual version is:

```text
function MATCH_CELLS_ACROSS_EVENT(
    vd_before,
    vd_after,
    event
):

    unaffected_before =
        CELLS_FAR_FROM_EVENT(vd_before, event)

    unaffected_after =
        CELLS_FAR_FROM_EVENT(vd_after, event)

    matches = []

    for cell_before in unaffected_before:

        cell_after =
            FIND_SAME_REGION_AFTER_EVENT(
                cell_before,
                vd_after
            )

        if cell_after exists:
            matches.append(
                (cell_before, cell_after)
            )

    local_match =
        MATCH_AFFECTED_CELLS_TO_LOCAL_CELLS(
            vd_before,
            vd_after,
            event
        )

    matches.extend(local_match)

    return matches
```

The most practical implementation is to give every 2D cell a stable geometric signature or choose an interior sample point.

For example:

```text
function CELL_MATCH_KEY(cell):

    p = REPRESENTATIVE_POINT(cell)

    determine:
        containing arrangement face
        incident vertical walls
        adjacent plane ids

    return canonical signature
```

During development, an exact geometric point-in-cell test is preferable to relying solely on IDs.

---

# 14. Starting a new 3D cell

```text
function START_3D_CELL(
    cell2d,
    event_z,
    metadata
):

    cell3d = Cell3D()

    cell3d.floor = DETERMINE_FLOOR(cell2d)
    cell3d.ceiling = NONE

    cell3d.vertical_walls =
        EXTRACT_VERTICAL_WALLS(cell2d)

    cell3d.lower_z = event_z
    cell3d.upper_z = +infinity

    register cell3d

    return cell3d
```

---

# 15. Ending a 3D cell

```text
function END_3D_CELL(active_cell, event_z):

    cell3d =
        GET_CELL(active_cell.cell3d_id)

    cell3d.upper_z = event_z

    cell3d.ceiling =
        DETERMINE_CEILING(active_cell)

    FINALIZE_VERTICAL_WALLS(cell3d)

    VERIFY_CELL_INVARIANTS(cell3d)

    remove active_cell
```

---

# 16. Continue a 3D cell

```text
function CONTINUE_3D_CELL(
    active_cell,
    cell2d_after,
    event
):

    active_cell.current_2d_cell_id =
        cell2d_after.id

    active_cell.current_vertical_walls =
        EXTRACT_VERTICAL_WALLS(cell2d_after)

    UPDATE_CELL_BOUNDARY_INFORMATION(
        active_cell,
        event
    )
```

---

# 17. Process one sweep event

```text
function PROCESS_EVENT(
    P,
    event,
    active_cells
):

    old_vd = CURRENT_VD

    # During development:
    recomputed_before,
    recomputed_after =
        COMPUTE_VD_AROUND_EVENT(P, event)

    ASSERT(
        equivalent(
            old_vd,
            recomputed_before
        )
    )

    new_vd =
        UPDATE_2D_DECOMPOSITION(
            old_vd,
            event
        )

    # Extremely useful development check:
    ASSERT(
        equivalent(
            new_vd,
            recomputed_after
        )
    )

    matches =
        MATCH_CELLS_ACROSS_EVENT(
            old_vd,
            new_vd,
            event
        )

    matched_before = all left sides of matches
    matched_after  = all right sides of matches

    # 1. End cells that disappeared.
    for cell_before in old_vd.cells:

        if cell_before not in matched_before:

            active =
                FIND_ACTIVE_CELL(
                    cell_before.id,
                    active_cells
                )

            if active exists:
                END_3D_CELL(
                    active,
                    event.z
                )

    # 2. Continue cells that survived.
    for (cell_before, cell_after) in matches:

        active =
            FIND_ACTIVE_CELL(
                cell_before.id,
                active_cells
            )

        if active exists:

            CONTINUE_3D_CELL(
                active,
                cell_after,
                event
            )

    # 3. Start cells that did not exist before.
    for cell_after in new_vd.cells:

        if cell_after not in matched_after:

            new_cell =
                START_3D_CELL(
                    cell_after,
                    event.z,
                    event
                )

            active_cells[cell_after.id] =
                new_cell

    CURRENT_VD = new_vd
```

---

# 18. Main algorithm

```text
function VERTICAL_DECOMPOSITION_3D(P):

    # -------------------------------------------------
    # Phase 1: geometry preprocessing
    # -------------------------------------------------

    intersection_lines =
        COMPUTE_INTERSECTION_LINES(P)

    # -------------------------------------------------
    # Phase 2: event generation
    # -------------------------------------------------

    events =
        GENERATE_ALL_EVENTS(P)

    event_groups =
        GROUP_EVENTS_WITH_SAME_Z(events)

    # -------------------------------------------------
    # Phase 3: initial decomposition
    # -------------------------------------------------

    initial_vd =
        INITIAL_SLICE(P, events)

    CURRENT_VD = initial_vd

    active_cells = EMPTY_MAP()
    output_cells = []

    # Every initial 2D cell starts a 3D cell.
    z0 = CHOOSE_Z_BELOW_ALL_EVENTS(events)

    for cell2d in initial_vd.cells:

        cell3d =
            START_3D_CELL(
                cell2d,
                z0,
                metadata = NONE
            )

        active_cells[cell2d.id] = cell3d

    # -------------------------------------------------
    # Phase 4: sweep
    # -------------------------------------------------

    for group in event_groups:

        # If general position is used, group contains one event.
        # Otherwise process a simultaneous-event transaction.

        for event in group.events:

            PROCESS_EVENT(
                P,
                event,
                active_cells
            )

    # -------------------------------------------------
    # Phase 5: close cells at +infinity
    # -------------------------------------------------

    for active in active_cells:

        END_3D_CELL(
            active,
            +infinity
        )

        output_cells.append(
            GET_CELL(active.cell3d_id)
        )

    return output_cells
```

---

# 19. Recommended simplification for the first implementation

Do NOT initially implement the incremental update of the 2D decomposition.

Instead use this version:

```text
for every event e:

    z_before = e.z - epsilon
    z_after  = e.z + epsilon

    vd_before =
        COMPUTE_VERTICAL_DECOMPOSITION_AT_Z(z_before)

    vd_after =
        COMPUTE_VERTICAL_DECOMPOSITION_AT_Z(z_after)

    MATCH_CELLS(vd_before, vd_after)

    UPDATE_3D_CELLS(...)
```

Once this version passes all tests, replace:

```text
COMPUTE_VERTICAL_DECOMPOSITION_AT_Z(z_after)
```

with:

```text
UPDATE_2D_DECOMPOSITION_INCREMENTALLY(...)
```

and compare the incremental result against the brute-force result.

This gives you a very strong reference implementation.

---

# 20. Debugging representation

For every event, save a snapshot:

```text
SweepSnapshot:
    event_id
    event_z

    vd_before
    event
    vd_after

    cell_matches

    cells_started
    cells_ended
```

This makes failures reproducible.

For example:

```text
Event 17
z = 4.271

before:
    12 cells

event:
    VERTICAL_ALIGNMENT
    lines L3, L9

after:
    12 cells

continued:
    11

ended:
    1

started:
    1
```

---

# 21. Invariants to check after EVERY operation

## Geometry invariants

```text
INTERSECT_PLANES(P,Q):
    returned line lies on both planes

INTERSECT_THREE_PLANES(P,Q,R):
    returned point lies on all three planes

SLICE_PLANE_AT_Z(P,z):
    every returned 2D point maps to a 3D point lying on P
```

---

## Arrangement invariants

For every arrangement vertex:

```text
vertex lies on all incident edges
```

For every edge:

```text
edge lies on its supporting line
```

For every half-edge:

```text
twin(twin(e)) == e
```

For every face:

```text
boundary cycle is valid
```

---

## Vertical decomposition invariants

```text
Every decomposition cell has a valid boundary.

Every vertical decomposition wall belongs to
the underlying arrangement/decomposition.

No two decomposition cells overlap in their interiors.

The union of decomposition cells equals the
underlying arrangement domain.

number_of_vertical_walls(cell) <= 4
```

---

## Sweep invariants

At every open interval:

```text
CURRENT_VD == independently computed VD at a sample z
```

Every active 3D cell corresponds to exactly one active 2D cell.

```text
|active 3D cells| == |current 2D cells|
```

provided your chosen cell correspondence has this one-to-one interpretation.

---

## 3D cell invariants

At completion:

```text
cell.floor is NONE or exactly one plane
cell.ceiling is NONE or exactly one plane

len(cell.vertical_walls) <= 4
```

Additional useful invariant:

```text
floor.z < ceiling.z
```

for bounded cells.

---

# 22. Test suite

The tests should be layered.

## Level 1 — basic geometry

### Test 1. Plane-plane intersection

```text
P1: x=0
P2: y=0
```

Expected:

```text
intersection = z-axis
```

---

### Test 2. Three-plane intersection

```text
P1: x=0
P2: y=0
P3: z=1
```

Expected:

```text
(0,0,1)
```

---

### Test 3. Parallel planes

```text
x=0
x=2
```

Expected:

```text
NO INTERSECTION
```

---

### Test 4. Horizontal slice

For:

```text
P: x + y + z - 5 = 0
```

at `z=2`:

```text
x + y - 3 = 0
```

---

# 23. Level 2 — 2D arrangement tests

### Test 5. Two intersecting lines

Expected:

```text
1 vertex
4 rays / edge pieces
4 sectors
```

depending on representation.

---

### Test 6. Three concurrent lines

Expected one common vertex.

---

### Test 7. Three lines in general position

Expected:

```text
3 intersection vertices
```

and the known number of arrangement faces.

---

### Test 8. Randomized arrangement

Generate:

```text
N random nonparallel lines
```

Compare the arrangement against a simple brute-force implementation.

---

# 24. Level 3 — vertical decomposition tests

### Test 9. Single line

Expected trivial decomposition.

---

### Test 10. Two intersecting lines

Check that the vertical rays generated from the intersection vertex produce the expected decomposition.

---

### Test 11. Triangle arrangement

Use three lines forming a triangle.

Check:

```text
all decomposition cells
all vertical walls
adjacency
no overlap
full coverage
```

---

### Test 12. Random arrangement

For 10–100 random small arrangements:

```text
VD = COMPUTE_VERTICAL_DECOMPOSITION(arrangement)

assert invariants

compare cell adjacency against a brute-force geometric test
```

---

# 25. Level 4 — zone tests

This should be a dedicated test suite.

### Test 13. Simple zone

Arrangement:

```text
5 nonparallel lines
```

Query line:

```text
L
```

Compute:

```text
zone = COMPUTE_ZONE(A,L)
```

Then independently:

```text
intersections = all intersections of L with arrangement edges
sort by position along L
```

Verify that the sequence of crossed faces/edges is identical.

---

### Test 14. Zone boundary cases

Test when L:

```text
- passes far from all vertices
- passes through a vertex
- overlaps an arrangement edge
- is parallel to some arrangement edges
- is extremely close to a vertex
```

The last three are especially useful for exposing robustness problems.

---

# 26. Level 5 — event-generation tests

## Test 15. Triple intersection

Construct:

```text
P1, P2, P3
```

with known triple intersection:

```text
(x,y,z) = (1,2,3)
```

Expected:

```text
exactly one triple event at z=3
```

---

## Test 16. No triple intersection

Generate random planes in general position where no triple intersection occurs.

Expected:

```text
zero triple events
```

---

## Test 17. Known vertical alignment

Construct two pairwise intersection lines whose projections become aligned at a known `z`.

Expected:

```text
one alignment event at z = z_expected
```

This test is important because it isolates the hardest event-generation component.

---

## Test 18. Zone algorithm vs brute force

For small `n`:

```text
events_zone =
    GENERATE_ALIGNMENT_EVENTS_USING_ZONES(P)

events_brute_force =
    GENERATE_ALL_ALIGNMENT_EVENTS_BRUTE_FORCE(P)
```

Then compare the canonicalized event sets.

This is probably the single most valuable correctness test for the alignment implementation.

---

# 27. Level 6 — sweep tests

## Test 19. One plane

There should be exactly one 3D region decomposition induced by the plane, with:

```text
one ceiling/floor boundary as appropriate
zero or appropriate vertical walls
```

---

## Test 20. Two planes

Two intersecting planes should produce the expected subdivision into cells.

Check:

```text
number of cells
floor/ceiling planes
vertical walls
```

---

## Test 21. Three planes with one triple event

Use a configuration where all three planes intersect at a known point.

Sweep:

```text
below event
at event
above event
```

Verify that only the local 2D combinatorics change.

---

## Test 22. Alignment-only event

Construct a configuration with:

```text
no triple event near z0
one vertical alignment event at z0
```

Verify:

```text
VD(before) != VD(after)
```

in the expected local neighborhood, while everything outside that neighborhood remains unchanged.

---

# 28. Level 7 — differential testing

For small instances, implement two algorithms:

### Reference implementation

```text
for every interesting z interval:
    recompute everything from scratch
```

### Incremental implementation

```text
start at -infinity
process events
update locally
```

Then compare after every event:

```text
reference_vd
      ==
incremental_vd
```

and:

```text
reference_cells
      ==
incremental_cells
```

up to canonical renaming of IDs.

This should remain in the test suite permanently.

---

# 29. Level 8 — randomized property tests

Generate random planes with a fixed seed.

For each generated instance:

```text
events = GENERATE_ALL_EVENTS(P)

sort(events)

for each consecutive event pair:
    choose random z strictly between them

    vd1 =
        incremental/current decomposition

    vd2 =
        brute-force decomposition at z

    assert equivalent(vd1, vd2)
```

Additionally:

```text
for every output cell:
    assert floor/ceiling <= 1
    assert vertical_walls <= 4
```

This kind of randomized differential testing will likely find bugs that hand-designed examples miss.

---

# 30. Strong global correctness test

For every open z interval:

```text
choose z_sample

compute 2D VD at z_sample
```

For every resulting 2D cell `C`:

```text
lift C through the interval in z
```

The corresponding 3D cell should be the one tracked by the sweep.

Then verify a partition property on random 3D points:

```text
for random point p:

    determine its 3D decomposition cell
    independently determine which cell it should belong to

    assert same cell
```

For bounded test configurations, you can sample a finite bounding box.

---

# 31. Implementation plan

The safest implementation strategy is evolutionary.

## Step 1 — Robust geometric kernel

Implement:

```text
Plane
Point3D
Line3D
Point2D
Line2D

INTERSECT_PLANES
INTERSECT_THREE_PLANES
SLICE_PLANE_AT_Z
```

Tests:

```text
Tests 1–4
```

Do not continue until these are reliable.

---

## Step 2 — 2D line arrangement

Implement:

```text
BUILD_LINE_ARRANGEMENT
```

using a DCEL or another explicit planar representation.

Tests:

```text
Tests 5–8
```

At this point you have a completely independent, testable 2D arrangement package.

---

## Step 3 — 2D vertical decomposition

Implement:

```text
COMPUTE_VERTICAL_DECOMPOSITION
```

and make it operate only on `Arrangement2D`.

Tests:

```text
Tests 9–12
```

Do not involve 3D cells yet.

---

## Step 4 — Zone computation

Implement:

```text
COMPUTE_ZONE
```

as a standalone module.

Tests:

```text
Tests 13–14
```

Add randomized differential tests against a brute-force implementation.

This is worth doing before touching alignment events.

---

## Step 5 — Pairwise intersection lines and triple events

Implement:

```text
COMPUTE_INTERSECTION_LINES
GENERATE_TRIPLE_EVENTS
```

Tests:

```text
Tests 15–16
```

Now the sweep event infrastructure exists except for alignment events.

---

## Step 6 — Vertical-alignment event generator

Implement exactly the modular chain:

```text
BUILD_VERTICAL_WALL
    ->
SLICE_PLANES_BY_WALL
    ->
BUILD_2D_ARRANGEMENT
    ->
COMPUTE_ZONE
    ->
EXTRACT_ALIGNMENT_EVENTS_FROM_ZONE
```

Then:

```text
GENERATE_ALIGNMENT_EVENTS
```

Tests:

```text
Tests 17–18
```

This should be developed against a brute-force event generator.

Do not optimize it initially.

---

## Step 7 — Reference sweep

Before implementing incremental updates, implement:

```text
for each event:
    recompute VD immediately before event
    recompute VD immediately after event
    match cells
```

Then implement the 3D cell lifecycle:

```text
START_3D_CELL
CONTINUE_3D_CELL
END_3D_CELL
```

Tests:

```text
Tests 19–22
```

At this stage the complete algorithm can already be validated even though the 2D decomposition is recomputed.

---

## Step 8 — Differential-test incremental updates

Implement:

```text
UPDATE_FOR_TRIPLE_INTERSECTION
UPDATE_FOR_VERTICAL_ALIGNMENT
```

but keep:

```text
COMPUTE_VERTICAL_DECOMPOSITION_AT_Z
```

as the reference oracle.

After every event:

```text
vd_incremental =
    UPDATE_2D_DECOMPOSITION(...)

vd_reference =
    COMPUTE_VERTICAL_DECOMPOSITION_AT_Z(...)

assert equivalent(vd_incremental, vd_reference)
```

This is the crucial transition from a correct-but-slow implementation to the intended sweep algorithm.

---

## Step 9 — Optimize

Only after correctness is established, optimize:

```text
- event generation
- line arrangement construction
- zone computation
- local VD updates
- cell matching
- memory representation
```

Use profiling rather than optimizing all modules simultaneously.

---

# 32. Recommended software architecture

A clean package structure could look like:

```text
geometry/
    point3d
    line3d
    plane
    point2d
    line2d
    intersections

arrangement2d/
    arrangement
    dcel
    line_arrangement
    face
    edge

vertical_decomposition/
    vd_cell
    vertical_decomposition
    vertical_rays

zone/
    zone
    compute_zone

events/
    event
    triple_events
    vertical_wall
    alignment_events
    event_sorting

sweep/
    slice
    event_update
    cell_matching
    active_cells
    sweep

cells3d/
    cell3d
    cell_lifecycle
    validation

tests/
    geometry_tests
    arrangement_tests
    vd_tests
    zone_tests
    event_tests
    sweep_tests
    randomized_tests
    differential_tests
```

The most important architectural boundary is:

```text
2D geometry/decomposition
        ^
        |
        | used by
        |
3D sweep/event logic
```

The 2D modules should not know that a sweep exists.

---

# 33. A particularly useful testing strategy

Keep three implementations where practical:

```text
1. Simple/reference implementation
2. Incremental implementation
3. Optional brute-force geometric oracle
```

For example:

```text
reference_VD(z)
incremental_VD_after_event(e)
```

should always agree.

Similarly:

```text
zone_algorithm(L)
brute_force_zone(L)
```

should agree.

And:

```text
alignment_events_zone(P)
alignment_events_bruteforce(P)
```

should agree.

This creates a chain of independently verifiable components:

```text
                 geometry
                    |
                    v
             2D arrangement
                    |
            +-------+-------+
            |               |
            v               v
       vertical VD         zone
            |               |
            +-------+-------+
                    |
                    v
             event generator
                    |
                    v
              sweep updater
                    |
                    v
               3D cells
```

That structure means a failure at the final 3D-cell level can usually be traced back to one of a small number of modules instead of debugging the entire algorithm at once.

