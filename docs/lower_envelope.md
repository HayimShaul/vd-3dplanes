# Lower envelope and hypograph VD

Standalone construction of the 3D vertical-decomposition cells that are
**unbounded from below** (`Cell3D` with `floor=None`), without running the
full z-sweep.

## Duality

Write each non-vertical plane `a x + b y + c z + d = 0` (`c ≠ 0`) as
`z = A x + B y + C` with `A = -a/c`, `B = -b/c`, `C = -d/c`, and map it to
the dual point `(A, B, C)`. The lower envelope

```text
f(x, y) = min_i ⟨p_i, (x, y, 1)⟩
```

is determined by the **lower convex hull** of the dual points:

| Hull (lower) | Envelope |
| --- | --- |
| vertex `p_i` | face (plane `i` lowest) |
| edge `(p_i, p_j)` | edge (planes `i` and `j` meet) |
| facet `(p_i, p_j, p_k)` | vertex (triple intersection) |

## Pipeline

1. **Float envelope** (`envelope_float.py`): scipy `ConvexHull`, keep facets
   whose outward normal has negative z.
2. **Float hypograph VD** (`vd_float.py`): minimization diagram in xy →
   y-parallel trapezoids → `Cell3D(floor=None, ceiling=plane_i, …)`.
3. **Exact envelope** (`hull.py`, `envelope.py`): same duality over
   `Fraction`; enumerate lower-hull triangles.
4. **Exact hypograph VD** (`vd.py`): same trapezoid construction with exact
   arithmetic; compare to the float reference.

Vertical walls of each trapezoid include `x = const` sides and the
top/bottom diagram edges (planes parallel to the z-axis), so each prism
stays under its envelope face.

## CLI

```bash
# float envelope
python -m vd3d.lower_envelope examples/three_planes.txt --gui

# exact envelope vs float
python -m vd3d.lower_envelope examples/three_planes.txt --exact --gui

# float hypograph VD
python -m vd3d.lower_envelope examples/three_planes.txt --vd --gui

# exact hypograph VD vs float
python -m vd3d.lower_envelope examples/three_planes.txt --exact --vd --gui

python -m vd3d.lower_envelope --random-planes 4 --seed 42 --exact --vd --gui
```

## Tests

```bash
pytest tests/unit/test_lower_envelope.py tests/unit/test_architecture.py
```

`*_float.py` modules may use numpy/scipy; the rest of `vd3d/lower_envelope`
stays exact and must not import viz.
