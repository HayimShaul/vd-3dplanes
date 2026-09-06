# Conventions

These rules are locked. Later phases must not invent a different meaning of
"vertical", a different scalar type, or a different sweep axis.

## Exact arithmetic

Every coordinate and coefficient in the kernel is a `fractions.Fraction`
(`vd3d.geometry.Scalar`).

- `vd3d.geometry.as_scalar` accepts `int`, `Fraction`, or a ratio string (`"3/2"`).
- It rejects `float` and `bool`.
- Convert to `float` only through `vd3d.viz.convert.to_float`, and only when drawing.

## Axes

- The horizontal plane is `xy`. Height is `z`.
- The sweep coordinate is `z`.
- At a fixed `z = z0`, each input plane is sliced by the horizontal plane
  `H(z0) = {(x, y, z0)}` to obtain a 2D line.

## What "vertical" means

- **2D vertical** means parallel to the **y-axis** (`x = const`).
  A 2D vertical-decomposition cell is a trapezoid with at most two vertical
  walls and at most two supporting-line sides.
- **3D vertical wall** means a plane **parallel to the z-axis**.
  A 3D cell is a prism: at most one floor plane, at most one ceiling plane,
  and at most four vertical walls.

## Events

- **Triple event:** three planes meet. In the `z`-slice, three 2D lines become concurrent.
- **Alignment event:** two arrangement vertices get the same `x` at that `z`
  (their 2D vertical walls coincide). Predicate: `x_L1(z) = x_L2(z)` for
  intersection lines `L1` and `L2`. Visible if and only if the open vertical
  segment between those two points hits no other arrangement feature.

## General position (first implementation)

Assume, until a later robustness phase:

- no four planes meet in a common point
- no two independent events share exactly the same `z`
- no vertical input plane (`c = 0` in `ax + by + cz + d = 0`)
- no vertical intersection line (`dx = dy = 0`)
- no plane pair is parallel unless a test explicitly asks for that case

## Module boundary

2D packages (`arrangement2d`, `vertical_decomposition`, `zone`) must not import
`vd3d.events`, `vd3d.sweep`, or `vd3d.cells3d`. Kernel packages must not import
`vd3d.viz` or `matplotlib`.
