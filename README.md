# Vertical Decomposition of Planes in 3D

This library computes the 3D vertical decomposition of an arrangement of planes in 3D.
It also comes with a simple visualization tool to test the vert. decom.


## Dependencies

- SymPy: For geometric computations and primitives
- NumPy: For numerical computations and array operations
- Matplotlib: For visualization
- SciPy: For Delaunay triangulation

## Installation

1. Clone the repository
2. Install the required dependencies:
```bash
pip install sympy numpy matplotlib scipy
```

## Usage

### 2D Vertical Decomposition

```python
from vd import test_vd2d

# Run 2D vertical decomposition test with random segments and rays
test_vd2d()
```

### 3D Voronoi Diagram

```python
from vd import test_vd

# Run 3D Voronoi diagram test with random planes
test_vd()
```

## Theoretic Efficiency
The library is NOT claiming to be efficient.
While the vertical decomposition of 3D planes can be computed in $O(n^3)$. This library computes it in $O(n^6)$.


## Implementation Details

The project implements vertical decomposition using the following key components:

1. **Intersection Detection**: Computes intersections between geometric primitives
2. **Projection**: Projects points and lines onto planes and axes
3. **Cell Computation**: Determines cell boundaries and vertices
4. **Visualization**: Renders cells and boundaries using matplotlib

## Visualization

The project provides multiple visualization options:
- 2D cell visualization with boundaries and intersection points
- 3D cell visualization using Poly3DCollection
- Support for both regular and slanted plane visualizations

## Cells
Cells are given as a 6-tuple: (x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil)

Each component represents a boundary of the cell:
- x_floor, x_ceil: A scalar definining the minimum and maximum x-coordinates (can be None if unbounded)
- y_floor, y_ceil: A 3DLine defining the upper and lower vertical walls bounding the cell on the y direction (can be None if unbounded)
- z_floor, z_ceil: A Plane defining the upper and lower planes bounding the cell on the z direction (can be None if unbounded)

For 2D cells, only the first four components are used (x_floor, x_ceil, y_floor, y_ceil).

## Functions

- `vd.vd2d(p, p_segs)`: Computes 2D vertical decomposition on a plane
- `vd.vd(planes)`: Computes 3D Voronoi diagram for a set of planes
- `cells.get_cell_wall_surface(cell, p)`: Computes the polygon of a face of the cell for visualization. Assumes a bounding box of (-10, -10, -10) - (10,10,10). `p` needs to be one of (x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil). Also use `get_cell_x_floor_surface(cell)`, `get_cell_x_ceil_surface(cell)`, `get_cell_y_floor_surface(cell)`, `get_cell_y_ceil_surface(cell)`
- `cells.find_cell_vertices(cell)`: Finds all vertices of a cell by computing intersections of cell boundaries
- `cells.is_point_in_cell(p, cell)`: Checks if a point strictly lies within a cell (not on boundaries)
- `cells.is_point_in_cell_or_on_boundary(p, cell)`: Checks if a point lies within a cell or on its boundaries
- `cells.is_intersecting_cell(plane, cell, endpoints)`: Determines if a plane intersects a cell by checking if vertices lie on both sides of the plane. If `endpoints` are 'None' they will be computed from `cell`


## Testing

The project includes several test functions:
- `test_vd2d()`: Tests 2D vertical decomposition
- `test_vd()`: Tests 3D vertical decomposition

## Notes

- The implementation uses rational arithmetic through SymPy for exact geometric computations
- Visualization may require conversion from symbolic to floating-point values
- Random seed can be set for reproducible results

## License

If you use this library in your research please cite it.

The library is published under the MIT licence.


## Authors

Hayim Shaul (hayim.shaul@gmail.com)
