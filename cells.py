import math
from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
import numpy as np
import intersection
import project
import z_dist

def get_cell_wall_surface(cell, p):
    x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil = cell

    p1 = p.random_point()
    p2 = p.random_point()
    p3 = p.random_point()
    trans_x_axis = (p2 - p1)
    trans_x_axis /= math.sqrt(trans_x_axis.dot(trans_x_axis))
    trans_y_axis = (p3 - p1)
    trans_y_axis /= math.sqrt(trans_y_axis.dot(trans_y_axis))

    planes = []
    if x_floor is not None:
        planes.append(Plane(Point3D(x_floor, 0, 0), Point3D(x_floor, 1, 0), Point3D(x_floor, 0, 1)))
    if x_ceil is not None:
        planes.append(Plane(Point3D(x_ceil, 0, 0), Point3D(x_ceil, 1, 0), Point3D(x_ceil, 0, 1)))
    if y_floor is not None:
        planes.append(Plane(y_floor.p1, y_floor.p2, y_floor.p1 + Point3D(0, 0, 1)))
    if y_ceil is not None:
        planes.append(Plane(y_ceil.p1, y_ceil.p2, y_ceil.p1 + Point3D(0, 0, 1)))
    if z_floor is not None:
        planes.append(z_floor)
    if z_ceil is not None:
        planes.append(z_ceil)
    planes.append(Plane(Point3D(10,0,0), Point3D(10,10,0), Point3D(10,0,10)))
    planes.append(Plane(Point3D(-10,0,0), Point3D(-10,10,0), Point3D(-10,0,10)))
    planes.append(Plane(Point3D(0,10,0), Point3D(10,10,0), Point3D(0,10,10)))
    planes.append(Plane(Point3D(0,-10,0), Point3D(10,-10,0), Point3D(0,-10,10)))
    planes.append(Plane(Point3D(0,0,10), Point3D(10,0,10), Point3D(0,10,10)))
    planes.append(Plane(Point3D(0,0,-10), Point3D(10,0,-10), Point3D(0,-10,-10)))

    points = []
    for i in range(len(planes)):
        p1 = planes[i]
        if p == p1:
            continue
        for j in range(i+1, len(planes)):
            p2 = planes[j]
            if p == p2:
                continue
            line = intersection.intersect(p, p1)
            if len(line) == 0 or type(line[0]).__name__ != "Line3D":
                continue
            point = intersection.intersect(line[0], p2)
            if len(point) == 0 or type(point[0]).__name__ != "Point3D":
                continue
            pnt = point[0]
            # if is_point_in_cell_or_on_boundary(pnt, cell):
            if is_point_in_cell_or_on_boundary(pnt, cell) and abs(pnt.x) <= 10 and abs(pnt.y) <= 10 and abs(pnt.z) <= 10:
                points.append(pnt)


    # Remove duplicates
    unique_points = []
    for point in points:
        is_duplicate = False
        for existing in unique_points:
            if point == existing:
                is_duplicate = True
                break
        if not is_duplicate:
            unique_points.append(point)
    points = unique_points


    if len(points) < 3:
        return []

    center = Point3D(0, 0, 0)
    for p in points:
        center += p
    center /= len(points)

    # Sort points clockwise around the center
    def clockwise_sort_key(point):
        point -= center
        trans_x = point.x*trans_x_axis.x + point.y*trans_x_axis.y + point.z*trans_x_axis.z
        trans_y = point.x*trans_y_axis.x + point.y*trans_y_axis.y + point.z*trans_y_axis.z

        return math.atan2(trans_y, trans_x)
    
    sorted_points = sorted(points, key=clockwise_sort_key)
    return sorted_points


def get_cell_x_floor_surface(cell):
    if cell[0] is None:
        return []
    return get_cell_wall_surface(cell, Plane(Point3D(cell[0], 0, 0), Point3D(cell[0], 1, 0), Point3D(cell[0], 0, 1)))

def get_cell_x_ceil_surface(cell):
    if cell[1] is None:
        return []
    return get_cell_wall_surface(cell, Plane(Point3D(cell[1], 0, 0), Point3D(cell[1], 1, 0), Point3D(cell[1], 0, 1)))

def get_cell_y_floor_surface(cell):
    if cell[2] is None:
        return []
    return get_cell_wall_surface(cell, Plane(cell[2].p1, cell[2].p2, cell[2].p1 + Point3D(0, 0, 1)))

def get_cell_y_ceil_surface(cell):
    if cell[3] is None:
        return []
    return get_cell_wall_surface(cell, Plane(cell[3].p1, cell[3].p2, cell[3].p1 + Point3D(0, 0, 1)))




# def get_cell_y_floor_surface(cell):
#     return get_cell_y_surface(cell, cell[2])

# def get_cell_y_ceil_surface(cell):
#     return get_cell_y_surface(cell, cell[3])
    



def find_cell_vertices(cell) -> list[Point3D]:
    """
    Find all vertices of a cell defined by (x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil).
    Returns list of Point3D objects representing vertices.
    
    Handles infinite extents:
    - None for x_floor/x_ceil indicates infinite x extent
    - None for y_floor/y_ceil indicates infinite y extent
    - None for z_floor/z_ceil indicates infinite z extent
    """
    x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil = cell
    vertices = []
    
    # If all bounds in any dimension are None, no vertices exist
    if (x_floor is None and x_ceil is None) or \
       (y_floor is None and y_ceil is None) or \
       (z_floor is None and z_ceil is None):
        print("cell is infinite")
        return []
    
    # Create list of x-values to check
    x_values = []
    if x_floor is not None:
        x_values.append(x_floor)
    if x_ceil is not None:
        x_values.append(x_ceil)
        
    # Create list of y-bounds to check
    y_bounds = []
    if y_floor is not None:
        y_bounds.append(y_floor)
    if y_ceil is not None:
        y_bounds.append(y_ceil)
        
    # Create list of z-bounds to check
    z_bounds = []
    if z_floor is not None:
        z_bounds.append(z_floor)
    if z_ceil is not None:
        z_bounds.append(z_ceil)
    
    # Find intersections at x boundaries if they exist
    for x in x_values:
        print("x")
        point1 = Point3D(x, 0, 0)
        for y_bound in y_bounds:
            print("y")
            y_bound = project.project(y_bound, "xy", 'z')
            point2 = project.project(point1, y_bound, 'y')
            for z_bound in z_bounds:
                print("z")
                point3 = project.project(point2, z_bound, 'z')
            
                vertices.append(point3)
    
    return vertices

def is_point_in_cell(p: Point3D, cell) -> bool:
    """
    Check if point p lies within cell.
    Cell is defined as (x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil).
    
    Returns True if point lies within all bounds of the cell:
    - x_floor < p.x < x_ceil (if bounds exist)
    - y_floor < p.y < y_ceil (using height function)
    - z_floor < p.z < z_ceil (using height function)
    """
    x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil = cell
    
    # Check x bounds
    if x_floor is not None and p.x <= x_floor:
        return False
    if x_ceil is not None and p.x >= x_ceil:
        return False
        
    # Check y bounds
    if y_floor is not None:
        h_floor = z_dist.height(p, y_floor, 'y')
        if h_floor is False or h_floor <= 0:  # point is below or at y_floor
            return False
    if y_ceil is not None:
        h_ceil = z_dist.height(p, y_ceil, 'y')
        if h_ceil is False or h_ceil >= 0:  # point is above or at y_ceil
            return False
            
    # Check z bounds
    if z_floor is not None:
        h_floor = z_dist.height(p, z_floor, 'z')
        if h_floor is False or h_floor <= 0:  # point is below or at z_floor
            return False
    if z_ceil is not None:
        h_ceil = z_dist.height(p, z_ceil, 'z')
        if h_ceil is False or h_ceil >= 0:  # point is above or at z_ceil
            return False
            
    return True

def is_point_in_cell_or_on_boundary(p: Point3D, cell) -> bool:
    """
    Check if point p lies within cell.
    Cell is defined as (x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil).
    
    Returns True if point lies within all bounds of the cell:
    - x_floor < p.x < x_ceil (if bounds exist)
    - y_floor < p.y < y_ceil (using height function)
    - z_floor < p.z < z_ceil (using height function)
    """
    x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil = cell
    
    # Check x bounds
    if x_floor is not None and p.x < x_floor:
        return False
    if x_ceil is not None and p.x > x_ceil:
        return False
        
    # Check y bounds
    if y_floor is not None:
        h_floor = z_dist.height(p, y_floor, 'y')
        if h_floor is False or h_floor < 0:  # point is below or at y_floor
            return False
    if y_ceil is not None:
        h_ceil = z_dist.height(p, y_ceil, 'y')
        if h_ceil is False or h_ceil > 0:  # point is above or at y_ceil
            return False
            
    # Check z bounds
    if z_floor is not None:
        h_floor = z_dist.height(p, z_floor, 'z')
        if h_floor is False or h_floor < 0:  # point is below or at z_floor
            return False
    if z_ceil is not None:
        h_ceil = z_dist.height(p, z_ceil, 'z')
        if h_ceil is False or h_ceil > 0:  # point is above or at z_ceil
            return False
            
    return True
