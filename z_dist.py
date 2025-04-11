from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
import numpy as np
import intersection
import project

def height_point_plane(point: Point3D, plane: Plane, axis: str):
    """Returns the height of point above plane."""
    eq = plane.equation().as_coefficients_dict().items()
    coef = {str(k): v for k, v in eq}
    A = coef.get('x',0)
    B = coef.get('y',0)
    C = coef.get('z',0)
    D = coef.get('1',0)
    # Ax + By + Cz + 1 = 0
    x = point.x
    y = point.y
    z = point.z

    if (axis == 'z'):
        z = -(A*x+B*y+D)/C
        return point.z - z
    elif (axis == 'y'):
        y = -(A*x+C*z+D)/B
        return point.y - y
    else:
        raise ValueError(f"unknown axis {axis}")

def height_point_segment(point: Point3D, seg: Segment3D, axis: str):
    """Returns the height of point above segment."""
    if axis == 'y':
        if point.x < min(seg.p1.x, seg.p2.x) or point.x > max(seg.p1.x, seg.p2.x):
            return False
    else:
        raise ValueError(f"axis {axis} not supported in height(point, segment)")

    line = Plane(seg.p1, seg.p2)
    return height(point, line, 'y')

def height_point_ray(point: Point3D, ray: Ray3D, axis: str):
    """Returns the height of point above ray."""
    if axis == 'y':
        if point.x < ray.p1.x and ray.direction.x > 0:
            return False
        if point.x > ray.p1.x and ray.direction.x < 0:
            return False
    else:
        raise ValueError(f"axis {axis} not supported in height(point, ray)")

    line = Line3D(ray.p1, ray.p1 + ray.direction)
    return height(point, line, 'y')

def height_point_line(point: Point3D, line: Line3D, axis: str):
    """Returns the height of point above line."""
    
    if axis == 'z':
        # If line is vertical (parallel to y-axis), height is undefined
        if line.direction.x == 0 and line.direction.y == 0:
            raise ValueError("cannot compute height of point above a z-vertical line")
        
        # Find t where line.p1 + t*direction has x coordinate equal to point.x
        # line.p1.x + t*direction.x = point.x
        if line.direction.x != 0:
            t = (point.x - line.p1.x) / line.direction.x
        else:
            t = (point.y - line.p1.y) / line.direction.y
        
        # Get z coordinate at this t
        z = line.p1.z + t * line.direction.z
        
        # Return height (difference in z coordinates)
        return point.z - z
    elif axis == 'y':
        # If line is vertical (parallel to y-axis), height is undefined
        if line.direction.x == 0 and line.direction.z == 0:
            raise ValueError("cannot compute height of point above a y-vertical line")
        
        # Find t where line.p1 + t*direction has x coordinate equal to point.x
        # line.p1.x + t*direction.x = point.x
        if line.direction.x != 0:
            t = (point.x - line.p1.x) / line.direction.x
        else:
            t = (point.z - line.p1.z) / line.direction.z
        
        # Get y coordinate at this t
        y = line.p1.y + t * line.direction.y
        
        # Return height (difference in y coordinates)
        return point.y - y
    else:
        raise ValueError(f"axis {axis} not supported in height(point, line)")
    

def height_ray_plane(ray: Ray3D, plane: Plane, axis: str):
    """Returns the height of ray above plane.
    Returns False if ray intersects plane (not strictly above/below)."""
    if axis != 'z':
        raise ValueError(f"axis {axis} not supported in height(ray, plane)")
    
    # Get height of starting point
    h1 = height_point_plane(ray.p1, plane, axis)
    
    # Return height at starting point
    return h1

def height_segment_plane(seg: Segment3D, plane: Plane, axis: str):
    """Returns the height of segment above plane."""
    if axis != 'z':
        raise ValueError(f"axis {axis} not supported in height(segment, plane)")
    
    # Get height of both endpoints
    p = (seg.p1 + seg.p2) / 2
    h = height_point_plane(p, plane, axis)
    
    return h

def height(a, b, axis):
    """Main height function that dispatches to specific implementations."""
    type_a = type(a).__name__
    type_b = type(b).__name__

    if type_a == 'Point3D':
        if type_b == 'Plane':
            return height_point_plane(a, b, axis)
        elif type_b == 'Segment3D':
            return height_point_segment(a, b, axis)
        elif type_b == 'Ray3D':
            return height_point_ray(a, b, axis)
        elif type_b == 'Line3D':
            return height_point_line(a, b, axis)
    elif type_a == 'Ray3D' and type_b == 'Plane':
        return height_ray_plane(a, b, axis)
    elif type_a == 'Segment3D' and type_b == 'Plane':
        return height_segment_plane(a, b, axis)
    raise NotImplementedError(f"Height not implemented for {type_a} and {type_b} and {axis}")

def incident_point_plane(point: Point3D, plane: Plane):
    return height(point, plane, 'z') == 0

def incident_segment_plane(seg: Segment3D, plane: Plane):
    return height(seg.p1, plane, 'z') == 0 and height(seg.p2, plane, 'z') == 0

def incident_ray_plane(ray: Ray3D, plane: Plane):
    return height(ray.p1, plane, 'z') == 0 and height(ray.p1 + ray.direction, plane, 'z') == 0

def incident_line_plane(line: Line3D, plane: Plane):
    return height(line.p1, plane, 'z') == 0 and height(line.p2, plane, 'z') == 0

def incident_point_line(point: Point3D, line: Line3D):
    """Returns True if point lies on line."""
    # SymPy's Line3D has a contains method that checks if a point lies on the line
    return line.contains(point)

def incident_point_segment(point: Point3D, seg: Segment3D):
    """Returns True if point lies on segment."""
    # First check if point lies on the infinite line containing the segment
    line = Line3D(seg.p1, seg.p2)
    if not incident_point_line(point, line):
        return False
        
    # Then check if point lies within the segment bounds
    # We can use x-coordinate for this check
    return (min(seg.p1.x, seg.p2.x) <= point.x <= max(seg.p1.x, seg.p2.x))

def incident_point_ray(point: Point3D, ray: Ray3D):
    """Returns True if point lies on ray."""
    # First check if point lies on the infinite line containing the ray
    line = Line3D(ray.p1, ray.p1 + ray.direction)
    if not incident_point_line(point, line):
        return False
        
    # Then check if point lies in the direction of the ray
    # Get vector from ray start to point
    point_vector = point - ray.p1
    
    # Check if point_vector points in same direction as ray
    # For each non-zero component of ray direction, check if point_vector component has same sign
    if ray.direction.x != 0:
        if (point_vector.x > 0) != (ray.direction.x > 0):
            return False
    if ray.direction.y != 0:
        if (point_vector.y > 0) != (ray.direction.y > 0):
            return False
    if ray.direction.z != 0:
        if (point_vector.z > 0) != (ray.direction.z > 0):
            return False
            
    return True

def incident(a, b) -> bool:
    """Main incident function that dispatches to specific implementations."""
    type_a = type(a).__name__
    type_b = type(b).__name__

    if type_a == 'Point3D':
        if type_b == 'Plane':
            return incident_point_plane(a, b)
        elif type_b == 'Segment3D':
            return incident_point_segment(a, b)
        elif type_b == 'Line3D':
            return incident_point_line(a, b)
        elif type_b == 'Ray3D':
            return incident_point_ray(a, b)
    elif type_a == 'Segment3D' and type_b == 'Plane':
        return incident_segment_plane(a, b)
    elif type_a == 'Ray3D' and type_b == 'Plane':
        return incident_ray_plane(a, b)
    elif type_a == 'Line3D' and type_b == 'Plane':
        return incident_line_plane(a, b)
    elif type_a == 'Plane':
        return incident(b, a)  # Swap arguments for reverse cases
    elif type_b == 'Point3D':  # Handle reverse cases
        if type_a == 'Segment3D':
            return incident_point_segment(b, a)
        elif type_a == 'Line3D':
            return incident_point_line(b, a)
        elif type_a == 'Ray3D':
            return incident_point_ray(b, a)
    raise NotImplementedError(f"Incident not implemented for {type_a} and {type_b}")

def is_above_point_plane(point: Point3D, plane: Plane, axis: str):
    return height(point, plane, axis) > 0

def is_above_segment_plane(seg: Segment3D, plane: Plane, axis: str):
    return height(seg.p1, plane, axis) > 0 and height(seg.p2, plane, axis) > 0

def is_above_ray_plane(ray: Ray3D, plane: Plane, axis: str):
    if axis == 'z':
        return height(ray.p1, plane, axis) > 0 and ray.direction.z >= 0
    raise ValueError("axis must be z in is_above(Ray3D, Plane)")

def is_above(a, b, axis) -> bool:
    """Main is_above function that dispatches to specific implementations."""
    type_a = type(a).__name__
    type_b = type(b).__name__

    if type_b == 'Plane':
        if type_a == 'Point3D':
            return is_above_point_plane(a, b, axis)
        elif type_a == 'Segment3D':
            return is_above_segment_plane(a, b, axis)
        elif type_a == 'Ray3D':
            return is_above_ray_plane(a, b, axis)
    raise NotImplementedError(f"is_above not implemented for {type_a} and {type_b} with {axis}")

def is_below_point_plane(point: Point3D, plane: Plane, axis: str):
    return height(point, plane, axis) < 0

def is_below_segment_plane(seg: Segment3D, plane: Plane, axis: str):
    return height(seg.p1, plane, axis) < 0 and height(seg.p2, plane, axis) < 0

def is_below_ray_plane(ray: Ray3D, plane: Plane, axis: str):
    if axis == 'z':
        return height(ray.p1, plane, axis) < 0 and ray.direction.z <= 0
    raise ValueError("axis must be z in is_below(Ray3D, Plane)")

def is_below(a, b, axis) -> bool:
    """Main is_below function that dispatches to specific implementations."""
    type_a = type(a).__name__
    type_b = type(b).__name__

    if type_b == 'Plane':
        if type_a == 'Point3D':
            return is_below_point_plane(a, b, axis)
        elif type_a == 'Segment3D':
            return is_below_segment_plane(a, b, axis)
        elif type_a == 'Ray3D':
            return is_below_ray_plane(a, b, axis)
    raise NotImplementedError(f"is_below not implemented for {type_a} and {type_b} with {axis}")

def is_directly_above_point(point: Point3D, plane: Plane, planes, axis: str):
    if not axis == 'z':
        raise ValueError("axis must be z in is_directly_above")
    h = height(point, plane, axis)
    if h <= 0:
        return False
    for p in planes:
        h2 = height(point, p, axis)
        if h2 > 0 and h2 < h:
            return False
    return True

def is_directly_above_segment(seg: Segment3D, plane: Plane, planes, axis: str):
    return is_directly_above_point(seg.p1, plane, planes, axis) and is_directly_above_point(seg.p2, plane, planes, axis)

def is_directly_above_ray(ray: Ray3D, plane: Plane, planes, axis: str):
    p2 = ray.p1 + ray.direction
    return is_directly_above_point(ray.p1, plane, planes, axis) and is_directly_above_point(p2, plane, planes, axis)

def is_directly_above(a, plane, planes, axis) -> bool:
    """Main is_directly_above function that dispatches to specific implementations."""
    type_a = type(a).__name__

    if type_a == 'Point3D':
        return is_directly_above_point(a, plane, planes, axis)
    elif type_a == 'Segment3D':
        return is_directly_above_segment(a, plane, planes, axis)
    elif type_a == 'Ray3D':
        return is_directly_above_ray(a, plane, planes, axis)
    raise NotImplementedError(f"is_directly_above not implemented for {type_a}")

def find_directly_above(a, bs, axis):
    best_h = None
    best_b = None
    for b in bs:
        h = height(a, b, axis)
        if h == False:
            continue
        if h >= 0:
            continue
        if best_b == None or best_h < h:
            best_h = h
            best_b = b
    return best_b

def find_directly_below(a, bs, axis):
    best_h = None
    best_b = None
    for b in bs:
        h = height(a, b, axis)
        if h == False:
            continue
        if h <= 0:
            continue
        if best_b == None or best_h > h:
            best_h = h
            best_b = b
    return best_b
