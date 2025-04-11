# projection functionalities


from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
import numpy as np


# project a onto b along an axis


def project_point3D_plane(point: Point3D, plane: Plane, axis: str):
    """Project a point onto a plane along the specified axis."""
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
    elif (axis == 'y'):
        y = -(A*x+C*z+D)/B
    else:
        raise ValueError(f"unknown axis {axis}")

    return Point3D(x,y,z)


def project_segment3D_plane(seg: Segment3D, plane: Plane, axis: str):
    """Project a segment onto a plane along the specified axis."""
    p1 = project(seg.p1, plane, axis)
    p2 = project(seg.p2, plane, axis)
    return Segment3D(p1, p2)


def project_ray3D_plane(ray: Ray3D, plane: Plane, axis: str):
    """Project a ray onto a plane along the specified axis."""
    p1 = project(ray.p1, plane, axis)
    # For the second point, we'll project a point further along the ray
    p2 = project(ray.p1 + ray.direction, plane, axis)
    # Create a new ray from the projected points
    return Ray3D(p1, p2)


def project_onto_special_plane(thing, plane: str, axis: str):
    """Project onto special planes like 'xy'."""
    if plane == 'xy' and axis == 'z':
        return project(thing, Plane(Point3D(0,0,0), (0,0,1)), 'z')
    else:
        raise ValueError(f"Invalid projection onto {plane}")


def project_line_plane(line: Line3D, plane: Plane, axis: str):
    """Project a line onto a plane along the specified axis."""
    # Project two points on the line and create a new line
    p1 = project(line.p1, plane, axis)
    p2 = project(line.p2, plane, axis)
    return Line3D(p1, p2)


def project(a, b, axis):
    """Main project function that dispatches to specific implementations based on types."""
    type_a = type(a).__name__
    type_b = type(b).__name__

    if type_a == 'Point3D':
        if type_b == 'Plane':
            return project_point3D_plane(a, b, axis)
        elif type_b == 'Segment3D':
            return project_point_segment(a, b, axis)
        elif type_b == 'Line3D':
            return project_point_line(a, b, axis)
        elif type_b == 'Ray3D':
            return project_point_ray(a, b, axis)
    elif type_a == 'Segment3D' and type_b == 'Plane':
        return project_segment3D_plane(a, b, axis)
    elif type_a == 'Ray3D' and type_b == 'Plane':
        return project_ray3D_plane(a, b, axis)
    elif type_a == 'Line3D' and type_b == 'Plane':
        return project_line_plane(a, b, axis)
    elif type_b == 'str':
        return project_onto_special_plane(a, b, axis)
    else:
        raise NotImplementedError(f"Projection not implemented for {type_a} and {type_b} and {axis}")


def project_point_line(point: Point3D, line: Line3D, axis: str):
    """Project a point onto a line along the specified axis."""
    
    if axis == 'y':
        # Verify objects lie on XY plane
        if point.z != 0 or line.p1.z != 0 or line.p2.z != 0:
            raise ValueError("Objects must lie on XY plane")
        
        # Get line direction
        direction = line.direction
        
        # If line is vertical (parallel to y-axis), return point with same x as line
        if direction.x == 0:
            raise ValueError("cannot project onto a vertical line")
        
        # Find t where line.p1 + t*direction has x coordinate equal to point.x
        # line.p1.x + t*direction.x = point.x
        t = (point.x - line.p1.x) / direction.x
        
        # Get y coordinate at this t
        y = line.p1.y + t * direction.y
        
        return Point3D(point.x, y, 0)

    elif axis == 'z':
        # Get line direction
        direction = line.direction
        
        if direction.x == 0 and direction.y == 0:
            raise ValueError("cannot project onto a z-vertical line")

        if direction.x != 0:
            # Find t where line.p1 + t*direction has x coordinate equal to point.x
            # line.p1.x + t*direction.x = point.x
            t = (point.x - line.p1.x) / direction.x
        else:   
            # Find t where line.p1 + t*direction has y coordinate equal to point.y
            # line.p1.y + t*direction.y = point.y
            t = (point.y - line.p1.y) / direction.y
        
        # Get z coordinate at this t
        z = line.p1.z + t * direction.z
        
        return Point3D(point.x, point.y, z)
    else:
        raise ValueError("project_point_line only supports y or z axis")
        


def project_point_segment(point: Point3D, seg: Segment3D, axis: str):
    """Project a point onto a segment along the specified axis."""
    # Project onto infinite line first
    line = Line3D(seg.p1, seg.p2)
    projected = project_point_line(point, line, axis)
    
    # Check if projection falls within segment bounds
    if projected.x < min(seg.p1.x, seg.p2.x) or projected.x > max(seg.p1.x, seg.p2.x):
        print(f"projected.x={float(projected.x)}, seg.p1.x={float(seg.p1.x)}, seg.p2.x={float(seg.p2.x)}")
        return None
        
    return projected


def project_point_ray(point: Point3D, ray: Ray3D, axis: str):
    # Project onto infinite line first
    line = Line3D(ray.p1, ray.p1 + ray.direction)
    projected = project_point_line(point, line, axis)
    
    # Check if projection falls within ray bounds
    if ray.direction.x > 0 and projected.x < ray.p1.x:
        return None
    if ray.direction.x < 0 and projected.x > ray.p1.x:
        return None
        
    return projected

