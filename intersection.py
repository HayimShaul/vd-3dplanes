from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
import numpy as np


def get_all_intersection_points(planes):
    """
    Takes a list of Planes and computes all intersection points between any 3 planes.

    Args:
        planes: List of sympy Plane objects

    Returns:
        List of Point3D objects representing intersection points
    """
    points = []
    n = len(planes)
    
    # Check all combinations of 3 planes
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                try:
                    # Get intersection point of 3 planes
                    point = planes[i].intersection(planes[j], planes[k])[0]
                    if isinstance(point, Point3D):
                        points.append(point)
                except:
                    # Skip if no intersection exists
                    continue
                    
    return points


def intersect_plane_plane(p1: Plane, p2: Plane):
    """Intersect two planes to get their intersection line."""
    # Use sympy's built-in intersection method
    intersection = p1.intersection(p2)
    
    # Return the first element which should be a Line3D
    # If no intersection exists, this will raise an exception
    if len(intersection) == 0:
        return []
    return [intersection[0]]


def intersect_line3D_plane(l: Line3D, p: Plane):
    """Intersect a line with a plane to get their intersection point."""
    # Use sympy's built-in intersection method
    intersection = l.intersection(p)
    
    # Return the first element which should be a Point3D
    # If no intersection exists, this will raise an exception
    if len(intersection) == 0:
        return []
    return [intersection[0]]


def intersect_segment3D_segment3D(s1: Segment3D, s2: Segment3D):
    l1 = Line3D(s1.p1, s1.p2)
    l2 = Line3D(s2.p1, s2.p2)
    p = intersect(l1, l2)[0]
    if p.x < min(s1.p1.x, s1.p2.x) or p.x > max(s1.p1.x, s1.p2.x):
        return []
    if p.x < min(s2.p1.x, s2.p2.x) or p.x > max(s2.p1.x, s2.p2.x):
        return []
    return [p]


def intersect_ray3D_segment3D(r: Ray3D, s: Segment3D):
    l1 = Line3D(r.p1, r.p2)
    l2 = Line3D(s.p1, s.p2)
    p = intersect(l1, l2)[0]
    if p.x < min(s.p1.x, s.p2.x) or p.x > max(s.p1.x, s.p2.x):
        return []
    if p.x < r.p1.x and r.direction.x > 0:
        return []
    if p.x > r.p1.x and r.direction.x < 0:
        return []
    return [p]


def intersect_line3D_line3D(l1: Line3D, l2: Line3D):
    """Intersect two lines to get their intersection point."""
    intersection = l1.intersection(l2)
    if len(intersection) == 0:
        return []
    return [intersection[0]]


def intersect_ray3D_ray3D(r1: Ray3D, r2: Ray3D):
    """Intersect two rays to get their intersection point."""
    l1 = Line3D(r1.p1, r1.p1 + r1.direction)
    l2 = Line3D(r2.p1, r2.p1 + r2.direction)
    p = intersect_line3D_line3D(l1, l2)
    if p == []:
        return p
    p = p[0]
    
    # Check if intersection point is in the direction of both rays
    if p.x < r1.p1.x and r1.direction.x > 0:
        return []
    if p.x > r1.p1.x and r1.direction.x < 0:
        return []
    if p.x < r2.p1.x and r2.direction.x > 0:
        return []
    if p.x > r2.p1.x and r2.direction.x < 0:
        return []
    return [p]


def intersect(a, b):
    """Main intersection function that dispatches to specific implementations based on types."""
    type_a = type(a).__name__
    type_b = type(b).__name__

    if type_a == 'Plane' and type_b == 'Plane':
        return intersect_plane_plane(a, b)
    elif type_a == 'Line3D' and type_b == 'Plane':
        return intersect_line3D_plane(a, b)
    elif type_a == 'Plane' and type_b == 'Line3D':
        return intersect_line3D_plane(b, a)
    elif type_a == 'Line3D' and type_b == 'Line3D':
        return intersect_line3D_line3D(a, b)
    elif type_a == 'Segment3D' and type_b == 'Segment3D':
        return intersect_segment3D_segment3D(a, b)
    elif type_a == 'Ray3D' and type_b == 'Segment3D':
        return intersect_ray3D_segment3D(a, b)
    elif type_a == 'Segment3D' and type_b == 'Ray3D':
        return intersect_ray3D_segment3D(b, a)
    elif type_a == 'Ray3D' and type_b == 'Ray3D':
        return intersect_ray3D_ray3D(a, b)
    else:
        raise NotImplementedError(f"Intersection not implemented for {type_a} and {type_b}")

def _get_direction_ratios(dir1, dir2):
    """Helper function to compare direction ratios."""
    # Handle zero components carefully
    ratios = []
    
    # For each component pair
    for c1, c2 in [(dir1.x, dir2.x), (dir1.y, dir2.y), (dir1.z, dir2.z)]:
        # If both components are 0, they're parallel in this dimension
        if c1 == 0 and c2 == 0:
            continue
        # If only one is 0, vectors aren't parallel
        if c1 == 0 or c2 == 0:
            return None
        # Otherwise add ratio to list
        ratios.append(c1 / c2)
    
    # If we have no ratios (all components were 0,0), vectors are parallel
    if not ratios:
        return True
        
    # Check if all ratios are equal
    return all(r == ratios[0] for r in ratios)

def parallel_seg_seg(seg1: Segment3D, seg2: Segment3D) -> bool:
    """Returns True if two segments are parallel."""
    dir1 = seg1.p2 - seg1.p1
    dir2 = seg2.p2 - seg2.p1
    result = _get_direction_ratios(dir1, dir2)
    return result is True or result is not None

def parallel_seg_line(seg: Segment3D, line: Line3D) -> bool:
    """Returns True if segment is parallel to line."""
    dir1 = seg.p2 - seg.p1
    dir2 = line.direction
    result = _get_direction_ratios(dir1, dir2)
    return result is True or result is not None

def parallel_seg_ray(seg: Segment3D, ray: Ray3D) -> bool:
    """Returns True if segment is parallel to ray."""
    dir1 = seg.p2 - seg.p1
    dir2 = ray.direction
    result = _get_direction_ratios(dir1, dir2)
    return result is True or result is not None

def parallel_line_line(line1: Line3D, line2: Line3D) -> bool:
    """Returns True if two lines are parallel."""
    dir1 = line1.direction
    dir2 = line2.direction
    result = _get_direction_ratios(dir1, dir2)
    return result is True or result is not None

def parallel_line_ray(line: Line3D, ray: Ray3D) -> bool:
    """Returns True if line is parallel to ray."""
    dir1 = line.direction
    dir2 = ray.direction
    result = _get_direction_ratios(dir1, dir2)
    return result is True or result is not None

def parallel_ray_ray(ray1: Ray3D, ray2: Ray3D) -> bool:
    """Returns True if two rays are parallel."""
    dir1 = ray1.direction
    dir2 = ray2.direction
    result = _get_direction_ratios(dir1, dir2)
    return result is True or result is not None

def parallel(a, b) -> bool:
    """Main parallel function that dispatches to specific implementations."""
    type_a = type(a).__name__
    type_b = type(b).__name__

    if type_a == 'Segment3D':
        if type_b == 'Segment3D':
            return parallel_seg_seg(a, b)
        elif type_b == 'Line3D':
            return parallel_seg_line(a, b)
        elif type_b == 'Ray3D':
            return parallel_seg_ray(a, b)
    elif type_a == 'Line3D':
        if type_b == 'Segment3D':
            return parallel_seg_line(b, a)
        elif type_b == 'Line3D':
            return parallel_line_line(a, b)
        elif type_b == 'Ray3D':
            return parallel_line_ray(a, b)
    elif type_a == 'Ray3D':
        if type_b == 'Segment3D':
            return parallel_seg_ray(b, a)
        elif type_b == 'Line3D':
            return parallel_line_ray(b, a)
        elif type_b == 'Ray3D':
            return parallel_ray_ray(a, b)
            
    raise NotImplementedError(f"Parallel not implemented for {type_a} and {type_b}")

