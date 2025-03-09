from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
from functools import singledispatch
from typing import Tuple
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


# Define a generic function
@singledispatch
def intersect(a, b) -> Tuple:
    """Generic project function (default case)."""
    raise NotImplementedError(f"Intersection not implemented for {type(a)} and {type(b)}")


@intersect.register
def _(p1: Plane, p2: Plane):
    """Intersect two planes to get their intersection line."""
    # Use sympy's built-in intersection method
    intersection = p1.intersection(p2)
    
    # Return the first element which should be a Line3D
    # If no intersection exists, this will raise an exception
    return intersection[0]

@intersect.register
def _(l: Line3D, p: Plane):
    """Intersect a line with a plane to get their intersection point."""
    # Use sympy's built-in intersection method
    intersection = l.intersection(p)
    
    # Return the first element which should be a Point3D
    # If no intersection exists, this will raise an exception
    return intersection[0]


@intersect.register
def _(p: Plane, l: Line3D):
    return intersect(l, p)


@intersect.register
def _(s1: Segment3D, s2: Segment3D):
    l1 = Line3D(s1.p1, s1.p2)
    l2 = Line3D(s2.p1, s2.p2)
    p = intersect(l1, l2)
    if p.x < min(s1.p1.x, s1.p2.x) or p.x > max(s1.p1.x, s1.p2.x):
        return []
    if p.x < min(s2.p1.x, s2.p2.x) or p.x > max(s2.p1.x, s2.p2.x):
        return []
    return [p]

@intersect.register
def _(r: Ray3D, s: Segment3D):
    l1 = Line3D(r.p1, r.p2)
    l2 = Line3D(s.p1, s.p2)
    p = intersect(l1, l2)
    if p.x < min(s.p1.x, s.p2.x) or p.x > max(s.p1.x, s.p2.x):
        return []
    if p.x < r.p1.x and r.direction.x > 0:
        return []
    if p.x > r.p1.x and r.direction.x < 0:
        return []
    return [p]

@intersect.register
def _(s: Segment3D, r: Ray3D):
    return intersect(r, s)

