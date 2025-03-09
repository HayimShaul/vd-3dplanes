from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
from functools import singledispatch
from typing import Tuple
import numpy as np
import intersection


# return the height of b above a in the axis direction
@singledispatch
def height(a, b, axis) -> Tuple:
    """Generic height function (default case)."""
    raise NotImplementedError(f"Height not implemented for {type(a)} and {type(b)} and {type(axis)}")

# returns the height of point above plane.
# if the point is below the plane the returned value will be negative.
@height.register
def _(point : Point3D, plane: Plane, axis):
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

# returns the height of point above segment.
# the point needs to be within the bounds og the segment
@height.register
def _(point : Point3D, seg: Segment3D, axis):
    if axis == 'y':
        if point.x < min(seg.p1.x, seg.p2.x) or point.x > max(seg.p1.x, seg.p2.x):
            return False
    else:
        raise ValueError(f"axis {axis} not supported in height(point, segment)")


    plane = Plane(seg.p1, seg.p2, seg.p1 + Point3D(0,0,1))
    return height(point, plane, 'y')

# returns the height of point above segment.
# the point needs to be within the bounds og the segment
@height.register
def _(point : Point3D, ray: Ray3D, axis):
    if axis == 'y':
        if point.x < ray.p1.x and ray.direction.x > 0:
            return False
        if point.x > ray.p1.x and ray.direction.x < 0:
            return False
    else:
        raise ValueError(f"axis {axis} not supported in height(point, ray)")


    plane = Plane(ray.p1, ray.p1 + ray.direction, ray.p1 + Point3D(0,0,1))
    return height(point, plane, 'y')



# check whether a lies on b
@singledispatch
def incident(a, b) -> Tuple:
    """Generic project function (default case)."""
    raise NotImplementedError(f"Incident not implemented for {type(a)} and {type(b)}")

@incident.register
def _(point: Point3D, plane: Plane):
    return height(point, Plane, 'z') == 0

@incident.register
def _(seg: Segment3D, plane: Plane):
    return height(seg.p1, plane, 'z') == 0 and height(seg.p2, plane, 'z') == 0

@incident.register
def _(ray: Ray3D, plane: Plane):
    return height(ray.p1, plane, 'z') == 0 and height(ray.p1 + ray.direction, plane, 'z') == 0

@incident.register
def _(line: Line3D, plane: Plane):
    return height(line.p1, plane, 'z') == 0 and height(line.p2, plane, 'z') == 0

@incident.register
def _(plane: Plane, point: Point3D):
    return incident(point, plane)

@incident.register
def _(plane: Plane, seg: Segment3D):
    return incident(seg, plane)

@incident.register
def _(plane: Plane, ray: Ray3D):
    return incident(ray, plane)

@incident.register
def _(plane: Plane, line: Line3D):
    return incident(line, plane)



# check whether a is above b on the axis direction
@singledispatch
def is_above(a, b, axis) -> Tuple:
    """Generic project function (default case)."""
    raise NotImplementedError(f"Incident not implemented for {type(a)} and {type(b)}, with {type(axis)}")

@is_above.register
def _(point: Point3D, plane: Plane, axis : str):
    return height(point, Plane, axis) > 0

@is_above.register
def _(seg: Segment3D, plane: Plane, axis : str):
    return height(seg.p1, plane, axis) > 0 and height(seg.p2, plane, axis) > 0

@is_above.register
def _(ray: Ray3D, plane: Plane, axis : str):
    if axis == 'z':
        return height(ray.p1, plane, axis) > 0 and ray.direction.z >= 0
    raise ValueError("axis must be z in is_above(Ray3D, Plane)")


# check whether a is below b on the axis direction
@singledispatch
def is_below(a, b, axis) -> Tuple:
    """Generic project function (default case)."""
    raise NotImplementedError(f"Incident not implemented for {type(a)} and {type(b)}, with {type(axis)}")

@is_below.register
def _(point: Point3D, plane: Plane, axis : str):
    return height(point, Plane, axis) < 0

@is_below.register
def _(seg: Segment3D, plane: Plane, axis : str):
    return height(seg.p1, plane, axis) < 0 and height(seg.p2, plane, axis) < 0

@is_below.register
def _(ray: Ray3D, plane: Plane, axis : str):
    if axis == 'z':
        return height(ray.p1, plane, axis) < 0 and ray.direction.z <= 0
    raise ValueError("axis must be z in is_below(Ray3D, Plane)")


# check whether a is below b on the axis direction
@singledispatch
def is_directly_above(a, plane, planes, axis) -> Tuple:
    """Generic project function (default case)."""
    raise NotImplementedError(f"Incident not implemented for {type(a)}, {type(plane)}, {type(planes)}, {type(axis)}")

@is_directly_above.register
def _(point: Point3D, plane: Plane, planes, axis: str):
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

@is_directly_above.register
def _(seg: Segment3D, plane: Plane, planes, axis: str):
    return is_directly_above(seg.p1, plane, planes, axis) and is_directly_above(seg.p2, plane, planes, axis)

@is_directly_above.register
def _(ray: Ray3D, plane: Plane, planes, axis: str):
    p2 = ray.p1 + ray.direction
    return is_directly_above(ray.p1, plane, planes, axis) and is_directly_above(p2, plane, planes, axis)


# find the element form bs that is directly above a
def find_directly_above(a, bs, axis):
    best_h = None
    best_b = None
    for b in bs:
        h = height(a, b, axis)
        if h == False:
            continue
        if h <= 0:
            continue
        if best_b == None or best_h < h:
            best_h = h
            best_b = b
    return best_b

# find the element form bs that is directly above a
def find_directly_below(a, bs, axis):
    best_h = None
    best_b = None
    for b in bs:
        h = height(a, b, axis)
        if h == False:
            continue
        if h >= 0:
            continue
        if best_b == None or best_h > h:
            best_h = h
            best_b = b
    return best_b
