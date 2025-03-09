from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
from functools import singledispatch
from typing import Tuple

@singledispatch
def break_element(a, v, axis) -> Tuple:
    """Generic break_element function (default case)."""
    raise NotImplementedError(f"break_element not implemented for {type(a)}, {type(v)}, {type(axis)}")

@break_element.register
def _(seg: Segment3D, v, axis: str):
    if str != 'y':
        raise ValueError("break_element is implemented only when breaking along y-axis")

    if v < min(seg.p1.y, seg.p2.y) or v < max(seg.p1.y, seg.p2.y):
        raise ValueError("break_element v is not in y-range of segment")

    # Find point p on segment with y = v
    direction = seg.p2 - seg.p1
    t = (v - seg.p1.y) / direction.y
    p = seg.p1 + t * direction

    # Break segment into two segments at point p
    seg1 = Segment3D(seg.p1, p)
    seg2 = Segment3D(p, seg.p2)
    return [seg1, seg2]


@break_element.register
def _(ray: Ray3D, v, axis: str):
    if str != 'y':
        raise ValueError("break_element is implemented only when breaking along y-axis")

    # For a ray, we only need to check if v is in the correct direction
    if ray.direction.y > 0 and v < ray.p1.y:
        raise ValueError("break_element v is below ray starting point")
    if ray.direction.y < 0 and v > ray.p1.y:
        raise ValueError("break_element v is above ray starting point")

    # Find point p on ray with y = v
    t = (v - ray.p1.y) / ray.direction.y
    p = ray.p1 + t * ray.direction

    # Break ray into segment and new ray at point p
    seg = Segment3D(ray.p1, p)
    new_ray = Ray3D(p, ray.direction)
    return [seg, new_ray]





@singledispatch
def endpoints(a) -> Tuple:
    """Generic endpoints function (default case)."""
    raise NotImplementedError(f"endpoints not implemented for {type(a)}")

@endpoints.register
def _(seg: Segment3D):
    return [seg.p1, seg.p2]

@endpoints.register 
def _(ray: Ray3D):
    return [ray.p1]
