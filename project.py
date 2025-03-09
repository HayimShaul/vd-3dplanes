# projection functionalities


from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
from functools import singledispatch
from typing import Tuple
import numpy as np


# project a onto b along an axis


# Define a generic function
@singledispatch
def project(a, b, axis) -> Tuple:
    """Generic project function (default case)."""
    raise NotImplementedError(f"Intersection not implemented for {type(a)} and {type(b)}  and  {type(axis)}")


# project a point onto a plane along the z axis
@project.register
def _(point: Point3D, plane: Plane, axis: str):
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



# project a point onto a plane along the z axis
@project.register
def _(seg: Segment3D, plane: Plane, axis: str):
    p1 = project(seg.p1, plane, axis)
    p2 = project(seg.p2, plane, axis)
    return Segment3D(p1, p2)



# project a point onto a plane along the z axis
@project.register
def _(thing, plane: str, axis: str):
    if plane == 'xy' and axis == 'z':
        return project(thing, Plane(Point3D(0,0,0), (0,0,1)), 'z')
    else:
        raise("Invalid projection onto {plane}")


# project a ray onto a plane along the z axis
@project.register
def _(ray: Ray3D, plane: Plane, axis: str):
    p1 = project(ray.p1, plane, axis)
    # For the second point, we'll project a point further along the ray
    p2 = project(ray.p1 + ray.direction, plane, axis)
    # Create a new ray from the projected points
    return Ray3D(p1, direction=(p2 - p1))

