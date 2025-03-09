from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
from functools import singledispatch
from typing import Tuple
import numpy as np

# Define geometric objects
class Segment:
    def __init__(self, start: Point3D, end: Point3D):
        self.start = Point3D(start)
        self.end = Point3D(end)

    def __repr__(self):
        return f"Segment(start={self.start}, end={self.end})"

class Ray:
    def __init__(self, start: Point3D, direction: Point3D):
        self.start = Point3D(start)
        self.direction = Point3D(direction)

    def __repr__(self):
        return f"Ray(start={self.start}, direction={self.direction})"


class Ray:
    def __init__(self, start: Point3D, direction: Point3D):
        self.start = Point3D(start)
        self.direction = Point3D(direction)

    def __repr__(self):
        return f"Ray(start={self.start}, direction={self.direction})"


class Plane:
    def __init__(self, plane: Plane):
        self.plane = Plane(plane)

    def __repr__(self):
        return f"Plane(plane={self.plane})"


