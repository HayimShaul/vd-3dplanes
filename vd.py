from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
from functools import singledispatch
from typing import Tuple
import numpy as np
import intersection
import project
import z_dist
import primitives

# compute the 2d vertical decomposition on a plane p.
# intersections is the list of all intersection segments of all pairs of planes.
def vd2d(p : Plane, intersections, planes):
    # compute all cells that p is a floor to
    incidents = []
    above = []
    p_segs = []
    for i in intersections:
        if z_dist.is_incident(i, p):
            incidents.append(i)
            p_segs.append(i)
    
        if z_dist.is_above(i, p, 'z'):
            if z_dist.is_directly_above(i, p, planes, 'z'):
                above.append(i)
                p_segs.append(project.project(i, p, 'z'))

    p_points = []
    for i in p_segs:
        p_points.extend(primitives.endpoints(i))

    for i in range(len(p_segs)):
        for j in range(i+1, len(p_segs)):
            p_points.extend(intersection.intersect(p_segs[i], p_segs[j]))


    for p in p_points:
        # find the segment that is y-above p
        seg = z_dist.find_directly_above(p, p_segs, 'y')
        if not seg == None:
            # break seg at p.y
            p_segs.remove(seg)
            new_segs = primitives.break_element(seg, p.y, 'y')
            p_segs.extend(new_segs)

        # find the segment that is y-above p
        seg = z_dist.find_directly_below(p, p_segs, 'y')
        if not seg == None:
            # break seg at p.y
            p_segs.remove(seg)
            new_segs = primitives.break_element(seg, p.y, 'y')
            p_segs.extend(new_segs)


    cells = []
    for s in p_segs:
        if isinstance(s, Segment3D):
            x_floor = min(s.p1.x, s.p2.x)
            x_ceil = max(s.p1.x, s.p2.x)
            mid_point = (s.p1 + s.p2) / 2
        elif isinstance(s, Ray3D):
            mid_point = s.p1 + s.direction
            if s.direction.x > 0:
                x_floor = s.p1.x
                x_ceil = None
            else:
                x_floor = None 
                x_ceil = s.p1.x

        y_floor = Line3D(Point3D(s.p1.x, s.p1.y, 0), Point3D(s.p2.x, s.p2.y, 0))
        y_ceil_seg = z_dist.find_directly_above(mid_point, p_segs, 'y')
        y_ceil = Line3D(Point3D(y_ceil_seg.p1.x, y_ceil_seg.p1.y, 0), Point3D(y_ceil_seg.p2.x, y_ceil_seg.p2.y, 0))

        center_point = (s.p1 + s.p2 + y_ceil_seg.p1 + y_ceil_seg.p2) / 4
        z_floor = p
        z_ceil = z_dist.find_directly_above(center_point, planes, 'z')

        cells.append( [x_floor, x_ceil, y_floor, y_ceil, z_floor, z_ceil])

        # check if s is on the lower envelope induced on p
        y_abyss = z_dist.find_directly_below(mid_point, p_segs, 'y')
        # check if p is on the lower envelope
        z_abyss = z_dist.find_directly_below(center_point, planes, 'z')

        # if p is on the lower envelope
        if z_abyss == None:
            cells.append( [x_floor, x_ceil, y_floor, y_ceil, False, p])

        # if s is on the lower envelope induced on p, add another cell
        if y_abyss == None:
            center_point = mid_point + (0, -1, 0)
            center_point = project.project(center_point, p, 'z')
            z_ceil = z_dist.find_directly_above(center_point, planes, 'z')
            cells.append( [x_floor, x_ceil, False, y_floor, p, z_ceil])

            if z_abyss == None:
                cells.append( [x_floor, x_ceil, False, y_floor, False, p])

    return cells

def vd(planes):
    intersection_lines = []
    n = len(planes)
    for i in range(n):
        for j in range(i+1, n):
            try:
                line = intersection.intersect(planes[i], planes[j])
                intersection_lines.append(line)
            except:
                continue


    # all intersection of 2 planes broken to segments and rays
    intersections = []
    for l in intersection_lines:
        # Get direction vector of line l
        direction = l.direction
        # Scale direction vector to x=1
        direction = direction / direction[0]

        intersection_points = [intersection.intersect(l,p) for p in planes]
        intersection_points.sort(key=lambda p: p.x)

        for i in range(len(intersection_points) - 1):
            # Calculate point on line l at x = intersection_points[i].x
            x = intersection_points[i].x
            y = l.p1.y + direction[1] * (x - l.p1.x) 
            z = l.p1.z + direction[2] * (x - l.p1.x)
            lx1 = Point3D(x, y, z)

            x = intersection_points[i+1].x
            y = l.p1.y + direction[1] * (x - l.p1.x) 
            z = l.p1.z + direction[2] * (x - l.p1.x)
            lx2 = Point3D(x, y, z)

            intersections.append(Segment3D(lx1, lx2))

        intersections.append(Ray3D(intersection_points[0], direction))
        intersections.append(Ray3D(intersection_points[len(intersection_points) - 1], -direction))

    
    cells = []
    for p in planes:
        cells2d = vd2d(...)
        for each c in cells2d:
            add a cell that is c is its z-floor - find its z-ceil

            if c is on the lower envelope:
                add a cell that is c is its z-ceil - with no z-floor

