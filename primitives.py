from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols

def break_segment3D(seg: Segment3D, v, axis: str):
    if axis != 'x':
        raise ValueError("break_element is implemented only when breaking along x-axis")

    if v < min(seg.p1.x, seg.p2.x) or v > max(seg.p1.x, seg.p2.x):
        raise ValueError("break_element v is not in x-range of segment")

    # Find point p on segment with x = v
    direction = seg.p2 - seg.p1
    t = (v - seg.p1.x) / direction.x
    p = seg.p1 + t * direction

    if p == seg.p1 or p == seg.p2:
        return [seg]

    # Break segment into two segments at point p
    seg1 = Segment3D(seg.p1, p)
    seg2 = Segment3D(p, seg.p2)
    return [seg1, seg2]

def break_ray3D(ray: Ray3D, v, axis: str):
    if axis != 'x':
        raise ValueError("break_element is implemented only when breaking along x-axis")

    # For a ray, we only need to check if v is in the correct direction
    if ray.direction.x > 0 and v < ray.p1.x:
        raise ValueError("break_element v is below ray starting point")
    if ray.direction.x < 0 and v > ray.p1.x:
        raise ValueError("break_element v is above ray starting point")

    # Find point p on ray with y = v
    t = (v - ray.p1.x) / ray.direction.x
    p = ray.p1 + t * ray.direction

    if p == ray.p1:
        return [ray]
 
    # Break ray into segment and new ray at point p
    seg = Segment3D(ray.p1, p)
    new_ray = Ray3D(p, p + ray.direction)
    return [seg, new_ray]

def break_line3D(line: Line3D, v, axis: str):
    """Break a line into a ray and a ray at point v along the axis."""
    if axis != 'x':
        raise ValueError("break_element is implemented only when breaking along x-axis")

    # Get direction vector
    direction = line.direction
    
    # If line is vertical (parallel to x-axis), we can't break it at x=v
    if direction.x == 0:
        raise ValueError("Cannot break vertical line along x-axis")

    # Find point p on line with x = v
    t = (v - line.p1.x) / direction.x
    p = line.p1 + t * direction

    # Break line into two rays at point p
    ray1 = Ray3D(p, p - direction)  # Ray going in negative direction
    ray2 = Ray3D(p, p + direction)   # Ray going in positive direction

    return [ray1, ray2]

def break_element(a, v, axis: str):
    """Main break_element function that dispatches to specific implementations based on types."""
    type_a = type(a).__name__

    if type_a == 'Segment3D':
        return break_segment3D(a, v, axis)
    elif type_a == 'Ray3D':
        return break_ray3D(a, v, axis)
    elif type_a == 'Line3D':
        return break_line3D(a, v, axis)
    else:
        raise NotImplementedError(f"break_element not implemented for {type_a}")

def endpoints_segment3D(seg: Segment3D):
    return [seg.p1, seg.p2]

def endpoints_ray3D(ray: Ray3D):
    return [ray.p1]

def endpoints(a):
    """Main endpoints function that dispatches to specific implementations based on types."""
    type_a = type(a).__name__

    if type_a == 'Segment3D':
        return endpoints_segment3D(a)
    elif type_a == 'Ray3D':
        return endpoints_ray3D(a)
    else:
        raise NotImplementedError(f"endpoints not implemented for {type_a}")

def mid_point_segment(seg: Segment3D):
    """Returns the midpoint of a segment."""
    return (seg.p1 + seg.p2) / 2

def mid_point_ray(ray: Ray3D):
    """Returns a point one unit along the ray from its start point."""
    return ray.p1 + ray.direction

def mid_point(a):
    """Main mid_point function that dispatches to specific implementations based on types."""
    type_a = type(a).__name__

    if type_a == 'Segment3D':
        return mid_point_segment(a)
    elif type_a == 'Ray3D':
        return mid_point_ray(a)
    else:
        raise NotImplementedError(f"mid_point not implemented for {type_a}")
