from sympy import Point3D, Plane, Line3D, Ray3D, Segment3D, oo, solve, symbols
from functools import singledispatch
from typing import Tuple
from scipy.spatial import Delaunay
import numpy as np
import intersection
import project
import z_dist
import primitives
import cells

import random
import matplotlib.pyplot as plt

# find the center point in a 2d cell
def find_center_point(c):
    if c[0] == None:
        x = c[1] - 1
    elif c[1] == None:
        x = c[0] + 1
    else:
        x = (c[0] + c[1]) / 2

    if c[2] == None:
        return project.project(Point3D(x, 0, 0), c[3], 'y') + Point3D(0, -1, 0)
    elif c[3] == None:
        return project.project(Point3D(x, 0, 0), c[2], 'y') + Point3D(0, 1, 0)
    else:
        y1 = project.project(Point3D(x, 0, 0), c[2], 'y')
        y2 = project.project(Point3D(x, 0, 0), c[3], 'y')
        return Point3D(x, (y1.y + y2.y) / 2, 0)

def find_x_values_of_cell(s):
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

    return x_floor, mid_point, x_ceil

def break_segment_at_points(big_s, points):
    # break the segment into smaller segments at points projected onto it 
    segs = [big_s]
    for pa in points:
        if not z_dist.incident(pa, big_s):
            raise ValueError("pa is not incident to big_s")
        # Find segment s in segs_above that contains point pa
        s = None
        for seg in segs:
            # if not z_dist.incident(seg.p1, big_s):
            #     print("seg is not a part of big_s. seg =", seg, "big_s =", big_s)
            #     raise ValueError("seg.p1 is not incident to big_s")
            # if not z_dist.incident(seg.p2, big_s):
            #     raise ValueError("seg.p2 is not incident to big_s")

            if z_dist.incident(pa, seg):
                s = seg
                break

        if s == None:
            raise IndexError("pa is not incident to any s")
        
        # Remove the found segment from segs_above
        segs.remove(s)
        new_segs = primitives.break_element(s, pa.x, 'x')
        segs.extend(new_segs)
    return segs

# compute the 2d vertical decomposition on a plane p.
# p_segs is a list of segments and rays.
# it is assumed that the segments and rays are non intersecting.
# If there were any intersections in the input it is assumed that lines and segments were broken down into smaller segments and rays that do not intersect
def vd2d(p : Plane, p_segs):
    global points_above
    global points_below

    p_points = []
    for i in p_segs:
        p_points.extend(primitives.endpoints(i))

    # points_above[s] are the projections of points of events above s
    points_above = {}
    # points_below[s] are the projections of points of events below s
    points_below = {}

    for s in p_segs:
        points_above[s] = []
        points_below[s] = []

    for i in range(len(p_segs)):
        for j in range(i+1, len(p_segs)):
            if intersection.parallel(p_segs[i], p_segs[j]):
                continue
            int = intersection.intersect(p_segs[i], p_segs[j])
            p_points.extend(int)
            points_above[p_segs[i]].extend(int)
            points_above[p_segs[j]].extend(int)
            points_below[p_segs[i]].extend(int)
            points_below[p_segs[j]].extend(int)

    for p in p_points:
        # find the segment that is y-above p
        s = z_dist.find_directly_above(p, p_segs, 'y')
        if not s == None:
            projp = project.project(p, s, 'y')
            if projp == None:
                raise ValueError("no projection found")
            if not z_dist.incident(projp, s):
                raise ValueError("projection is not incident to segment")
            points_below[s].append(projp)

        # find the segment that is y-above p
        s = z_dist.find_directly_below(p, p_segs, 'y')
        if not s == None:
            projp = project.project(p, s, 'y')
            if projp == None:
                raise ValueError("no projection found")
            if not z_dist.incident(projp, s):
                raise ValueError("projection is not incident to segment")
            points_above[s].append(projp)


    cells2d = []
    for big_s in p_segs:
        segs_above = break_segment_at_points(big_s, points_above[big_s])

        for s in segs_above:
            x_floor, mid_point, x_ceil = find_x_values_of_cell(s)

            y_floor = Line3D(Point3D(s.p1.x, s.p1.y, 0), Point3D(s.p2.x, s.p2.y, 0))
            y_ceil_seg = z_dist.find_directly_above(mid_point, p_segs, 'y')
            y_ceil = None
            if y_ceil_seg != None:
                y_ceil = Line3D(Point3D(y_ceil_seg.p1.x, y_ceil_seg.p1.y, 0), Point3D(y_ceil_seg.p2.x, y_ceil_seg.p2.y, 0))

            cells2d.append( [x_floor, x_ceil, y_floor, y_ceil] )


        # break the segment into smaller segments at points projected onto it 
        segs_below = break_segment_at_points(big_s, points_below[big_s])
        for s in segs_below:
            x_floor, mid_point, x_ceil = find_x_values_of_cell(s)

            # check if s is on the lower envelope induced on p
            y_abyss = z_dist.find_directly_below(mid_point, p_segs, 'y')

            # if s is on the lower envelope induced on p, add another cell
            if y_abyss == None:
                y_ceil = Line3D(Point3D(s.p1.x, s.p1.y, 0), Point3D(s.p2.x, s.p2.y, 0))
                cells2d.append( [x_floor, x_ceil, None, y_ceil] )

    return cells2d

def vd(planes):
    # the intersection lines of the planes when broken into segments by points projected from above and below
    # the segments on the upper face of the plane. These segments are the intersection lines on this plane as well as  intersection segments of other planes projected onto it.
    intersect_segs_above = {}
    intersect_segs_below = {}
    segs_above = {}
    segs_below = {}
    intersection_lines = {}
    for p in planes:
        intersection_lines[p] = []
        intersect_segs_above[p] = []
        intersect_segs_below[p] = []
        segs_above[p] = []
        segs_below[p] = []

    n = len(planes)
    for i in range(n):
        for j in range(i+1, n):
            line = intersection.intersect(planes[i], planes[j])
            if type(line[0]) != Line3D:
                raise ValueError("line is not a Line3D")

            intersection_lines[planes[i]].append(line[0])
            intersection_lines[planes[j]].append(line[0])

    # compute intersect_segs_above for each plane
    for i in range(len(planes)):
        # find the points where intersect_segs_above should break
        for s_focus in intersection_lines[planes[i]]:
            s_focus_proj = project.project(s_focus, "xy", 'z')
            break_points_above = [] 
            break_points_below = [] 

            for j in range(len(planes)):
                if i == j:
                    continue
                int_point = intersection.intersect(s_focus, planes[j])
                if type(int_point[0]) != Point3D:
                    continue
                int_point = int_point[0]
                # add intersection points of 3 planes
                break_points_above.append(int_point)
                break_points_below.append(int_point) 

                # find 2 intersection lines that intersect on their projection onto the xy plane
                for s_peer in intersection_lines[planes[j]]:
                    # if they actually intersect then their intersection point was already added as an intersection point of 3 planes
                    if intersection.intersect(s_peer, s_focus) != []:
                        continue
                    s_peer_proj = project.project(s_peer, "xy", 'z')
                    int_point_proj = intersection.intersect(s_focus_proj, s_peer_proj)
                    if int_point_proj == []:
                        continue
                    int_point_proj = int_point_proj[0]
                    int_point = project.project(int_point_proj, s_focus, 'z')
                    s_focus_height = z_dist.height(int_point_proj, s_focus, 'z')
                    s_peer_height = z_dist.height(int_point_proj, s_peer, 'z')

                    # check that no other plane is between s_focus and s_peer
                    vertically_visible = True
                    k = 0
                    while k < len(planes) and vertically_visible:
                        if k == i or k == j:
                            k += 1
                            continue
                        k_height = z_dist.height(int_point_proj, planes[k], 'z')
                        if k_height > min(s_focus_height, s_peer_height) and k_height < max(s_focus_height, s_peer_height):
                            vertically_visible = False
                            break
                        k += 1

                    if not vertically_visible:
                        continue

                    if s_focus_height < s_peer_height:
                        break_points_above.append(int_point)
                    else:
                        break_points_below.append(int_point)

            intersect_segs_above[planes[i]].extend(break_segment_at_points(s_focus, break_points_above))
            intersect_segs_below[planes[i]].extend(break_segment_at_points(s_focus, break_points_below))


    # project each segment in intersect_segs_above to be a seg_below of some other plane
    for p in planes:
        for s in intersect_segs_above[p]:
            other_p = z_dist.find_directly_above(s, planes, 'z')
            if not other_p == None:
                s_proj = project.project(s, other_p, 'z')
                segs_below[other_p].append(s_proj)

        for s in intersect_segs_below[p]:
            other_p = z_dist.find_directly_below(s, planes, 'z')
            if not other_p == None:
                s_proj = project.project(s, other_p, 'z')
                segs_above[other_p].append(s_proj)


    for p in planes:
        segs_above[p].extend(intersect_segs_above[p])
        segs_below[p].extend(intersect_segs_below[p])


    cells_list = []
    for p in planes:
        # compute the 2d vertical decomposition on the upper face of the plane
        proj_segs_above = []
        for s in segs_above[p]:
            proj_s = project.project(s, "xy", 'z')
            proj_segs_above.append(proj_s)
        cells2d = vd2d(p, proj_segs_above)

        for c in cells2d:
            center_point = find_center_point(c)
            center_point = project.project(center_point, p, 'z')

            plane_above = z_dist.find_directly_above(center_point, planes, 'z')
            cells_list.append( [c[0], c[1], c[2], c[3], p, plane_above ] )

        # compute the cells that are  below the arrangement.
        # To do that, compute the 2d vertical decomposition on the lower face of the plane and find cells that are unbounded from below
        proj_segs_below = []
        for s in segs_below[p]:
            proj_s = project.project(s, "xy", 'z')
            proj_segs_below.append(proj_s)
        cells2d = vd2d(p, proj_segs_below)

        for c in cells2d:
            center_point = find_center_point(c)
            center_point = project.project(center_point, p, 'z')

            plane_below = z_dist.find_directly_below(center_point, planes, 'z')
            if plane_below == None:
                cells_list.append( [c[0], c[1], c[2], c[3], None, p ] )

    return cells_list

def test_vd2d():
    global points_above
    global points_below

    # Create XY plane at z=0
    xy_plane = Plane(Point3D(0,0,0), Point3D(1,0,0), Point3D(0,1,0))
    
    # Initialize empty list for segments
    segs = []

    # Generate 10 random segments on XY plane
    while len(segs) < 3:
        # Generate random points on XY plane
        p1 = Point3D(random.uniform(-10,10), random.uniform(-10,10), 0)
        p2 = Point3D(random.uniform(-10,10), random.uniform(-10,10), 0)
        seg = Segment3D(p1, p2)
        segs.append(seg)
            
    # Generate 10 random rays on XY plane
    for i in range(3):
        # Generate random point and direction on XY plane
        p = Point3D(random.uniform(-10,10), random.uniform(-10,10), 0)
        dir = Point3D(random.uniform(-1,1), random.uniform(-1,1), 0)
        ray = Ray3D(p, dir)
        segs.append(ray)

    # Calculate Voronoi diagram cells
    cells_list = []
    points_above = {}
    points_below = {}
    cells_list = vd2d(xy_plane, segs)

    # Draw visualization
    plt.figure(figsize=(10,10))
    
    # Draw segments and rays
    for s in segs:
        if isinstance(s, Segment3D):
            plt.plot([float(s.p1.x), float(s.p2.x)], 
                    [float(s.p1.y), float(s.p2.y)], 'b-')
        else: # Ray3D
            # Convert direction components to float for plotting
            dx = float(s.direction.x)
            dy = float(s.direction.y)
            # Normalize direction vector to fixed length for visualization
            length = 2.0  # Adjust this value to change arrow length
            norm = (dx * dx + dy * dy) ** 0.5
            dx = dx * length / norm
            dy = dy * length / norm
            plt.arrow(float(s.p1.x), float(s.p1.y), dx, dy,
                     head_width=0.3, head_length=0.5, fc='b', ec='b')
    
    # Draw cell boundaries
    for cell in cells_list:
        x_floor, x_ceil = cell[0], cell[1]
        y_floor, y_ceil = cell[2], cell[3]
        
        # Draw vertical lines for x boundaries
        if x_floor is not False and x_floor is not None:
            if y_floor is None:
                y_bottom = -10
            else:
                y_bottom = project.project(Point3D(x_floor, 0, 0), y_floor, 'y').y
            if y_ceil is None:
                y_top = 10
            else:
                y_top = project.project(Point3D(x_floor, 0, 0), y_ceil, 'y').y

            plt.vlines(x=float(x_floor + 0.1), ymin=y_bottom, ymax=y_top, color='r', linestyle='--', alpha=0.5)
            
        if x_ceil is not False and x_ceil is not None:
            if y_floor is None:
                y_bottom = -10
            else:
                y_bottom = project.project(Point3D(x_ceil, 0, 0), y_floor, 'y').y
            if y_ceil is None:
                y_top = 10
            else:
                y_top = project.project(Point3D(x_ceil, 0, 0), y_ceil, 'y').y

            plt.vlines(x=float(x_ceil - 0.1), ymin=y_bottom, ymax=y_top, color='r', linestyle='--', alpha=0.5)

        # Draw horizontal lines for y boundaries
        if y_floor is not False and y_floor is not None:
            # For Line3D, get y coordinate at x_floor and x_ceil
            if isinstance(y_floor, Line3D):
                plt.plot([float(y_floor.p1.x), float(y_floor.p2.x)], 
                        [float(y_floor.p1.y + 0.1), float(y_floor.p2.y + 0.1)], 'r--', alpha=0.5)
                
        if y_ceil is not False and y_ceil is not None:
            # For Line3D, get y coordinate at x_floor and x_ceil
            if isinstance(y_ceil, Line3D):
                plt.plot([float(y_ceil.p1.x), float(y_ceil.p2.x)], 
                        [float(y_ceil.p1.y - 0.1), float(y_ceil.p2.y - 0.1)], 'r--', alpha=0.5)

    # Plot points_above with small vertical offset
    for s in segs:
        if s in points_above:
            for p in points_above[s]:
                offset_point = p + Point3D(0, 0.1, 0)
                plt.plot(float(offset_point.x), float(offset_point.y), 'go', markersize=5)  # green dots for points_above

        if s in points_below:
            for p in points_below[s]:
                offset_point = p - Point3D(0, 0.1, 0)
                plt.plot(float(offset_point.x), float(offset_point.y), 'ro', markersize=5)  # green dots for points_above

    plt.grid(True)
    plt.axis('equal')
    plt.show()



def visualize_3d_cells(planes, cells_list):
    """
    Visualize 3D cells using matplotlib's Poly3DCollection.
    Each cell is displayed with a different color.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import numpy as np
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    print("plotting planes ...")
    for p in planes:
        X = np.arange(-10, 10, 0.5)
        Y = np.arange(-10, 10, 0.5)
        X, Y = np.meshgrid(X, Y)
        Z = np.zeros_like(X)
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                Z[i, j] = project.project(Point3D(X[i, j], Y[i, j], 0), p, 'z').z
        ax.plot_surface(X, Y, Z, color='grey', alpha=0.5)
    print("done")

    print("plotting intersections ...")
    for p in planes:
        for q in planes:
            if p == q:
                continue
            int = intersection.intersect(p, q)[0]
            points = []
            for wall in [
                        Plane(Point3D(10,0,0), Point3D(10,10,0), Point3D(10,0,10)),
                        Plane(Point3D(-10,0,0), Point3D(-10,10,0), Point3D(-10,0,10)),
                        Plane(Point3D(0,10,0), Point3D(10,10,0), Point3D(0,10,10)),
                        Plane(Point3D(0,-10,0), Point3D(10,-10,0), Point3D(0,-10,10))]:
                    point = intersection.intersect(int, wall)[0]
                    if point.x <= 10 and point.x >= -10 and point.y <= 10 and point.y >= -10 and point.z <= 10 and point.z >= -10:
                        points.append(point)
            if points != []:
                poly = Poly3DCollection([points], alpha=0.5, color='cyan', edgecolor='k')
                ax.add_collection3d(poly)
    print("done")
    
    # Generate distinct colors for cells
    n_cells = len(cells_list)
    colors = plt.cm.rainbow(np.linspace(0, 1, n_cells))
    
    for cell, color in zip(cells_list, colors):
        print("plotting cell ...")
        wall_x_floor = []
        wall_x_ceil = []
        wall_y_floor = []
        wall_y_ceil = []
        # wall_x_floor = cells.get_cell_x_floor_surface(cell)
        # wall_x_ceil = cells.get_cell_x_ceil_surface(cell)
        wall_y_floor = cells.get_cell_y_floor_surface(cell)
        wall_y_ceil = cells.get_cell_y_ceil_surface(cell)
        # wall_z_floor = cells.get_cell_z_floor_surface(cell)
        # wall_z_ceil = cells.get_cell_z_ceil_surface(cell)
        for wall in [wall_x_floor, wall_x_ceil, wall_y_floor, wall_y_ceil]:
            if len(wall) == 0:
                continue
            print("drawing wall with", len(wall), "points")
            poly = Poly3DCollection([wall], alpha=0.5, color='cyan', edgecolor='k')
            ax.add_collection3d(poly)
        print("done")
    
    # Set axis labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # Auto-scale the axes
    ax.auto_scale_xyz([-10, 10], [-10, 10], [-10, 10])
    
    plt.show()

def test_vd2d_slanted():
    global points_above
    global points_below

    # Create XY plane at z=0
    slanted_plane = Plane(Point3D(0,0,1), Point3D(1,0,2), Point3D(0,1,3))
    
    # Initialize empty list for segments
    segs = []

    # Generate 10 random segments on XY plane
    while len(segs) < 3:
        # Generate random points on XY plane
        p1 = Point3D(random.uniform(-10,10), random.uniform(-10,10), 0)
        p1 = project.project(p1, slanted_plane, 'z')
        p2 = Point3D(random.uniform(-10,10), random.uniform(-10,10), 0)
        p2 = project.project(p2, slanted_plane, 'z')
        seg = Segment3D(p1, p2)
        segs.append(seg)
            
    # Generate 10 random rays on XY plane
    for i in range(3):
        # Generate random point and direction on XY plane
        p = Point3D(random.uniform(-10,10), random.uniform(-10,10), 0)
        p = project.project(p, slanted_plane, 'z')
        dir = Point3D(random.uniform(-1,1), random.uniform(-1,1), 0)
        dir = project.project(dir, slanted_plane, 'z')
        ray = Ray3D(p, dir)
        segs.append(ray)

    # Calculate Voronoi diagram cells
    cells_list = []
    points_above = {}
    points_below = {}
    cells_list = vd2d(slanted_plane, segs)

    # Draw visualization
    plt.figure(figsize=(10,10))
    
    # Draw segments and rays
    for s in segs:
        if isinstance(s, Segment3D):
            plt.plot([float(s.p1.x), float(s.p2.x)], 
                    [float(s.p1.y), float(s.p2.y)], 'b-')
        else: # Ray3D
            # Convert direction components to float for plotting
            dx = float(s.direction.x)
            dy = float(s.direction.y)
            # Normalize direction vector to fixed length for visualization
            length = 2.0  # Adjust this value to change arrow length
            norm = (dx * dx + dy * dy) ** 0.5
            dx = dx * length / norm
            dy = dy * length / norm
            plt.arrow(float(s.p1.x), float(s.p1.y), dx, dy,
                     head_width=0.3, head_length=0.5, fc='b', ec='b')
    
    # Draw cell boundaries
    for cell in cells_list:
        x_floor, x_ceil = cell[0], cell[1]
        y_floor, y_ceil = cell[2], cell[3]
        
        # Draw vertical lines for x boundaries
        if x_floor is not False and x_floor is not None:
            if y_floor is None:
                y_bottom = -10
            else:
                y_bottom = project.project(Point3D(x_floor, 0, 0), y_floor, 'y').y
            if y_ceil is None:
                y_top = 10
            else:
                y_top = project.project(Point3D(x_floor, 0, 0), y_ceil, 'y').y

            plt.vlines(x=float(x_floor + 0.1), ymin=y_bottom, ymax=y_top, color='r', linestyle='--', alpha=0.5)
            
        if x_ceil is not False and x_ceil is not None:
            if y_floor is None:
                y_bottom = -10
            else:
                y_bottom = project.project(Point3D(x_ceil, 0, 0), y_floor, 'y').y
            if y_ceil is None:
                y_top = 10
            else:
                y_top = project.project(Point3D(x_ceil, 0, 0), y_ceil, 'y').y

            plt.vlines(x=float(x_ceil - 0.1), ymin=y_bottom, ymax=y_top, color='r', linestyle='--', alpha=0.5)

        # Draw horizontal lines for y boundaries
        if y_floor is not False and y_floor is not None:
            # For Line3D, get y coordinate at x_floor and x_ceil
            if isinstance(y_floor, Line3D):
                plt.plot([float(y_floor.p1.x), float(y_floor.p2.x)], 
                        [float(y_floor.p1.y + 0.1), float(y_floor.p2.y + 0.1)], 'r--', alpha=0.5)
                
        if y_ceil is not False and y_ceil is not None:
            # For Line3D, get y coordinate at x_floor and x_ceil
            if isinstance(y_ceil, Line3D):
                plt.plot([float(y_ceil.p1.x), float(y_ceil.p2.x)], 
                        [float(y_ceil.p1.y - 0.1), float(y_ceil.p2.y - 0.1)], 'r--', alpha=0.5)

    # Plot points_above with small vertical offset
    for s in segs:
        if s in points_above:
            for p in points_above[s]:
                offset_point = p + Point3D(0, 0.1, 0)
                plt.plot(float(offset_point.x), float(offset_point.y), 'go', markersize=5)  # green dots for points_above

        if s in points_below:
            for p in points_below[s]:
                offset_point = p - Point3D(0, 0.1, 0)
                plt.plot(float(offset_point.x), float(offset_point.y), 'ro', markersize=5)  # green dots for points_above

    plt.grid(True)
    plt.axis('equal')
    plt.show()



def test_vd():
    global points_above
    global points_below

    # Create XY plane at z=0
    planes = []
    while len(planes) < 4:
        p1 = Point3D(random.uniform(-10,10), random.uniform(-10,10), random.uniform(-10,10))
        p2 = Point3D(random.uniform(-10,10), random.uniform(-10,10), random.uniform(-10,10))
        p3 = Point3D(random.uniform(-10,10), random.uniform(-10,10), random.uniform(-10,10))
        plane = Plane(p1, p2, p3)
        planes.append(plane)

    cells_list = []
    cells_list = vd(planes)

    print("There are", len(cells_list), "cells")

    # checks cells are disjoint
    # for i in range(len(cells_list)):
    #     endpoints = cells.find_cell_vertices(cells_list[i])
    #     for j in range(len(cells_list)):
    #         if i == j:
    #             continue
    #         for p in endpoints:
    #             if cells.is_point_in_cell(p, cells_list[j]):
    #                 raise ValueError("cell", i, "contains point", p, "of cell", j)

    visualize_3d_cells(planes, cells_list)

            

if __name__ == "__main__":
    random.seed(30)
    # test_vd2d()
    # test_vd2d_slanted()
    test_vd()