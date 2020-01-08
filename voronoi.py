"""
Voronoi diagram from Delaunay triangulation, Jonas Nockert 2020

1.  Deluanay triangulate points
1.1 Triangulate using sweep-line algorithm (some triangles will not be legal
    delaunay).
1.2 Legalize triangles using Lawson legalization.
2. Create Voronoi diagram from deluanay triangulation since it is its dual.

Mainly based on:

Yonghe et al. (2013) A Simple Sweep-line Delaunay Triangulation Algorithm
http://www.academicpub.org/jao/paperInfo.aspx?PaperID=15630

Reference material:

Borut Zalik (2005) An efficient sweep-line Delaunay triangulation algorithm

Lawson (1977) Software for C1 Surface Interpolation.
https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770025881.pdf

https://www.cse.iitb.ac.in/~cs749/spr2017/handouts/floater_mesh_simplify_optimize.pdf

https://www.redblobgames.com/x/1721-voronoi-alternative/

https://leatherbee.org/index.php/2018/10/06/terrain-generation-3-voronoi-diagrams/
"""
from collections import deque
from random import gammavariate, gauss, randint, random, shuffle

from pygame.math import Vector2

from bidirectional_list import BidirectionalList as BiList, BidirectionalNode as BiNode


class VoronoiDiagram:
    def __init__(self, points):
        self.points = points
        self.triangles, self.triangle_edges = self.delaunay_triangulation()
        self.edges, self.cells = self.create_voronoi_cells()

    def create_maze(self):
        """Use a variation of the Growing Tree algorithm to form a maze

        The Recursive Backtracker algorithm consists of the following steps

        1. Adds a starting cell to an empty stack, marking the cell as visited,
        2. pops a cell of the top of the stack,
        3. finds a non-visited neighbor to that cell (or step 7 in case none exists),
        4. removes the edge between this neighbor and cell
        5. re-adds the cell to the top of the stack in case there are more neighbors,
        6. adds the neighbor to the top of the stack, marking it as visited,
        7. goes to step 2 until the queue is empty

        The Growing Tree algorithm is a variation of the above in regards to
        how the next cell is picked in step 2. If the growing tree algorithm is
        set to always choose the most recently added cell, it is equivalent to
        Recursive Backtracking. If cells are chosen at random, it is equivalent
        to Prim's algorithm (a minimum spanning tree of a graph can be viewed
        as a maze). Another option is to pick the oldest cells in the queue,
        which tends to create long straight corridors and not very good mazes.

        It seems that often a ratio is picked between two or more of the above cases,
        e.g. randomly picking most recent or by random in a 50/50 split.

        Here, for variation, choices are sampled according to the gamma distribution
        which can be seen as the picks are done with some regularity in terms
        of rate of occurence. I've chosen to pick the oldest cells with a mean
        of every sixth time with most recent cells picked in between.

        The idea (hope) is that this will lead to longer winding dead-end branches
        off the "main" path but with a guarantee that there's always some winding
        segments between branches as the gamma distribution is positive valued.
        The branches should not compete too much for space in order for both
        branches and main path to achieve a degree of windingness.

        The benefit is longer dead-end corridors but the downside is that the
        paths generated are shorter than what Recursive Backtracker algorithm
        would produce. That is, this is a trade-off between DFS and BFS.

        Sources:
        https://en.wikipedia.org/wiki/Maze_generation_algorithm#Recursive_backtracker
        https://en.wikipedia.org/wiki/Gamma_distribution
        http://weblog.jamisbuck.org/2011/1/27/maze-generation-growing-tree-algorithm
        """
        assert self.triangles
        assert self.triangle_edges
        assert self.cells
        assert self.edges

        for cell in self.cells:
            cell.maze_edges = [True] * len(cell.edges)

        # Use a stack instead of recursion.
        q = deque()
        # Keep track of visited cells.
        visited = set()

        # Add a starting cell. In this case cell 0 is on the edge of the graph
        # but not in any predefined place. It depends on the initial point placements
        # and the subsequent triangulation.
        start_ix = 0
        # Color initial cell.
        self.cells[start_ix].set_color_rgb(0, 125, 0)
        visited.add(start_ix)
        q.append(start_ix)

        alpha = 3
        beta = 2
        delay = gammavariate(alpha, beta)
        while q:
            # Pick a cell from the left or right.
            if delay <= 0:
                # print("left pop")
                cell_ix = q.popleft()
                delay = gammavariate(alpha, beta) + 1
            else:
                # print("right pop")
                cell_ix = q.pop()
            delay -= 1
            current_cell = self.cells[cell_ix]

            # Look for an unvisited neighbor to the cell.
            neighbor_indices = list(range(len(current_cell.neighbors)))
            # The graph is generally not at all regular but the neighbors
            # are defined clockwise so perhaps checking neighbors in random
            # order helps avoid symmetries a little bit.
            shuffle(neighbor_indices)
            for i in neighbor_indices:
                n_ix = current_cell.neighbors[i]
                if n_ix is None or n_ix in visited:
                    continue

                # Found not visited neighbor.
                ce0, ce1 = current_cell.edges[i]
                neighbor_cell = self.cells[n_ix]

                # Remove edge from neighbor to cell.
                edge_indices = list(range(len(neighbor_cell.edges)))
                for j in edge_indices:
                    n_edge = neighbor_cell.edges[j]
                    if not n_edge:
                        continue
                    nce0, nce1 = n_edge
                    if (ce0, ce1) == (nce1, nce0):
                        neighbor_cell.maze_edges[j] = False
                    if (ce1, ce0) == (nce0, nce1):
                        neighbor_cell.maze_edges[j] = False
                # Remove edge from cell to neighbor.
                current_cell.maze_edges[i] = False

                # Since we found an unvisited neighbor, perhaps there are more
                # so add the cell to the top of the stack again.
                q.append(cell_ix)

                # Add neighbor to top of stack. In standard recursive backtracking
                # this cell will be the one to process next so it is essentially
                # a form of depth-first search.
                visited.add(n_ix)
                q.append(n_ix)
                break

        return start_ix

    def create_voronoi_cells(self):
        """Connect edges into cells and connect cells with their neighbors"""
        assert self.triangles
        assert self.triangle_edges
        vedges = self._create_voronoi_edges()

        cell_index = 0
        cells = []
        for edge, edge_data in vedges.items():
            if "used" in edge_data and edge_data["used"]:
                continue
            vedges[edge]["used"] = True

            ti0, ti1 = edge
            t0 = self.triangles[ti0]
            t1 = self.triangles[ti1]
            cell_triangles = [t0, t1]
            cell_edges = [edge]

            # Starting with edge between the circumcenters of t0 and t1, follow
            # edges counterclockwise until reaching (t0, t1) again.
            while True:
                p0 = t0.circumcenter
                p1 = t1.circumcenter
                dest_neighbors = t1.neighbors - set([ti0])
                p2 = None
                ti2 = None
                if dest_neighbors:
                    ti2_test = dest_neighbors.pop()
                    p2_test = self.triangles[ti2_test].circumcenter
                    if not self.point_on_left_side(p0, p1, p2_test):  # >= 0:
                        ti2 = ti2_test
                        # print("Found left neighbor (1)", ti2)
                if dest_neighbors:
                    ti2_test = dest_neighbors.pop()
                    p2_test = self.triangles[ti2_test].circumcenter
                    if not self.point_on_left_side(p0, p1, p2_test):  # >= 0:
                        if ti2:
                            print("Already have a left neighbor!!!")
                        ti2 = ti2_test
                        # print("Found left neighbor (2)", ti2)
                if not ti2:
                    print("No left neighbor")
                    break
                next_edge = (ti1, ti2)
                if "used" in vedges[next_edge] and vedges[next_edge]["used"]:
                    if len(cell_triangles) > 2:
                        for ce0, ce1 in cell_edges:
                            vedges[(ce0, ce1)]["neighbors"].append(cell_index)
                            vedges[(ce1, ce0)]["neighbors"].append(cell_index)
                        vcell = self.VoronoiCell(cell_edges, cell_triangles)
                        cells.append(vcell)
                        cell_index += 1
                        # print("complete cell", cell_triangles)
                    else:
                        print("Cell not complete")
                    break
                vedges[next_edge]["used"] = True
                cell_triangles.append(self.triangles[ti2])
                cell_edges.append(next_edge)
                ti0 = ti1
                ti1 = ti2
                t0 = self.triangles[ti0]
                t1 = self.triangles[ti1]

        # For each shared edge, get the neighbouring cell.
        for i, cell in enumerate(cells):
            neighbors = []
            for ce in cell.edges:
                ns = vedges[ce]["neighbors"]
                if len(ns) == 1:
                    neighbors.append(None)
                    continue
                if len(ns) < 1:
                    raise Exception("Less than one cell neighbor!")
                if len(ns) > 2:
                    raise Exception("More than two cell neighbors!")
                n1, n2 = ns
                if n1 == i:
                    neighbors.append(n2)
                else:
                    neighbors.append(n1)
            cell.set_neighbors(neighbors)

        return vedges, cells

    def _create_voronoi_edges(self):
        """Return edges as ordered pairs of triangle indices

        Note that each edge is added in both directions, e.g. (ti0, ti1) as
        well as (ti1, ti0). This as we want all polygons to be defined in CCW
        order.

        Note also that, technically, the voronoi edges are between the
        circumcenters of the triangle pair.
        """
        assert self.triangles
        assert self.triangle_edges
        vedges = {}
        for i, tri in enumerate(self.triangles):
            if tri.flagged_for_deletion:
                continue

            indices = tri.point_indices
            e0 = self.edge_tuple(indices[0], indices[1])
            e1 = self.edge_tuple(indices[1], indices[2])
            e2 = self.edge_tuple(indices[2], indices[0])
            neighbors = (
                set(self.triangle_edges[e0])
                .union(set(self.triangle_edges[e1]))
                .union(set(self.triangle_edges[e2]))
            )
            neighbors.remove(i)
            tri.set_neighbors(neighbors)
            for n in neighbors:
                nt = self.triangles[n]
                if nt.flagged_for_deletion:
                    continue
                if (n, i) in vedges:
                    vedges[(i, n)] = {"first": False, "neighbors": []}
                else:
                    vedges[(i, n)] = {"first": True, "neighbors": []}
        return vedges

    def delaunay_triangulation(self):
        """Create Delaunay triangulation from a list of points.

        Legalization of triangles (Lawson legalization), loosely based on
        Lawson (1977) "Software for C1 Surface Interpolation"
        (https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770025881.pdf):

                 *                   *
                /|\                 / \
            e1 / | \ e4         e1 /   \ e4
              /  |t2\             / t3  \
             /   |   \           /       \
            *    e    *   ->    *----e'---*
             \   |   /           \       /
              \t1|  /             \ t4  /
            e2 \ | / e3         e2 \   / e3
                \|/                 \ /
                 *                   *

        Legalization of the pair of triangles (t1, t2) sharing edge e results
        in the removal of e and addition of e' with t3 and t4 sharing the latter
        edge (t1 and t2 are removed). Edges e1-e4 have not changed but the triangles
        sharing these edges might now need to be legalized (in a recursive fashion).

        Initially, all edges (pairs of point indices) are put in a FIFO queue together
        with an associated value c. We also initialize a hash associating each edge
        with the same value c.

        Then, until the queue is empty
        1. Pop edge e from the queue.
        2. If the associated value c is not the same as in the hash, the edge has been
        superceded and is skipped.
        3. Retrieve the triangles (t1, t2) sharing the edge.
        4. If (t1, t2) are not legal,
        a. Remove edge e and add edge e'. There is no need to add e' to queue as its
            triangles are legal.
        b. Remove triangles t1 and t2, adding t3 and t4.
        c. Add all edges for t3 and t4 except e' to the queue. Since these edges
            might already be in the queue (but we don't know and it would probably
            be more expensive to find out) we just add them again with a new associated
            value c. This way, earlier duplicates will e skipped and the operations are
            ordered (not absolutely sure they need to be but we definitely want to
            avoid costly legalization checking if we can).
        """
        assert self.points
        # Create initial triangulation.
        triangles, triangle_edges = self._pretriangulate()

        # Add all edges to FIFO queue.
        counter = 1
        q = deque()
        qc = {}  # Keep track of latest in order to skip earlier duplicates.
        for key in triangle_edges.keys():
            q.append((key, counter))
            qc[key] = counter
            counter += 1

        # Go through edges one by one until all pairs of triangles sharing
        # the edge are legal.
        while True:
            try:
                edge, c = q.popleft()
            except IndexError:
                break

            # Has the edge been superceded/made irrelevant?
            if qc[edge] != c:
                continue

            # Edge point indices are not in any specific order.
            try:
                if len(triangle_edges[edge]) < 2:
                    continue
            except KeyError:
                # Just a precaution.
                print(triangle_edges)
                raise
            e0, e1 = edge
            ti0, ti1 = triangle_edges[edge]
            t0 = triangles[ti0]
            t1 = triangles[ti1]

            # Get point in t0 that is not in t1.
            e0_ix0 = t0.point_indices.index(e0)
            if t0.point_indices[(e0_ix0 + 1) % 3] == e1:
                pt0 = t0.point_indices[(e0_ix0 + 2) % 3]
            else:
                pt0 = t0.point_indices[(e0_ix0 + 1) % 3]

            o, r = self.Triangle.circumcircle(
                *[self.points[ix] for ix in t1.point_indices]
            )
            if o.distance_to(self.points[pt0]) <= r:
                # Flip edge.

                # Build ccw quad.
                e0_ix1 = t1.point_indices.index(e0)
                if t1.point_indices[(e0_ix1 + 1) % 3] == e1:
                    # e1 immediately after e0. Insert point between e0 and e1.
                    quad = [
                        t1.point_indices[e0_ix1],  # e0.
                        pt0,
                        t1.point_indices[(e0_ix1 + 1) % 3],  # e1.
                        t1.point_indices[(e0_ix1 + 2) % 3],
                    ]
                else:
                    # e1 must be immediately before e0. Insert point between e1 and e0.
                    quad = [
                        t1.point_indices[e0_ix1],  # e0.
                        t1.point_indices[(e0_ix1 + 1) % 3],
                        t1.point_indices[(e0_ix1 + 2) % 3],  # e1.
                        pt0,
                    ]
                triangles[ti0].flag_for_deletion()
                triangles[ti1].flag_for_deletion()
                qc[edge] = 1e20  # Essentially remove non-lD edge from queue.
                del triangle_edges[edge]

                triangles.append(self.Triangle(self.points, quad[0], quad[1], quad[3]))
                triangles.append(self.Triangle(self.points, quad[1], quad[2], quad[3]))

                new_edge = self.edge_tuple(quad[1], quad[3])
                new_ti0 = len(triangles) - 2
                new_ti1 = new_ti0 + 1
                triangle_edges[new_edge] = (new_ti0, new_ti1)

                q1 = self.edge_tuple(quad[0], quad[1])
                q2 = self.edge_tuple(quad[1], quad[2])
                q3 = self.edge_tuple(quad[2], quad[3])
                q4 = self.edge_tuple(quad[3], quad[0])

                tq1 = set(triangle_edges[q1])
                tq1.difference_update([ti0, ti1])
                tq1.add(new_ti0)
                tq2 = set(triangle_edges[q2])
                tq2.difference_update([ti0, ti1])
                tq2.add(new_ti1)
                tq3 = set(triangle_edges[q3])
                tq3.difference_update([ti0, ti1])
                tq3.add(new_ti1)
                tq4 = set(triangle_edges[q4])
                tq4.difference_update([ti0, ti1])
                tq4.add(new_ti0)

                del triangle_edges[q1]
                del triangle_edges[q2]
                del triangle_edges[q3]
                del triangle_edges[q4]
                triangle_edges[q1] = tuple(tq1)
                triangle_edges[q2] = tuple(tq2)
                triangle_edges[q3] = tuple(tq3)
                triangle_edges[q4] = tuple(tq4)

                # Add all edges in quad to queue.
                q.append((q1, counter))
                q.append((q2, counter + 1))
                q.append((q3, counter + 2))
                q.append((q4, counter + 3))
                qc[q1] = counter
                qc[q2] = counter + 1
                qc[q3] = counter + 2
                qc[q4] = counter + 3
                counter += 4

        return triangles, triangle_edges

    def _pretriangulate(self):
        """ 1.1 Initial triangulation of points

        Return triangles and edges.

        Yonghe et al. (2013) A Simple Sweep-line Delaunay Triangulation Algorithm
        http://www.academicpub.org/jao/paperInfo.aspx?PaperID=15630
        """
        # Sort points by x coordinate to prepare for sweep-line algorithm.
        self.points.sort(key=lambda v: v.x)

        triangles = []
        # An edge is a tuple of indices for a pair of points, always ordered so the
        # lower index comes before the higher. An edge is associated with a list of
        # triangles sharing it (one or two).
        triangle_edges = {}
        # The hull is the current outer border of the triangulation. It is a
        # bidirectional list of edges (a CCW polygon).
        hull = BiList()

        # Add the first triangle.
        if self.point_on_left_side(self.points[0], self.points[1], self.points[2]):
            # print("p2 on left side of p0p1")
            e0 = (0, 1)
            e1 = (1, 2)
            e2 = (2, 0)
            triangles.append(self.Triangle(self.points, 0, 1, 2))
        else:
            # print("p2 on right side of p0p1")
            e0 = (1, 0)
            e1 = (0, 2)
            e2 = (2, 1)
            triangles.append(self.Triangle(self.points, 1, 0, 2))
        hull.append(e0)
        hull.append(e1)
        hull.append(e2)
        triangle_edges[self.edge_tuple(e0)] = [0]
        triangle_edges[self.edge_tuple(e1)] = [0]
        triangle_edges[self.edge_tuple(e2)] = [0]

        pts_index = 3
        for pts_index in range(3, len(self.points)):
            next_pt = self.points[pts_index]

            ix = 0
            while ix < len(hull):
                node = hull[ix]
                i0 = node.value[0]
                i1 = node.value[1]
                p0 = self.points[i0]
                p1 = self.points[i1]
                if not self.point_on_left_side(p0, p1, next_pt):
                    tl = len(triangles)

                    triangles.append(self.Triangle(self.points, i1, i0, pts_index))

                    ab = self.edge_tuple(i0, i1)
                    if ab not in triangle_edges:
                        triangle_edges[ab] = []
                    triangle_edges[ab].append(tl)

                    ab = self.edge_tuple(i1, pts_index)
                    if ab not in triangle_edges:
                        triangle_edges[ab] = []
                    triangle_edges[ab].append(tl)

                    ab = self.edge_tuple(pts_index, i0)
                    if ab not in triangle_edges:
                        triangle_edges[ab] = []
                    triangle_edges[ab].append(tl)

                    # Add in reverse order.
                    succ = node.succ
                    pred = node.pred

                    if not succ:
                        succ = hull[0]
                    if not pred:
                        pred = hull[-1]

                    if succ.value == (i1, pts_index):
                        hull.remove(succ)
                    else:
                        hull.insert_after(node, (pts_index, i1))
                        ix += 1

                    if pred.value == (pts_index, i0):
                        hull.remove(pred)
                        ix -= 1
                    else:
                        hull.insert_after(node, (i0, pts_index))
                        ix += 1

                    hull.remove(node)
                    ix -= 1
                ix += 1

        return triangles, triangle_edges

    #
    # Helpers
    #

    def point_on_left_side(self, p0, p1, p2):
        """If v≤0, p2 is on the left side of p0p1, otherwise it is on the right side

        Yonghe et al. (2013) A Simple Sweep-line Delaunay Triangulation Algorithm
        http://www.academicpub.org/jao/paperInfo.aspx?PaperID=15630
        """
        return ((p2.x - p0.x) * (p1.y - p0.y) - (p1.x - p0.x) * (p2.y - p0.y)) <= 0

    def edge_tuple(self, i0, i1=None):
        if type(i0) is tuple:
            return tuple(sorted(i0))
        return tuple(sorted([i0, i1]))

    #
    # Inner classes
    #

    class VoronoiCell:
        def __init__(self, cell_edges, cell_triangles):
            self.edges = cell_edges
            self.triangles = cell_triangles
            self.polygon = [
                Vector2(t.circumcenter.x, t.circumcenter.y) for t in self.triangles
            ]
            self.color = None

        def set_color_rgb(self, r, g, b):
            self.color = (r, g, b)

        def set_color_random(self):
            self.color = (randint(0, 255), randint(0, 255), randint(0, 255))

        def set_neighbors(self, neighbors):
            """Same order as edges"""
            assert len(neighbors) == len(self.edges)
            self.neighbors = neighbors

        @property
        def center(self):
            """Mass center of voronoi polygon

            Source: https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
            """
            a = 0
            cx = 0
            cy = 0
            if self.polygon[0] != self.polygon[-1]:
                # Some edge cells are extra problematic.
                # print("NOT", self.edges, self.triangles, self.polygon)
                # self.set_color_rgb(255, 0, 0)
                pass
            polylen = len(self.polygon)
            if polylen == 2:
                # Just a precaution.
                print(self.edges, self.triangles, self.polygon)
                raise Exception("Two point polygon")
            for i in range(polylen):
                p = self.polygon[i]
                p_next = self.polygon[(i + 1) % polylen]
                a += p.x * p_next.y - p_next.x * p.y
                shared_factor = p.x * p_next.y - p_next.x * p.y
                cx += (p.x + p_next.x) * shared_factor
                cy += (p.y + p_next.y) * shared_factor
            a /= 2
            cx /= 6 * a
            cy /= 6 * a
            return Vector2(cx, cy)

    class Triangle:
        def __init__(self, points, pi0, pi1, pi2):
            self.point_indices = [pi0, pi1, pi2]
            o, r = self.circumcircle(points[pi0], points[pi1], points[pi2])
            self.circumcenter = o
            self.circumradius = r
            self.flagged_for_deletion = False
            self.neighbors = None

        def flag_for_deletion(self):
            self.flagged_for_deletion = True

        def set_neighbors(self, neighbors):
            self.neighbors = neighbors

        @classmethod
        def circumcircle(self, p0, p1o, p2o):
            """Return the circumcenter and circumradius given three Vector2 points

            Source: https://en.wikipedia.org/wiki/Circumscribed_circle
            """
            # Use p0 as origin.
            p1 = p1o - p0
            p2 = p2o - p0
            p1m = p1.magnitude_squared()
            p2m = p2.magnitude_squared()
            d = 2 * (p1.x * p2.y - p1.y * p2.x)
            ux = (1 / d) * (p2.y * p1m - p1.y * p2m)
            uy = (1 / d) * (p1.x * p2m - p2.x * p1m)
            u = Vector2(ux, uy)
            r = u.magnitude()
            return u + p0, r
