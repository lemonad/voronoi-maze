"""
Voronoi maze generation using recursive backtracking.

Jonas Nockert 2020
"""
from random import gauss, randint, random, sample

import pygame
import pygame.gfxdraw
from pygame.math import Vector2

from voronoi import VoronoiDiagram


# Test data approximated from paper:
#
# points = [
#     Vector2(100, 150),  # 0
#     Vector2(120, 230),  # 1
#     Vector2(140, 160),  # 2
#     Vector2(170, 120),  # 3
#     Vector2(230, 170),  # 4
#     Vector2(250, 270),  # 5
#     Vector2(260, 230),  # 6
#     Vector2(330, 171),  # 7
# ]
#
# Note that Pygame's coordinate system has the origin in upper
# left corner but I'm using a coordinate system with the origin
# in the lower left corner so the Y coordinates needs to be
# reversed in order to reproduce the above figure (approximation
# of the figure in the sweep-line delaunay triangulation paper)

# Place points on a semi-regular grid.
WIDTH = 1010
HEIGHT = 900
N_X = 38
N_Y = 38
X_SD = 3
Y_SD = 3
X_OFFSET = 0
Y_OFFSET = 0
X_DELTA = 27
Y_DELTA = 25
X_SHIFT = X_DELTA / 2
Y_SHIFT = 0

points = []
for y in range(N_Y):
    x_offset = X_OFFSET + (y % 2) * X_SHIFT
    for x in range(N_X):
        y_offset = Y_OFFSET + ((x + 1) % 2) * Y_SHIFT
        points.append(
            Vector2(
                x * X_DELTA + gauss(0, X_SD) + x_offset,
                y * Y_DELTA + gauss(0, Y_SD) + y_offset,
            )
        )

# Instead place points by uniform random sampling?
if False:
    points = []
    for i in range(N_X * N_Y):
        x = randint(0, 1000)
        y = randint(0, 900)
        points.append(Vector2(x, y))

voronoi = VoronoiDiagram(points)
start_ix = voronoi.create_maze()


def main():
    draw_points = False
    draw_triangulation = False
    draw_circumcircles1 = False
    draw_circumcircles2 = False
    draw_voronoi_cells = False
    draw_voronoi_connections = True
    draw_maze_edges = True

    pygame.init()
    screen = pygame.display.set_mode((WIDTH, HEIGHT))

    s = pygame.Surface(screen.get_size(), pygame.SRCALPHA, 32)
    s = s.convert()
    s.fill((0, 0, 0))

    screen.blit(s, (0, 0))
    pygame.display.flip()

    try:
        while True:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    return
                if event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_ESCAPE or event.unicode == "q":
                        return
                    if event.key == pygame.K_1:
                        draw_points = not draw_points
                    if event.key == pygame.K_2:
                        draw_triangulation = not draw_triangulation
                    if event.key == pygame.K_3:
                        draw_circumcircles1 = not draw_circumcircles1
                    if event.key == pygame.K_4:
                        draw_circumcircles2 = not draw_circumcircles2
                    if event.key == pygame.K_5:
                        draw_voronoi_cells = not draw_voronoi_cells
                    if event.key == pygame.K_6:
                        draw_voronoi_connections = not draw_voronoi_connections
                    if event.key == pygame.K_7:
                        draw_maze_edges = not draw_maze_edges

            s.fill((0, 0, 0))

            if draw_triangulation:
                for tri in voronoi.triangles:
                    if tri.flagged_for_deletion:
                        continue
                    pygame.draw.aalines(
                        s,
                        (55, 55, 155),
                        True,
                        [
                            (
                                voronoi.points[tri.point_indices[0]].x,
                                voronoi.points[tri.point_indices[0]].y,
                            ),
                            (
                                voronoi.points[tri.point_indices[1]].x,
                                voronoi.points[tri.point_indices[1]].y,
                            ),
                            (
                                voronoi.points[tri.point_indices[2]].x,
                                voronoi.points[tri.point_indices[2]].y,
                            ),
                        ],
                        1,
                    )

            if draw_points:
                for p in voronoi.points:
                    pygame.gfxdraw.aacircle(
                        s, round(p.x), round(p.y), 0, (200, 200, 200)
                    )

            if draw_circumcircles2:
                for tri in voronoi.triangles:
                    if tri.flagged_for_deletion:
                        continue
                    if (
                        abs(tri.circumcenter.x) > 32767
                        or abs(tri.circumcenter.y) > 32767
                    ):
                        continue
                    pygame.gfxdraw.aacircle(
                        s,
                        round(tri.circumcenter.x),
                        round(tri.circumcenter.y),
                        round(tri.circumradius),
                        (100, 100, 100),
                    )

            if draw_voronoi_cells:
                # voronoi_count = 0
                for i, cell in enumerate(voronoi.cells):
                    # voronoi_count += 1
                    if not cell.color:
                        pygame.draw.aalines(
                            s, (0, 100, 0), True, [(p.x, p.y) for p in cell.polygon], 0
                        )
                    else:
                        pygame.gfxdraw.filled_polygon(
                            s, [(p.x, p.y) for p in cell.polygon], cell.color
                        )

                # print("Voronoi count:", voronoi_count)

                # triangle_count = 0
                # for i, tri in enumerate(voronoi.triangles):
                #     if tri.flagged_for_deletion:
                #         continue
                #     triangle_count += 1
                # print("Triangle count:", triangle_count)
                # print("Point count:", len(points))

            if draw_voronoi_connections:
                for i, cell in enumerate(voronoi.cells):
                    if i == start_ix:
                        pygame.gfxdraw.filled_polygon(
                            s, [(p.x, p.y) for p in cell.polygon], cell.color
                        )
                    for j, n_ix in enumerate(cell.neighbors):
                        if n_ix is not None and not cell.maze_edges[j] and i < n_ix:
                            ncell = voronoi.cells[n_ix]
                            t0, t1 = cell.edges[j]
                            p0 = voronoi.triangles[t0].circumcenter
                            p1 = voronoi.triangles[t1].circumcenter
                            p = p0 + (p1 - p0) / 2.0

                            pygame.draw.aalines(
                                s,
                                (0, 178, 178),
                                False,
                                [
                                    (cell.center.x, cell.center.y),
                                    (p.x, p.y),
                                    (ncell.center.x, ncell.center.y),
                                ],
                                1,
                            )
                    if abs(cell.center.x) > 32767 or abs(cell.center.y) > 32767:
                        continue
                    pygame.gfxdraw.aacircle(
                        s, round(cell.center.x), round(cell.center.y), 0, (0, 255, 255)
                    )

            if draw_circumcircles1:
                for tri in voronoi.triangles:
                    if tri.flagged_for_deletion:
                        continue
                    if (
                        abs(tri.circumcenter.x) > 32767
                        or abs(tri.circumcenter.y) > 32767
                    ):
                        continue
                    pygame.gfxdraw.aacircle(
                        s,
                        round(tri.circumcenter.x),
                        round(tri.circumcenter.y),
                        0,
                        (255, 0, 255),
                    )

            if draw_maze_edges:
                for cell in voronoi.cells:
                    for i in range(len(cell.edges)):
                        has_edge = cell.maze_edges[i]
                        t0, t1 = cell.edges[i]
                        p0 = voronoi.triangles[t0].circumcenter
                        p1 = voronoi.triangles[t1].circumcenter
                        if not has_edge:
                            continue
                        pygame.draw.aaline(
                            s, (255, 255, 255), (p0.x, p0.y), (p1.x, p1.y), 2
                        )

            screen.blit(s, (0, 0))
            pygame.display.flip()

    finally:
        pygame.quit()


if __name__ == "__main__":
    main()
