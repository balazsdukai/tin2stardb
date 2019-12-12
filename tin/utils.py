# -*- coding: utf-8 -*-

"""Various utility functions for handling geometry etc."""

import math
from statistics import mean
from typing import Tuple, Union, Iterable, Generator, Mapping
import logging

MODULE_MATPLOTLIB_AVAILABLE = True
try:
    import matplotlib.pyplot as plt
    import matplotlib.lines as lines
except ImportError as e:
    MODULE_MATPLOTLIB_AVAILABLE = False

log = logging.getLogger(__name__)


def __ccw__(vertices, star, link):
    """Sort the link in CounterClockWise order around the star"""
    x, y, z = 0, 1, 2
    localized = [(vertices[v][x] - vertices[star][x],
                  vertices[v][y] - vertices[star][y]) for v in link]
    rev_lookup = {localized[i]: a for i, a in enumerate(link)}
    return rev_lookup, sorted(localized, key=lambda p: math.atan2(p[1], p[0]))


def sort_ccw(vertices, stars) -> Generator:
    """Sort vertices in counter-clockwise order."""
    for star, link in stars.items():
        rev_lookup, ccw = __ccw__(vertices, star, link)
        yield star, [rev_lookup[co] for co in ccw]


def link_is_ccw(vertices, stars) -> Generator:
    """Check if the link of the star is ordered CounterClockWise."""
    for star, link in stars.items():
        rev_lookup, ccw = __ccw__(vertices, star, link)
        yield star, all(rev_lookup[co]==link[i] for i,co in enumerate(ccw))


def link_is_consistent(stars) -> Generator:
    """Checks if the links are consistent, thus vertex A is also present in the
    link of vertex B, if vertex B is in the link of vertex A."""
    for star, link in stars.items():
        yield star, all(star in stars[_star] for _star in link)


def triangle_is_consistent(stars, triangles) -> Generator:
    """Check that each adjacent triangle's vertices are consistent in it's star
    with the current one."""
    def __stars_are_consistent(tri, stars):
        for i, star in enumerate(tri):
            link = stars[star]
            try:
                # star is vertex 1 (v1) of the triangle (tri)
                idx_v2 = link.index(tri[i-2])
                idx_v3 = link.index(tri[i-1])
                _idx_v3 = idx_v2+1  if idx_v2+1 < len(link) else 0
                if _idx_v3 !=idx_v3:
                    pass
                yield _idx_v3 == idx_v3
            except ValueError:
                yield False
    for tri in triangles:
        yield tri, all(__stars_are_consistent(tri, stars))


def distance(a,b) -> float:
    """Distance between point a and point b"""
    x,y = 0,1
    return math.sqrt((a[x] - b[x])**2 + (a[y] - b[y])**2)


def orientation(a: Tuple[float, float], b: Tuple[float, float],
                c: Tuple[float, float]):
    """
    Determine if point (p) is LEFT, RIGHT, COLLINEAR with line segment (ab).

    :param a: Point 1
    :param b: Point 2
    :param c: Point which orientation to is determined with respect to (a,b)
    :return: 1 if (a,b,p) is CCW, 0 if p is collinear, -1 if (a,b,p) is CW

    >>> orientation((0.0, 0.0), (1.0, 0.0), (2.0, 0.0))
    0
    >>> orientation((0.0, 0.0), (1.0, 0.0), (0.5, 0.0))
    0
    >>> orientation((0.0, 0.0), (1.0, 0.0), (0.5, 1.0))
    1
    >>> orientation((0.0, 0.0), (1.0, 0.0), (0.5, -1.0))
    -1
    """
    x,y = 0,1
    re = ((a[x] - c[x]) * (b[y] - c[y])) - ((a[y] - c[y]) * (b[x] - c[x]))
    if re > 0:
        return 1
    elif re == 0:
        return 0
    else:
        return -1


def is_between(a,c,b) -> bool:
    """Return True if point c is on the segment ab

    Ref.: https://stackoverflow.com/a/328193
    """
    return math.isclose(distance(a,c) + distance(c,b), distance(a,b))


def in_bbox(tri: Tuple, bbox: Tuple) -> bool:
    """Evaluates if a triangle is in the provided bounding box.

    A triangle is in the BBOX if it's centorid is either completely within
    the BBOX, or overlaps with the South (lower) or West (left) boundaries
    of the BBOX.

    :param tri: A triangle defined as a tuple of three cooridnates of (x,y,z)
    :param bbox: Bounding Box as (minx, miny, maxx, maxy)
    """
    if not bbox or not tri:
        return False
    x,y,z = 0,1,2
    minx, miny, maxx, maxy = bbox
    # mean x,y,z coordinate of the triangle
    centroid = (mean(v[x] for v in tri),
                mean(v[y] for v in tri))
    within = ((minx < centroid[x] < maxx) and
              (miny < centroid[y] < maxy))
    on_south_bdry = is_between((minx, miny), centroid, (maxx, miny))
    on_west_bdry = is_between((minx, miny), centroid, (minx, maxy))
    return any((within, on_south_bdry, on_west_bdry))


def bbox(polygon) -> Tuple[float, float, float, float]:
    """Compute the Bounding Box of a polygon.

    :param polygon: List of coordinate pairs (x,y)
    """
    x,y = 0,1
    vtx = polygon[0]
    minx, miny, maxx, maxy = vtx[x], vtx[y], vtx[x], vtx[y]
    for vtx in polygon[1:]:
        if vtx[x] < minx:
            minx = vtx[x]
        elif vtx[y] < miny:
            miny = vtx[y]
        elif vtx[x] > maxx:
            maxx = vtx[x]
        elif vtx[y] > maxy:
            maxy = vtx[y]
    return minx, miny, maxx, maxy


def get_polygon(feature):
    """Get the polygon boundaries from a GeoJSON feature."""
    if not feature['geometry']['type'] == 'Polygon':
        log.warning(f"Feature ID {feature['properties']['id']} is not a Polygon")
    else:
        return feature['geometry']['coordinates'][0]


def find_side(polygon: Iterable[Tuple[float, ...]],
              neighbor: Iterable[Tuple[float, ...]],
              abs_tol: float = 0.0) ->\
        Union[Tuple[None, None],
              Tuple[str, Tuple[Tuple[float, float], Tuple[float, float]]]]:
    """Determines on which side does the neighbor polygon is located.

    .. warning::

        Assumes touching BBOXes of equal dimensions.

    :param polygon: The base polygon. A list of coordinate tuples.
    :param neighbor: The neighbor polygon.
    :param abs_tol: Absolute coordinate tolerance. Passed on to `:math.isclose`
    :returns: One of ['E', 'N', 'W', 'S'], the touching line segment
    """
    minx, miny, maxx, maxy = 0,1,2,3
    bbox_base = bbox(polygon)
    bbox_nbr = bbox(neighbor)
    if math.isclose(bbox_nbr[minx], bbox_base[maxx], abs_tol=abs_tol) \
        and math.isclose(bbox_nbr[miny], bbox_base[miny], abs_tol=abs_tol):
        return 'E', ((bbox_base[maxx], bbox_base[miny]), (bbox_base[maxx], bbox_base[maxy]))
    elif math.isclose(bbox_nbr[minx], bbox_base[minx], abs_tol=abs_tol) \
        and math.isclose(bbox_nbr[miny], bbox_base[maxy], abs_tol=abs_tol):
        return 'N', ((bbox_base[maxx], bbox_base[maxy]), (bbox_base[minx], bbox_base[maxy]))
    elif math.isclose(bbox_nbr[maxx], bbox_base[minx], abs_tol=abs_tol) \
        and math.isclose(bbox_nbr[maxy], bbox_base[maxy], abs_tol=abs_tol):
        return 'W', ((bbox_base[minx], bbox_base[maxy]), (bbox_base[minx], bbox_base[miny]), )
    elif math.isclose(bbox_nbr[maxx], bbox_base[maxx], abs_tol=abs_tol) \
        and math.isclose(bbox_nbr[maxy], bbox_base[miny], abs_tol=abs_tol):
        return 'S', ((bbox_base[minx], bbox_base[miny]), (bbox_base[maxx], bbox_base[miny]), )
    else:
        return None,None


def plot_star(vid, stars, vertices):
    """Plots the location of a vertex and its incident vertices in its link.

    :Example: plot_star(1, stars, vertices)

    :param vid: Vertex ID
    :param stars: List with the Link of the vertex
    :param vertices: List with vertex coordinates (used as lookup)
    :return: Plots a plot on screen
    """
    if not MODULE_MATPLOTLIB_AVAILABLE:
        raise ModuleNotFoundError("matplotlib is not installed, cannot plot")
    plt.clf()
    pts = [vertices[vid]] + [vertices[v] for v in stars[vid]]
    r = list(zip(*pts))
    plt.scatter(*r[0:2])
    labels = [vid] + stars[vid]
    # zip joins x and y coordinates in pairs
    for i, e in enumerate(labels):
        if e == vid:
            plt.annotate(e,  # this is the text
                         (pts[i][0], pts[i][1]),  # this is the point to label
                         textcoords="offset points",  # how to position the text
                         xytext=(0, 10),  # distance from text to points (x,y)
                         ha='center',
                         # horizontal alignment can be left, right or center
                         color='red')
        else:
            plt.annotate(e,  # this is the text
                         (pts[i][0], pts[i][1]),
                         textcoords="offset points",
                         xytext=(0, 10),
                         ha='center')
    plt.show()


def mean_coordinate(points: Iterable[Tuple]) -> Tuple[float, float]:
    """Compute the mean x- and y-coordinate from a list of points.

    :param points: An iterable of coordinate tuples where the first two elements
        of the tuple are the x- and y-coordinate respectively.
    :returns: A tuple of (mean x, mean y) coordinates
    """
    mean_x = mean(pt[0] for pt in points)
    mean_y = mean(pt[1] for pt in points)
    return mean_x, mean_y

# Computing Morton-code. Reference: https://github.com/trevorprater/pymorton ---

def __part1by1_64(n):
    """64-bit mask"""
    n &= 0x00000000ffffffff                  # binary: 11111111111111111111111111111111,                                len: 32
    n = (n | (n << 16)) & 0x0000FFFF0000FFFF # binary: 1111111111111111000000001111111111111111,                        len: 40
    n = (n | (n << 8))  & 0x00FF00FF00FF00FF # binary: 11111111000000001111111100000000111111110000000011111111,        len: 56
    n = (n | (n << 4))  & 0x0F0F0F0F0F0F0F0F # binary: 111100001111000011110000111100001111000011110000111100001111,    len: 60
    n = (n | (n << 2))  & 0x3333333333333333 # binary: 11001100110011001100110011001100110011001100110011001100110011,  len: 62
    n = (n | (n << 1))  & 0x5555555555555555 # binary: 101010101010101010101010101010101010101010101010101010101010101, len: 63

    return n


def __unpart1by1_64(n):
    n &= 0x5555555555555555                  # binary: 101010101010101010101010101010101010101010101010101010101010101, len: 63
    n = (n ^ (n >> 1))  & 0x3333333333333333 # binary: 11001100110011001100110011001100110011001100110011001100110011,  len: 62
    n = (n ^ (n >> 2))  & 0x0f0f0f0f0f0f0f0f # binary: 111100001111000011110000111100001111000011110000111100001111,    len: 60
    n = (n ^ (n >> 4))  & 0x00ff00ff00ff00ff # binary: 11111111000000001111111100000000111111110000000011111111,        len: 56
    n = (n ^ (n >> 8))  & 0x0000ffff0000ffff # binary: 1111111111111111000000001111111111111111,                        len: 40
    n = (n ^ (n >> 16)) & 0x00000000ffffffff # binary: 11111111111111111111111111111111,                                len: 32
    return n


def interleave(*args):
    """Interleave two integers"""
    if len(args) != 2:
        raise ValueError('Usage: interleave2(x, y)')
    for arg in args:
        if not isinstance(arg, int):
            print('Usage: interleave2(x, y)')
            raise ValueError("Supplied arguments contain a non-integer!")

    return __part1by1_64(args[0]) | (__part1by1_64(args[1]) << 1)


def deinterleave(n):
    if not isinstance(n, int):
        print('Usage: deinterleave2(n)')
        raise ValueError("Supplied arguments contain a non-integer!")

    return __unpart1by1_64(n), __unpart1by1_64(n >> 1)


def morton_code(x: float, y: float):
    """Takes an (x,y) coordinate tuple and computes their Morton-key.

    Casts float to integers by multiplying them with 100 (millimeter precision).
    """
    return interleave(int(x * 100), int(y * 100))


def rev_morton_code(morton_key: int) -> Tuple[float, float]:
    """Get the coordinates from a Morton-key"""
    x,y = deinterleave(morton_key)
    return float(x)/100.0, float(y)/100.0


# Compute tile range -----------------------------------------------------------
def tilesize(tin_paths) -> Tuple[float, float]:
    """Compute the tile size from Morton-codes for the input TINs.

    .. note:: Assumes regular grid.

    :returns: The x- and y-dimensions of a tile
    """
    centroids = []
    for i, morton_code in enumerate(tin_paths):
        if i == 2:
            break
        else:
            centroids.append(rev_morton_code(morton_code))
    return abs(centroids[0][0] - centroids[1][0]), abs(centroids[0][1] - centroids[1][1])


def __in_bbox__(point, range):
    """Check if a point is within a BBOX."""


def compute_8neighbors(tin_paths: Mapping, tilesize: Tuple) -> Generator:
    """Computes the 8 neighbors for the tiles, using a predefined search range."""
    for mc, filepath in tin_paths.items():
        center = rev_morton_code(mc)
        neighbours = set()
        # Search range
        minx = center[0] - tilesize[0]
        miny = center[1] - tilesize[1]
        maxx = center[0] + tilesize[0]
        maxy = center[1] + tilesize[1]
        for mc_nbr, filepath_nbr in tin_paths.items():
            if mc_nbr != mc:
                center_nbr = rev_morton_code(mc_nbr)
                within = ((minx < center_nbr[0] < maxx) and
                          (miny < center_nbr[1] < maxy))
                if within:
                    neighbours.update([filepath_nbr.with_suffix('.pickle'),])
        yield mc, (filepath, neighbours)
