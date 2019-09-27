# -*- coding: utf-8 -*-

"""Various utility functions for handling geometry etc."""

import math
from statistics import mean
from typing import Tuple
import logging

MODULE_MATPLOTLIB_AVAILABLE = True
try:
    import matplotlib.pyplot as plt
except ImportError as e:
    MODULE_MATPLOTLIB_AVAILABLE = False

log = logging.getLogger(__name__)


def distance(a,b) -> float:
    """Distance between point a and point b"""
    x,y = 0,1
    return math.sqrt((a[x] - b[x])**2 + (a[y] - b[y])**2)


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


def bbox(polygon) -> Tuple:
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
