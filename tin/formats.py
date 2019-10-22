# -*- coding: utf-8 -*-

"""File format handlers for in-database TIN storage."""

from __future__ import annotations
import re
import logging
from pathlib import Path
from typing import Tuple

from tin.star import Star, StarDb
from tin import utils

log = logging.getLogger(__name__)


class FormatFactory:
    """Registers and instantiates a Formatter.

    Formatters provide data format specific methods for importing/exporting
    a TIN from the Star-structure.
    """

    def __init__(self):
        self._formatters = {}

    def register_formatter(self, key, formatter):
        """Register a formatter for use.

        :param key: Name of the formatter
        :param formatter: Can be a function, a class, or an object that
            implements .__call__()
        """
        self._formatters[key] = formatter

    def create(self, key, **kwargs):
        """Instantiate a Processor"""
        formatter = self._formatters.get(key)
        if not formatter:
            raise ValueError(key)
        return formatter(**kwargs)


class OBJ(object):

    @staticmethod
    def parse_obj(path: str, bbox: Tuple[float, float, float, float] = None):
        """Import from Wavefront OBJ.

        Returns a list of vertices and an adjacency table. The adjacency table
        is a dict, storing { vertex: [incident vertices] }, which is the same
        structure as the Star, with the difference that [incident vertices] is
        not ordered CCW.

        .. note: Assumes that the OBJ file stores the list of vertices in one
            block (not mixed with other elements such as faces), and that the
            vertices are stored before the faces.

        :param path: Path to the OBJ file
        :param bbox: If Bounding Box is provided as (minx, miny, maxx, maxy),
            only the triangles that are within the BBOX are parsed.
        :return: (vertex list, adjacency table)
        """
        log.info(f"Parsing {path}")
        v_pat = re.compile(r"^v\s[\s\S]*")  # vertex
        f_pat = re.compile(r"^f\s[\s\S]*")  # face
        vertices = []
        incident = {}

        if bbox:
            log.info(f"Using BBOX {bbox} for subset")
        with open(path, 'r') as f_in:
            for line in f_in:
                v = v_pat.match(line)
                f = f_pat.match(line)
                if v:
                    # v has value like 'v 96126.383 440267.281 5.42\n'
                    vertices.append(
                        tuple(float(co) for co in v.group().strip().split()[1:])
                    )
                elif f:
                    # This software uses a 0-based index, so the OBJ input
                    # is changed to a 0-based too.
                    # f has value like 'f 1 2 3\n'
                    tri = [int(tid)-1 for tid in f.group().strip().split()[1:]]
                    if bbox:
                        tri_deref = tuple(vertices[vid] for vid in tri)
                        if not utils.in_bbox(tri_deref, bbox):
                            tri = None
                    if tri:
                        for vid in tri:
                            if vid in incident:
                                x = [k for k in tri if
                                     k != vid and k not in incident[vid]]
                                incident[vid] += x
                            else:
                                x = [k for k in tri if k != vid]
                                incident[vid] = x
        # reindex the vertices if they are subset
        if bbox:
            new_vertices = []
            adjacency_table = {}
            lookup = {old: new for new,old in enumerate(incident)}
            for old, new in lookup.items():
                adjacency_table[new] = [lookup[vid] for vid in incident[old]]
                new_vertices.append(vertices[old])
        else:
            new_vertices = vertices
            adjacency_table = incident
        return new_vertices, adjacency_table


class OBJMem(OBJ, Star):

    def read(self, path: str, bbox: Tuple[float] = None) -> None:
        """Read an OBJ into a Star-structure in memory.

        :param path:
        :param bbox:
        """
        self.points, adjacency_table = self.parse_obj(path, bbox=bbox)
        self.stars = dict(utils.sort_ccw(self.points, adjacency_table))

    def write(self, path: Path) -> None:
        """Write to Wavefront OBJ.

        :param path: :py:class:Path to the output file
        """
        with path.open(mode='w') as fo:
            log.info(f"Writing {path}")
            fo.write('# Converted from Star TIN structure to OBJ.\n')
            fo.write('# Tool written by Balázs Dukai, b.dukai@tudelft.nl\n')
            for point in self.points:
                fo.write('v %s %s %s\n' % point)
            for tri in self.triangles():
                # Star uses a 0-based index
                _tri = [t+1 for t in tri]
                fo.write(f"f {_tri[0]} {_tri[1]} {_tri[2]}\n")


class OBJDb(OBJ, StarDb):

    def insert(self, path, epsg, bbox):
        """Insert an OBJ into a Star-structure in the database."""
        vertices, adjacency_table = self.parse_obj(path, bbox=bbox)
        stars = utils.sort_ccw(vertices, adjacency_table)
        insert_gen = self._insert_generator(vertices, stars, epsg)
        self._insert_query(insert_gen)

    def export(self, path: Path):
        """Export to Wavefront OBJ.

        :param path: :py:class:Path to the output file
        """
        with path.open(mode='w') as fo:
            log.info(f"Writing {path}")
            fo.write('# Converted from Star TIN structure to OBJ.\n')
            fo.write('# Tool written by Balázs Dukai, b.dukai@tudelft.nl\n')
            for resultset in self.points():
                for vtx in resultset:
                    fo.write('v %s %s %s\n' % vtx)
            for tri in self.triangles():
                fo.write('f %s %s %s\n' % tri)


factory = FormatFactory()
factory.register_formatter('objdb', OBJDb)
factory.register_formatter('objmem', OBJMem)