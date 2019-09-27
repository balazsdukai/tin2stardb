# -*- coding: utf-8 -*-

"""File format handlers for in-memory TIN storage."""

import re
import logging
from pathlib import Path
from typing import Tuple

from tin.star import Star
from tin.formats import FormatFactory
from tin import utils

log = logging.getLogger(__name__)


class OBJ(Star):

    def read(self, path, bbox=None):
        """Read an OBJ into a Star-structure in memory.

        :param path:
        :param bbox:
        """
        self.points, adjacency_table = self._parse_obj(path, bbox=bbox)
        self.stars = self._sort_ccw(self.points, adjacency_table)

    def write(self, path: Path):
        """Write to Wavefront OBJ.

        :param path: :py:class:Path to the output file
        """
        with path.open(mode='w') as fo:
            log.info(f"Writing {path}")
            fo.write('# Converted from Star TIN structure to OBJ.\n')
            fo.write('# Tool written by Balázs Dukai, b.dukai@tudelft.nl\n')
            # The first element in self.points is None, to fake a 1-based index
            for point in self.points[1:]:
                fo.write('v %s %s %s\n' % point)
            for tri in self.triangles():
                fo.write('f %s %s %s\n' % tri)

    @staticmethod
    def _parse_obj(path: str, bbox: Tuple = None):
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
        vertices = [None]  # because OBJ has 1-based indexing
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
                    # f has value like 'f 1 2 3\n'
                    tri = [int(tid) for tid in f.group().strip().split()[1:]]
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
            new_vertices = [None] # because OBJ has 1-based indexing
            adjacency_table = {}
            lookup = {old: new+1 for new,old in enumerate(incident)}
            for old, new in lookup.items():
                adjacency_table[new] = [lookup[vid] for vid in incident[old]]
                new_vertices.append(vertices[old])
        else:
            new_vertices = vertices
            adjacency_table = incident
        return new_vertices, adjacency_table


factory_memory = FormatFactory()
factory_memory.register_formatter('obj', OBJ)