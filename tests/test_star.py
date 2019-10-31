# -*- coding: utf-8 -*-

"""Tests the star module."""

import logging
from math import isclose
from statistics import mean
import csv

from tin import formats

log = logging.getLogger(__name__)

class TestStar:

    def test_pointlocation(self, obj_base):
        infile = obj_base / '37fz2_9.obj'
        point = (96663.12766666667, 439773.052) # centroid of a triangle
        obj = formats.factory.create('objmem')
        obj.read(infile)
        tri = obj.pointlocation(point)
        x_ctr = mean(obj.points[t][0] for t in tri)
        y_ctr = mean(obj.points[t][1] for t in tri)
        # check if the located triangle's centroid and the search point (same
        # centroid) match
        assert isclose(point[0], x_ctr)
        assert isclose(point[1], y_ctr)

    def test_straightwalk(self, obj_base):
        line = [(97246.0, 441430.0), (96123.0, 441430.0)]
        infile = obj_base / '37fz2_8.obj'
        infile_2 = obj_base / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        triangles = obj.straight_walk(*line)
        print(triangles)
        obj.add(obj_2)
        triangles_2 = obj.straight_walk(*line)
        print(triangles)
        assert triangles == triangles_2

    def test_add(self, obj_5m):
        base = obj_5m / '37fz2_8.obj'
        neighbor = obj_5m / '37fz2_9.obj'
        obj_base = formats.factory.create('objmem')
        obj_neighbor = formats.factory.create('objmem')
        obj_base.read(base)
        obj_neighbor.read(neighbor)
        base_pts = len(obj_base.points)
        base_stars = max(obj_base.stars)
        obj_base.add(obj_neighbor)
        assert base_pts + len(obj_neighbor.points) == len(obj_base.points)
        assert len(obj_base.points) == max(obj_base.stars)+1

    def test_merge(self, obj_5m):
        infile = obj_5m / '37fz2_8.obj'
        infile_2 = obj_5m / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        obj.merge(obj_2, strategy='deduplicate', precision=3)
        assert obj.is_valid()

    def test_merge_return(self, obj_5m,data_dir):
        infile = obj_5m / '37fz2_8.obj'
        infile_2 = obj_5m / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        common_pts = obj.merge(obj_2, strategy='deduplicate', precision=3,
                               ret_common_pts=True)
        with open(data_dir / '37fz2_8-9_merge_common-points-avg.csv', 'w') as fo:
            pointwriter = csv.writer(fo, delimiter='\t')
            for p in common_pts:
                pointwriter.writerow(p)
        assert obj.is_valid()

    def test_is_valid(self, obj_base):
        infile = obj_base / '37fz2_9.obj'
        obj = formats.OBJMem()
        obj.read(infile)
        assert obj.is_valid()