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
        point = (96663.12766666667, 439773.052)  # centroid of a triangle
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
        candidate = obj_5m / '37fz2_9.obj'
        obj_base = formats.factory.create('objmem')
        obj_candidate = formats.factory.create('objmem')
        obj_base.read(base)
        obj_candidate.read(candidate)
        base_pts = len(obj_base.points)
        base_stars = max(obj_base.stars)
        obj_base.add(obj_candidate)
        assert base_pts + len(obj_candidate.points) == len(obj_base.points)
        assert len(obj_base.points) == max(obj_base.stars) + 1

    def test_reindex(self, obj_5m):
        base = obj_5m / '37fz2_8.obj'
        obj_base = formats.factory.create('objmem')
        obj_base.read(base)
        base_pts = max(obj_base.stars)
        obj_base.reindex(start_id=base_pts)
        assert all(base_pts < v for v in obj_base.stars)

    def test_deduplicate_candidate(self, obj_5m):
        """Test deduplication on the candidate TIN.
        """
        def __point_id_is_sequential(points):
            start_id = min(points)
            for v in sorted(points):
                yield v == start_id
                start_id += 1

        base = obj_5m / '37fz2_8.obj'
        candidate = obj_5m / '37fz2_9.obj'
        obj_base = formats.factory.create('objmem')
        obj_candidate = formats.factory.create('objmem')
        obj_base.read(base)
        obj_candidate.read(candidate)
        base_pts = len(obj_base.points)
        base_stars = max(obj_base.stars)
        obj_candidate.reindex(base_stars)
        obj_base.deduplicate(obj_candidate, precision=3)
        # Expect that the candidate is reindexed so that the indices start
        # where the base indices end
        assert max(obj_base.stars) + 1 == min(obj_candidate.stars)
        # Expect that there is a 1-to-1 mapping between the points and the stars
        assert obj_base.stars.keys().isdisjoint(obj_base.points.keys()) is False
        assert obj_candidate.stars.keys().isdisjoint(obj_candidate.points.keys()) is False
        # Expect that the point indices are sequentially increment by 1
        assert all(__point_id_is_sequential(obj_base.points))
        assert all(__point_id_is_sequential(obj_candidate.points))

    def test_merge(self, obj_5m):
        infile = obj_5m / '37fz2_8.obj'
        infile_2 = obj_5m / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        obj.merge(obj_2, strategy='deduplicate', precision=3)
        assert obj.is_valid()

    def test_merge_return(self, obj_5m, data_dir):
        infile = obj_5m / '37fz2_8.obj'
        infile_2 = obj_5m / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        common_pts = obj.merge(obj_2, strategy='deduplicate', precision=3)
        with open(data_dir / '37fz2_8-9_merge_common-points-avg.csv',
                  'w') as fo:
            pointwriter = csv.writer(fo, delimiter='\t')
            for p in common_pts:
                pointwriter.writerow(p)
        assert obj.is_valid()

    def test_is_valid(self, obj_base):
        infile = obj_base / '37fz2_9.obj'
        obj = formats.OBJMem()
        obj.read(infile)
        assert obj.is_valid()

    def test_find_duplicate_points(self, obj_5m):
        infile = obj_5m / '37fz2_8.obj'
        infile_2 = obj_5m / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        duplicates = obj.find_duplicate_points(obj_2, precision=3)
        #TODO: better assertion could be done here
        assert len(duplicates) > 0

    def test_get_seam(self, obj_5m):
        infile = obj_5m / '37fz2_8.obj'
        infile_2 = obj_5m / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        duplicates = obj.find_duplicate_points(obj_2, precision=3)
        seam = obj.get_seam(duplicates)
        assert len(seam[0]) == len(seam[1])

    def test_add_seam(self, obj_5m):
        infile = obj_5m / '37fz2_8.obj'
        infile_2 = obj_5m / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        duplicates = obj.find_duplicate_points(obj_2, precision=3)
        seam = obj.get_seam(duplicates)
        maxid = max(obj.stars)
        candidate_stars = len(obj_2.stars)
        candidate_pts = len(obj_2.points)
        obj_2.add_seam(start_id=maxid, seam=seam)
        assert len(obj_2.stars) == len(seam[0]) + candidate_stars
        assert len(obj_2.points_dict) == len(seam[1]) + candidate_pts

