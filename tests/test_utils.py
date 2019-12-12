#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests the utils module."""

import pytest
import logging
from pathlib import Path

from tin import utils, formats

log = logging.getLogger(__name__)


class TestBBOX:
    @pytest.mark.parametrize('triangle, result', [
        [((5, 0, 0), (3, -1, 0), (3, 1, 0)), True],
        [((5, 1, 0), (3, 0, 0), (3, 2, 0)), True],
        [((10, 0, 0), (13, -1, 0), (13, 1, 0)), False],
        [((10, 5, 0), (11, 7, 0), (9, 7, 0)), False],
    ])
    def test_in_bbox(self, triangle, result):
        bbox = (0, 0, 10, 10)
        assert utils.in_bbox(triangle, bbox) == result

    @pytest.mark.parametrize('polygon, bbox', [
        [((1.0, 4.0), (3.0,1.0), (6.0, 2.0), (6.0, 6.0), (2.0, 7.0)), (1.0, 1.0, 6.0, 7.0)],
        [((1.0, 4.0), (3.0,1.0), (6.0, 2.0), (6.0, 6.0), (2.0, 7.0), (1.0, 4.0)), (1.0, 1.0, 6.0, 7.0)]
    ])
    def test_bbox(self, polygon, bbox):
        assert utils.bbox(polygon) == bbox


class TestCCW:
    # FIXME: vertices will become dict instead of list

    @pytest.fixture(scope='class')
    def obj(self, obj_base):
        infile = obj_base / '37fz2_9.obj'
        obj = formats.OBJMem()
        obj.read(infile)
        pts = {v:point for v,point in enumerate(obj.points)}
        obj.points = pts
        yield obj

    def test_sort_ccw(self, obj_base):
        infile = obj_base / '37fz2_9.obj'
        obj = formats.OBJMem()
        vertices, adjacency_table = obj.parse_obj(infile)
        stars = utils.sort_ccw(vertices, adjacency_table)
        for vid, star in stars:
            # check for duplicates in the star
            assert len(star) == len(set(star))

    def test_sort_ccw_dict(self, obj_base):
        """Test when the vertices are stored in a dict instead of a list."""
        infile = obj_base / '37fz2_9.obj'
        obj = formats.OBJMem()
        _v, adjacency_table = obj.parse_obj(infile)
        vertices = {i: pt for i, pt in enumerate(_v)}
        stars = utils.sort_ccw(vertices, adjacency_table)
        for vid, star in stars:
            # check for duplicates in the star
            assert len(star) == len(set(star))

    def test_link_is_ccw(self, obj):
        res = utils.link_is_ccw(obj.points, obj.stars)
        assert all(ccw for vid, ccw in res)

    def test_link_is_consistent(self, obj):
        res = utils.link_is_consistent(obj.stars)
        assert all(consistent for vid, consistent in res)

    def test_triangle_is_consistent(self, obj):
        res = utils.triangle_is_consistent(obj.stars, obj.triangles())
        assert all(consistent for tri, consistent in res)


class TestSide:
    @pytest.mark.parametrize('poly, result', [
        [[(0.9, 0.5), (1.6, 0.5), (1.6, 0.9), (0.9, 0.9)],
         ('E', ((0.9, 0.5), (0.9, 0.9)))],
        [[(0.9, 0.9), (0.9, 1.6), (0.5, 1.6), (0.5, 0.9)],
         ('N', ((0.9, 0.9), (0.5, 0.9)))]
    ])
    def test_find_side(self, poly, result):
        base = [(0.5, 0.5), (0.9, 0.5), (0.9, 0.9), (0.5, 0.9)]
        side, segment = utils.find_side(base, poly)
        assert side == result[0]
        assert segment == result[1]

    def test_find_side_rd(self):
        base = [(96123.0, 440270.0), (97246.0, 440270.0), (97246.0, 441430.0), (96123.0, 441430.0)]
        neighbor = [(96123.0, 439110.0), (97246.0, 439110.0), (97246.0, 440269.969), (96123.0, 440269.969)]
        side, segment = utils.find_side(base, neighbor, abs_tol=0.1)
        assert side == 'S'

class TestSorting:
    @pytest.mark.parametrize('point', [
        (0, 0),
        (0.0, 0.0),
        (1.0, 1.0),
        (96663.25590546813, 439718.94288361823),
        (252914.232, 608211.603)
    ])
    def test_morton_code(self, point):
        utils.morton_code(*point)

    def test_rev_morton_code(self):
        point = (252914.232, 608211.603)
        morton_key = utils.morton_code(*point)
        point_res = utils.rev_morton_code(morton_key)
        assert pytest.approx(point[0], point_res[0]) and \
               pytest.approx(point[1], point_res[1])

class TestRange:
    def test_tilesize(self, obj_5m):
        tin_paths = {}
        for child in obj_5m.iterdir():
            if child.suffix == '.obj':
                vertices = formats.OBJ.parse_vertices(child)
                log.debug(f"Computing TIN centroids and Morton-key")
                center = utils.mean_coordinate(vertices)
                morton_key = utils.morton_code(*center)
                tin_paths[morton_key] = child
        cellsize = utils.tilesize(tin_paths)

    def test_compute_8neighbors(self, obj_5m):
        tin_paths = {}
        for child in obj_5m.iterdir():
            if child.suffix == '.obj':
                vertices = formats.OBJ.parse_vertices(child)
                log.debug(f"Computing TIN centroids and Morton-key")
                center = utils.mean_coordinate(vertices)
                morton_key = utils.morton_code(*center)
                tin_paths[morton_key] = child
        tin_paths = dict(utils.compute_8neighbors(tin_paths, (1300.0, 1300.0)))
        expectation_for_37fz2_6 = {'37fz2_1', '37fz2_7', '37fz2_2','37fz2_11', '37fz2_12'}
        _37fz2_6 = [p[1] for p in tin_paths.values() if p[0] == Path('/home/balazs/Development/tin/tests/data/obj/densified_5m/37fz2_6.obj')][0]
        assert _37fz2_6 == expectation_for_37fz2_6