#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests the in-memory classes from the formats module."""

import pytest
import logging
from pathlib import Path
from statistics import mean
from math import isclose

from tin import formats, schema

log = logging.getLogger(__name__)


@pytest.fixture(scope='module', params=['37fz2_8.obj', '37fz2_9.obj',
                                  '37fz2_13.obj', '37fz2_14.obj'])
def infile(data_dir, request):
    yield data_dir / 'obj' / request.param

class TestOBJMem:
    def test_parse_obj(self, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        points, adjacency_table = obj.parse_obj(infile)
        for center, link in adjacency_table.items():
            # check for duplicates in the link
            assert len(link) == len(set(link))

    def test_read(self, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        for center, link in obj.stars.items():
            # check for duplicates in the link
            assert len(link) == len(set(link))

    def test_write(self, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj.write(Path('/tmp/tin_test.obj'))

    def test_add(self, data_dir):
        base = data_dir / 'obj' / '37fz2_8.obj'
        neighbor = data_dir / 'obj' / '37fz2_9.obj'
        obj_base = formats.factory.create('objmem')
        obj_neighbor = formats.factory.create('objmem')
        obj_base.read(base)
        obj_neighbor.read(neighbor)
        base_pts = len(obj_base.points)
        base_stars = max(obj_base.stars)
        obj_base.add(obj_neighbor)
        assert base_pts + len(obj_neighbor.points) - 1 == len(obj_base.points)
        assert len(obj_base.points) - 1 == max(obj_base.stars)
        assert base_stars + max(obj_neighbor.stars) == max(obj_base.stars)

    def test_pointlocation(self, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
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

    def test_straightwalk(self, data_dir):
        line = [(97246.0, 441430.0), (96123.0, 441430.0)]
        infile = data_dir / 'obj' / '37fz2_8.obj'
        infile_2 = data_dir / 'obj' / '37fz2_9.obj'
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

    def test_merge(self, data_dir):
        infile = data_dir / 'obj' / '37fz2_8.obj'
        infile_2 = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        obj.merge(obj_2, strategy='deduplicate')


class TestOBJDb:
    def test_import(self, tin_db, tin_schema, infile, cfg):
        # TODO: need to factor out the schema creation and set up fixtures for it
        schema.create_relations(tin_db, tin_schema, cfg['epsg'])
        obj = formats.factory.create('objdb', conn=tin_db, schema=tin_schema)
        obj.insert(infile, cfg['epsg'], bbox=None)

    def test_export(self, tin_db, tin_schema, cfg, data_dir):
        obj = formats.factory.create('objdb', conn=tin_db, schema=tin_schema)
        obj.export(data_dir / 'test_db_out.obj')