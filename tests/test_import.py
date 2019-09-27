#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests the formats module."""

import pytest
import logging
from itertools import chain

from tin import formats

log = logging.getLogger(__name__)


@pytest.fixture('module', params=['37fz2_8.obj', '37fz2_9.obj',
                                  '37fz2_13.obj', '37fz2_14.obj'])
def infile(data_dir, request):
    yield data_dir / 'obj' / request.param


class TestOBJ:
    def test_duplicate_points(self, tin_db, tin_schema, data_dir):
        """Check if there are duplicate points in the input file"""
        infile = data_dir / 'obj' / '37fz2_9.obj'
        # infile = '/home/balazs/Development/TIN_scale_up/37fz2_8.obj'
        obj = formats.OBJ(conn=tin_db, schema=tin_schema)
        vertices, adjacency_table = obj._parse_obj(infile)
        coordinate_hash = {}
        for center, star in adjacency_table.items():
            # if center == 67:
            #     formats.plot_star(67, adjacency_table, vertices)
            for v in adjacency_table[center]:
                if vertices[v] in coordinate_hash \
                        and coordinate_hash[vertices[v]] != v:
                    print(f"Found duplicate: {v} and "
                          f"{coordinate_hash[vertices[v]]} with {vertices[v]}")
                else:
                    coordinate_hash[vertices[v]] = v

    def test_parse_obj(self, tin_db, tin_schema, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats.OBJ(conn=tin_db, schema=tin_schema)
        vertices, adjacency_table = obj._parse_obj(infile)
        for center, star in adjacency_table.items():
            # check for duplicates in the star
            assert len(star) == len(set(star))

    def test_sort_ccw(self, tin_db, tin_schema, infile, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats.OBJ(conn=tin_db, schema=tin_schema)
        vertices, adjacency_table = obj._parse_obj(infile)
        stars = obj._sort_ccw(vertices, adjacency_table)
        for vid, star in stars:
            # check for duplicates in the star
            assert len(star) == len(set(star))

    def test_parse_obj_bbox(self, data_dir, tin_schema, tin_db):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        bbox = (96364.2115377, 439370.11094742,
                96939.83980984, 440030.90870881)
        obj = formats.OBJ(conn=tin_db, schema=tin_schema)
        vertices, adjacency_table = obj._parse_obj(infile, bbox=bbox)
        assert set(adjacency_table.keys()) == set(chain(*adjacency_table.values()))
        assert set(range(1, len(vertices))) == set(adjacency_table.keys())
