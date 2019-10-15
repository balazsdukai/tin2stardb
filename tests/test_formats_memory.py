#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests the in-memory classes from the formats module."""

import pytest
import logging
from pathlib import Path

from tin import formats

log = logging.getLogger(__name__)


@pytest.fixture('module', params=['37fz2_8.obj', '37fz2_9.obj',
                                  '37fz2_13.obj', '37fz2_14.obj'])
def infile(data_dir, request):
    yield data_dir / 'obj' / request.param


class TestOBJ:
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
        for center, link in obj.stars:
            # check for duplicates in the link
            assert len(link) == len(set(link))

    def test_write(self, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj.write(Path('/tmp/tin_test.obj'))

    def test_merge(self, data_dir):
        base = data_dir / 'obj' / '37fz2_8.obj'
        neighbor = data_dir / 'obj' / '37fz2_9.obj'
        obj_base = formats.factory.create('objmem')
        obj_neighbor = formats.factory.create('objmem')
        obj_base.read(base)
        obj_neighbor.read(neighbor)
        base_pts = len(obj_base.points)
        base_stars = max(obj_base.stars)
        obj_base.merge(obj_neighbor)
        assert base_pts + len(obj_neighbor.points) - 1 == len(obj_base.points)
        assert len(obj_base.points) - 1 == max(obj_base.stars)
        assert base_stars + max(obj_neighbor.stars) == max(obj_base.stars)
