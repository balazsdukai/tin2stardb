#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests the in-memory classes from the formats module."""

import pytest
import logging
from pathlib import Path

from tin import formats, schema

log = logging.getLogger(__name__)


@pytest.fixture(scope='module', params=['37fz2_8.obj', '37fz2_9.obj',
                                  '37fz2_13.obj', '37fz2_14.obj'])
def infile(obj_base, request):
    yield obj_base / 'obj' / request.param

class TestOBJMem:
    def test_parse_obj(self, obj_base):
        infile = obj_base / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        points, adjacency_table = obj.parse_obj(infile)
        for center, link in adjacency_table.items():
            # check for duplicates in the link
            assert len(link) == len(set(link))

    def test_read(self, obj_base):
        infile = obj_base / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        for center, link in obj.stars.items():
            # check for duplicates in the link
            assert len(link) == len(set(link))

    def test_write(self, obj_base):
        infile = obj_base / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj.write(Path('/tmp/tin_test.obj'))


class TestOBJDb:
    def test_import(self, tin_db, tin_schema, infile, cfg):
        # TODO: need to factor out the schema creation and set up fixtures for it
        schema.create_relations(tin_db, tin_schema, cfg['epsg'])
        obj = formats.factory.create('objdb', conn=tin_db, schema=tin_schema)
        obj.insert(infile, cfg['epsg'], bbox=None)

    def test_export(self, tin_db, tin_schema, cfg, obj_base):
        obj = formats.factory.create('objdb', conn=tin_db, schema=tin_schema)
        obj.export(obj_base / 'test_db_out.obj')