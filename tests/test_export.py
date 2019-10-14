#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests the exporter module."""

import pytest

from tin import formats


@pytest.fixture('module')
def outfile(data_dir):
    yield data_dir / 'test_out.obj'


class TestExporter:

    def test_generate_triangles(self, tin_db, tin_schema):
        stardb = formats.StarDb(conn=tin_db, schema=tin_schema)
        triangles = stardb.triangles()
        assert len(list(triangles)) > 0

    def test_get_coordinates(self, tin_db, tin_schema):
        stardb = formats.StarDb(conn=tin_db, schema=tin_schema)
        coordinates = stardb.points()
        assert len(list(coordinates)) > 0

    def test_export_obj(self, tin_db, tin_schema, outfile):
        obj = formats.OBJDb(conn=tin_db, schema=tin_schema)
        obj.export(outfile)