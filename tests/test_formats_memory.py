#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests the formats_memory module."""

import pytest
import logging
from pathlib import Path

from tin import formats_memory

log = logging.getLogger(__name__)


@pytest.fixture('module', params=['37fz2_8.obj', '37fz2_9.obj',
                                  '37fz2_13.obj', '37fz2_14.obj'])
def infile(data_dir, request):
    yield data_dir / 'obj' / request.param


class TestOBJ:
    def test_parse_obj(self, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats_memory.factory_memory.create('obj')
        points, adjacency_table = obj._parse_obj(infile)
        for center, link in adjacency_table.items():
            # check for duplicates in the link
            assert len(link) == len(set(link))

    def test_read(self, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats_memory.factory_memory.create('obj')
        obj.read(infile)
        for center, link in obj.stars:
            # check for duplicates in the link
            assert len(link) == len(set(link))

    def test_write(self, data_dir):
        infile = data_dir / 'obj' / '37fz2_9.obj'
        obj = formats_memory.factory_memory.create('obj')
        obj.read(infile)
        obj.write(Path('/tmp/tin_test.obj'))