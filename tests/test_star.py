# -*- coding: utf-8 -*-

"""Tests the star module."""

import pytest
import logging

from tin import formats

log = logging.getLogger(__name__)

class TestStar:
    def test_is_valid(self, obj_base):
        infile = obj_base / '37fz2_9.obj'
        obj = formats.OBJMem()
        obj.read(infile)
        assert obj.is_valid()