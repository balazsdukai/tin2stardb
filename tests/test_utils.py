#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests the utils module."""

import pytest
import logging

from tin import utils

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