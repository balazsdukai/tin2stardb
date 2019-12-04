# -*- coding: utf-8 -*-
# Copyright:    (C) 2019 by Balázs Dukai, Delft University of Technology
# Begin:        2019-12-03
# Email:        b.dukai@tudelft.nl

"""Test streaming strategies"""

import logging

from tin import formats

log = logging.getLogger(__name__)

class TestStreamingMerge:

    def test_find_duplicate_points(self, obj_5m):
        infile = obj_5m / '37fz2_8.obj'
        infile_2 = obj_5m / '37fz2_9.obj'
        obj = formats.factory.create('objmem')
        obj.read(infile)
        obj_2 = formats.factory.create('objmem')
        obj_2.read(infile_2)
        duplicates = obj.find_duplicate_points(obj_2)
        #TODO: better assertion could be done here
        assert len(duplicates) > 0

    def test_streaming_merge(self):
        assert False