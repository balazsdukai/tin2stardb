# -*- coding: utf-8 -*-
# Copyright:    (C) 2019 by Balázs Dukai, Delft University of Technology
# Begin:        30-10-19
# Email:        b.dukai@tudelft.nl

import pytest
import csv
import logging
from pathlib import Path
from click.testing import CliRunner

from tin import main
from tin import formats

log = logging.getLogger(__name__)

# @pytest.mark.balazs_local
class TestIntegration:

    @pytest.fixture(scope='class')
    def obj_5m_buffer(self):
        return Path('/data/TIN_scale_up/results/densified_5m_buff/clipped')

    @pytest.fixture(scope='class')
    def obj_5m(self):
        return Path('/data/TIN_scale_up/results/densified_5m/AHN3')

    def test_merge_densified_5m(self, obj_5m, data_dir):
        """Merge multiple TINs in a row-wise order. Densified tile polygons.
        Unbuffered point cloud.

        The input tiles have densified vertices on the tile boundary.
        Densification is done with QGIS with the *Densify by interval*
        processing toolbox, with an interval parameter of *5* (5m).
        The point cloud is not buffered, so it has the exact same extent as
        the tile.
        """
        common_points = []
        quality = []
        row_wise_order = [
            '37fz2_1.obj', '37fz2_6.obj', '37fz2_11.obj',
            '37fz2_2.obj', '37fz2_7.obj', '37fz2_12.obj',
            '37fz2_3.obj', '37fz2_8.obj', '37fz2_13.obj',
            '37fz2_4.obj', '37fz2_9.obj', '37fz2_14.obj',
        ]
        base_obj = formats.factory.create('objmem')
        base_file = row_wise_order.pop(0)
        base_obj.read(obj_5m / base_file)
        for neighbor_file in row_wise_order:
            neighbor = formats.factory.create('objmem')
            neighbor.read(obj_5m / neighbor_file)
            # Remember that Star.merge() adds the TIN to the current one
            cp,_qu = base_obj.merge(neighbor, strategy='deduplicate',
                                    precision=3, return_quality=True)
            qu = [neighbor_file,]
            qu.extend(_qu)
            common_points.extend(cp)
            quality.append(qu)
            if not base_obj.is_valid():
                log.warning(f"After merging {neighbor_file} the TIN is not valid")
        # Write results
        base_obj.write(data_dir / 'densified_5m_merge.obj')
        with open(data_dir / 'densified_5m_merge_common-points-avg.csv', 'w') as fo:
            pointwriter = csv.writer(fo, delimiter='\t')
            for p in common_points:
                pointwriter.writerow(p)
        with open(data_dir / 'densified_5m_merge_quality.csv', 'w') as fo:
            pointwriter = csv.writer(fo, delimiter='\t')
            for q in quality:
                pointwriter.writerow(q)

    def test_merge_densified_5m_buffer(self, obj_5m_buffer, data_dir):
        """Merge multiple TINs in a row-wise order. Densified tile polygons.
        Buffered point cloud with 200 meters.

        The input tiles have densified vertices on the tile boundary.
        Densification is done with QGIS with the *Densify by interval*
        processing toolbox, with an interval parameter of *5* (5m).
        The point cloud is buffered with 200 meters, so it has a larger extent
        than the tile.
        """
        common_points = []
        quality = []
        row_wise_order = [
            '37fz2_1.obj', '37fz2_6.obj', '37fz2_11.obj',
            '37fz2_2.obj', '37fz2_7.obj', '37fz2_12.obj',
            '37fz2_3.obj', '37fz2_8.obj', '37fz2_13.obj',
            '37fz2_4.obj', '37fz2_9.obj', '37fz2_14.obj',
        ]
        base_obj = formats.factory.create('objmem')
        base_file = row_wise_order.pop(0)
        base_obj.read(obj_5m_buffer / base_file)
        for neighbor_file in row_wise_order:
            neighbor = formats.factory.create('objmem')
            neighbor.read(obj_5m_buffer / neighbor_file)
            # Remember that Star.merge() adds the TIN to the current one
            cp, _qu = base_obj.merge(neighbor, strategy='deduplicate',
                                     precision=3, return_quality=True)
            qu = [neighbor_file, ]
            qu.extend(_qu)
            common_points.extend(cp)
            quality.append(qu)
            if not base_obj.is_valid():
                log.warning(
                    f"After merging {neighbor_file} the TIN is not valid")
        # Write out results
        base_obj.write(data_dir / 'densified_5m_buffer_merge.obj')
        with open(data_dir / 'densified_5m_buffer_merge_common-points-avg.csv',
                  'w') as fo:
            pointwriter = csv.writer(fo, delimiter='\t')
            for p in common_points:
                pointwriter.writerow(p)
        with open(data_dir / 'densified_5m_buffer_merge_quality.csv', 'w') as fo:
            pointwriter = csv.writer(fo, delimiter='\t')
            for q in quality:
                pointwriter.writerow(q)

    def test_merge_densified_5m_stream(self, obj_5m, data_dir):
        """Merge multiple TINs in a stream in a Morton-order. Densified tile polygons.
        Buffered point cloud with 200 meters.

        The input tiles have densified vertices on the tile boundary.
        Densification is done with QGIS with the *Densify by interval*
        processing toolbox, with an interval parameter of *5* (5m).
        The point cloud is buffered with 200 meters, so it has a larger extent
        than the tile.
        """
        runner = CliRunner()
        result = runner.invoke(main.main, [
            '-vv',
            str(data_dir / 'tin_config.yml'),
            'merge',
            str(data_dir / 'obj' / 'densified_5m'),
            'all_merge.obj'
        ])
        print(result.output)
