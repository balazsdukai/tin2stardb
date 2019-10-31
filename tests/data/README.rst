Contents
---------

The TINs in OBJ format are generated with 3dfier. The following settings were
used for the `Terrain` in `3dfier <https://github.com/tudelft3d/3dfier>_`:

.. code-block::yaml

  Terrain:
    simplification: 100
    simplification_tinsimp: 0.1
    inner_buffer: 1.0
    use_LAS_classes:
      - 2
      - 9

Each TIN is a rectangular *tile*. A collection of tiles is a continuous
partitioning of an area. A tile has an index (or name), and an extent that is
 represented by a polygon.

`./obj/base` : A simple rectangle with 4 vertices (the tile extent) was used as
    `input_polygons` for 3dfier.

`./obj/densified_5m` : The tile extent was densified at 5m interval and used
as `input_polygons` for 3dfier. The point cloud is not buffered, so it has
the exact same extent as the tile.