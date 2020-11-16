# -*- coding: utf-8 -*-

"""Main app."""

from os import listdir
from pathlib import Path
from sys import stdout
import logging
import json
import pickle

import yaml
import click

from tin import schema, db, formats, utils


def parse_config(config):
    return yaml.load(config, Loader=yaml.FullLoader)


def configure_logging(verbosity):
    """Configures the general logging in the application"""
    log_level = max(10, 30 - 10 * verbosity)
    logging.basicConfig(
        stream=stdout,
        level=log_level,
        format='%(asctime)s\t%(name)-24s\t%(lineno)s\t[%(levelname)-8s]\t%(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


@click.group()
@click.argument('configuration', type=click.File('r'))
@click.option(
    '--verbose', '-v',
    count=True,
    help="Increase verbosity. You can increment the level by chaining the "
         "argument, eg. -vvv")
@click.option(
    '--quiet', '-q',
    count=True,
    help="Decrease verbosity.")
@click.pass_context
def main(ctx, configuration, verbose, quiet):
    """Read/Write a TIN to a database from files.

    CONFIGURATION is a yaml file that contains the database credentials and
    table configuration options.

    """
    ctx.ensure_object(dict)
    ctx.obj['cfg'] = parse_config(configuration)
    # For logging from the click commands
    verbosity = verbose - quiet
    configure_logging(verbosity)
    ctx.obj['log'] = logging.getLogger(__name__)
    return 0


@click.command('import')
@click.argument('filename', type=str)
@click.option('--bbox', nargs=4, type=float, default=None,
              help='Import only a subset of the TIN in the 2D bbox: '
                   '(minx miny maxx maxy).')
@click.option('--bboxes', type=click.File(), default=None,
              help='A GeoJSON file containing polygons of bounding boxes. Each '
                   'Feature must have a "tile" property, which is used for '
                   'finding the TIN file with equivalent name. Only used when '
                   'a directory of TIN files are provided as input.')
@click.pass_context
def import_cmd(ctx, filename, bbox, bboxes):
    """Import TIN from a file into a database, using a Star-structure.

    FILENAME can be a path to single file, or a directory containing the files
    to be imported. The importer tries to guess the file format from its suffix.

    Supported file formats: Wavefront OBJ (.obj)
    """
    path = Path(filename).resolve()
    conn = db.Db(**ctx.obj['cfg']['database'])
    tin_schema = db.Schema(ctx.obj['cfg']['tin_schema'])
    log = ctx.obj['log']
    try:
        schema.create_relations(conn, tin_schema, ctx.obj['cfg']['epsg'])
        if path.is_file():
            # path.suffix returns '.obj' from some/path/file.obj
            f = path.suffix.lower().split('.')[1] # suffix without leading dot
            if f == 'obj':
                # FIXME: needs objdb
                fmt = 'objmem'
            else:
                raise click.exceptions.BadParameter(f"Unsupported format {f}")
            formatter = formats.factory.create(fmt, conn=conn,
                                               schema=tin_schema)
            formatter.insert(path, ctx.obj['cfg']['epsg'], bbox=bbox)
            formatter.create_index()
        elif path.is_dir():
            suffixes = set([Path(f).suffix.lower() for f in listdir(filename)
                            if (path / f).is_file()])
            files = [path / f for f in listdir(filename)
                     if (path / f).is_file()]
            if len(files) == 0:
                click.exceptions.FileError(f"Did not find any file in {path}.")
            if len(suffixes) > 1:
                click.exceptions.FileError(f"Found the {suffixes} filetypes "
                                           f"in {path}. It is not clear which"
                                           f" to import. Please use only "
                                           f"a single format in a directory.")
            suffix = suffixes.pop()
            f = suffix.lower().split('.')[1]
            if f == 'obj':
                # FIXME: needs objdb
                fmt = 'objdb'
            else:
                raise click.exceptions.BadParameter(f"Unsupported format {f}")
            formatter = formats.factory.create(fmt, conn=conn,
                                               schema=tin_schema)
            geojson = json.load(bboxes) if bboxes else None
            for file_path in files:
                tile = file_path.stem
                if geojson:
                    polygon = [utils.get_polygon(f) for f in geojson['features']
                               if f['properties']['tile'] == tile]
                    if len(polygon) == 0:
                        log.warning(f"Couldn't match any GeoJSON feature to "
                                    f"tile {tile}")
                        bbox = None
                    elif len(polygon) > 1:
                        log.warning(f"More than 1 matching feature to tile"
                                    f" {tile}. Only using the first.")
                        bbox = utils.bbox(polygon[0])
                    else:
                        bbox = utils.bbox(polygon[0])
                formatter.insert(file_path, ctx.obj['cfg']['epsg'], bbox=bbox)
            formatter.create_index()
        else:
            raise FileNotFoundError(path)
    except Exception as e:
        raise
    finally:
        conn.close()


@click.command('export')
@click.argument('filename', type=str)
@click.option('--format',
              type=click.Choice(['obj']),
              help="Export format.")
@click.pass_context
def export_cmd(ctx, filename, format):
    """Export a TIN from a Star TIN database model to a file."""
    path = Path(filename).resolve()
    if not Path(path.parent).exists():
        raise NotADirectoryError(f"Directory {path.parent} not exists")
    conn = db.Db(**ctx.obj['cfg']['database'])
    tin_schema = db.Schema(ctx.obj['cfg']['tin_schema'])
    try:
        formatter = formats.factory.create(format, conn=conn, schema=tin_schema)
        formatter.export(path)
    except:
        raise
    finally:
        conn.close()


@click.command('delete')
@click.pass_context
def delete_cmd(ctx):
    """Delete a TIN database."""
    conn = db.Db(**ctx.obj['cfg']['database'])
    tin_schema = db.Schema(ctx.obj['cfg']['tin_schema'])
    try:
        schema.drop_relations(conn, tin_schema)
    except:
        raise
    finally:
        conn.close()


@click.command('subset')
@click.argument('infile', type=str)
@click.argument('outfile', type=str)
@click.option('--bbox', nargs=4, type=float, default=None,
              help='Import only a subset of the TIN in the 2D bbox: '
                   '(minx miny maxx maxy).')
@click.option('--bboxes', type=click.File(), default=None,
              help='A GeoJSON file containing polygons of bounding boxes. Each '
                   'Feature must have a "tile" property, which is used for '
                   'finding the TIN file with equivalent name. Only used when '
                   'a directory of TIN files are provided as input.')
@click.pass_context
def subset_cmd(ctx, infile, outfile, bbox, bboxes):
    """Subset a TIN from a file and write it out to a file.

    INFILE can be a path to single file, or a directory containing the files
    to be imported. The importer tries to guess the file format from its suffix.

    OUTFILE can be a path to a single file, or a directory.

    Supported file formats: Wavefront OBJ (.obj)
    """
    path = Path(infile).resolve()
    outpath = Path(outfile).resolve()
    log = ctx.obj['log']
    try:
        if path.is_file():
            # path.suffix returns '.obj' from some/path/file.obj
            f = path.suffix.lower().split('.')[1] # suffix without leading dot
            if f == 'obj':
                fmt = 'objmem'
            else:
                raise click.exceptions.BadParameter(f"Unsupported format {f}")
            formatter = formats.factory.create(fmt)
            formatter.read(path, bbox=bbox)
            formatter.write(outpath)
        elif path.is_dir():
            outpath.mkdir(exist_ok=True)
            suffixes = set([Path(f).suffix.lower() for f in listdir(infile)
                            if (path / f).is_file()])
            files = [path / f for f in listdir(infile)
                     if (path / f).is_file()]
            if len(files) == 0:
                click.exceptions.FileError(f"Did not find any file in {path}.")
            if len(suffixes) > 1:
                click.exceptions.FileError(f"Found the {suffixes} filetypes "
                                           f"in {path}. It is not clear which"
                                           f" to import. Please use only "
                                           f"a single format in a directory.")
            suffix = suffixes.pop()
            f = suffix.lower().split('.')[1]
            if f == 'obj':
                fmt = 'objmem'
            else:
                raise click.exceptions.BadParameter(f"Unsupported format {f}")
            formatter = formats.factory.create(fmt)
            geojson = json.load(bboxes) if bboxes else None
            for file_path in files:
                tile = file_path.stem
                if geojson:
                    polygon = [utils.get_polygon(f) for f in geojson['features']
                               if f['properties']['tile'] == tile]
                    if len(polygon) == 0:
                        log.warning(f"Couldn't match any GeoJSON feature to "
                                    f"tile {tile}")
                        bbox = None
                    elif len(polygon) > 1:
                        log.warning(f"More than 1 matching feature to tile"
                                    f" {tile}. Only using the first.")
                        bbox = utils.bbox(polygon[0])
                    else:
                        bbox = utils.bbox(polygon[0])
                formatter.read(file_path, bbox=bbox)
                formatter.write(outpath / file_path.name)
        else:
            raise click.exceptions.FileError(str(path))
    except Exception as e:
        raise


@click.command('merge')
@click.argument('directory', type=click.Path(exists=True))
@click.argument('outfile', type=str)
@click.option('--in_memory', is_flag=True)
@click.pass_context
def merge_cmd(ctx, directory, outfile, in_memory):
    """Merge TINs in OBJ format.

    DIRECTORY is the folder that contains the TINs to be merged. All .obj files
    in this folder are merged together into a single file.

    OUTFILE is the name of the output file. The OUTFILE is written into
    DIRECTORY.
    """
    log = ctx.obj['log']
    tin_paths = {} # Will store { morton code : (file path, (range)) }
    dirpath = Path(directory).resolve()
    outpath = Path(outfile).resolve()
    maxid = 0

    if not dirpath.is_dir():
        raise click.exceptions.ClickException(f"{dirpath} is not a directory")
    for child in dirpath.iterdir():
        if child.suffix == '.obj':
            vertices = formats.OBJ.parse_vertices(child)
            log.debug(f"Computing TIN centroids and Morton-key")
            center = utils.mean_coordinate(vertices)
            morton_key = utils.morton_code(*center)
            tin_paths[morton_key] = Path(child)
    del vertices
    if len(tin_paths) > 1:
        click.echo(f"Found {len(tin_paths)} .obj files in {directory}")
    else:
        raise click.exceptions.ClickException(
            f"You need at least two .obj files in "
            f"{dirpath} in order to merge them")

    # Compute the tilesize. Tiles must be rectangular.
    # TODO: Compute the tilesize from the data itself. See utils for the draft
    tilesize = ctx.obj['cfg']['approximate_tilesize']
    assert isinstance(tilesize[0], float)
    tin_paths = dict(utils.compute_8neighbors(tin_paths, tilesize))

    # write centroids for testing
    log.debug("Writing centroids")
    ctr_path = dirpath / 'centroids.csv'
    with ctr_path.open('w') as cout:
        for i,morton_key in enumerate(sorted(tin_paths)):
            path = tin_paths[morton_key][0]
            x,y = utils.rev_morton_code(morton_key)
            cout.write(f"{i}\t{path.name}\t{x}\t{y}\n")

    # Loop through the TINs in Morton-order
    morton_order = sorted(tin_paths)
    processed_tiles = set() # Will store {pickle_path, ...}
    base = formats.factory.create('objmem')
    candidate = formats.factory.create('objmem')
    for i,candidate_key in enumerate(morton_order[1:], start=1):
        # In case of the first two tiles, the candidate is the second tile
        candidate_path = tin_paths[candidate_key][0]
        log.debug(f"Candidate={candidate_path.name}")

        if candidate_path.name == '37fz2_1.obj':
            a=1
            pass

        # Need to keep track which point belongs to which tile
        tile_lookup = {}  # Will store { point ID : tile name }
        # Collect all the base tiles that have been merged
        if len(processed_tiles) > 0:
            base_paths = processed_tiles.intersection(tin_paths[candidate_key][1])
            log.debug(f"Base={','.join(p.name for p in base_paths)}")
            # FIXME: Why is the Star .points and .stars are not initialized with a dict, but a None?
            base = formats.factory.create('objmem', points={}, stars={})
            for tile in base_paths:
                with tile.open(mode='rb') as fin:
                    try:
                        obj = pickle.load(fin, encoding='bytes')
                        base.points.update(obj.points)
                        base.stars.update(obj.stars)
                        tile_lookup.update((pid, obj.name) for pid in base.points)
                    except pickle.UnpicklingError as e:
                        log.error(
                            f"Cannot load pickle {tile}. Error:\n{e}")
        elif i == 1:
            # The base is the very first tile
            base_key = morton_order[0]
            base_path = tin_paths[base_key][0]
            log.debug(f"Base={base_path.name}")
            base.read(base_path)
            # FIXME: points need to be dicts
            _ = {v:point for v,point in enumerate(base.points)}
            base.points = _
            # D_base_valid = base.is_valid()
            # if not D_base_valid:
            #     log.warning(f"Base is not valid")
        else:
            raise click.exceptions.ClickException("Not finding the base tiles")

        candidate.read(candidate_path)
        # FIXME: points need to be dicts
        _ = {v: point for v, point in enumerate(candidate.points)}
        candidate.points = _

        # $ check validity
        # D_candidate_valid = candidate.is_valid()
        # if not D_candidate_valid:
        #     log.warning(f"Candidate is not valid")

        if in_memory:
            # FIXME: hack because points need to be dicks
            _m = max(base.stars)
            if maxid < _m:
                start_id = _m
                maxid = _m
            else:
                start_id = maxid
            # Reindex the candidate so that the point indices begin after
            # the base
            candidate.reindex(start_id=start_id)
            base.merge(candidate, strategy='deduplicate', precision=3)
        else:
            _m = max(base.stars)
            if maxid < _m:
                start_id = _m
                maxid = _m
            else:
                start_id = maxid
            # Reindex the candidate so that the point indices begin after
            # the base
            candidate.reindex(start_id=start_id)
            # Remove the co-located points from the candidate, merge the
            # stars of the co-located points and keep these stars only in the
            # base
            base.deduplicate(candidate=candidate, precision=3)

            # Pickle 'em for later...
            # would be cleaner to swap the if-else but i want the first tile
            # to be written first to the file

            if i == 1:
                pickle_path = base_path.with_suffix('.pickle')
                base.pickle_star(pickle_path)
                processed_tiles.update([pickle_path, ])

                pickle_path = candidate_path.with_suffix('.pickle')
                candidate.pickle_star(pickle_path)
                processed_tiles.update([pickle_path, ])
            else:
                pickle_path = candidate_path.with_suffix('.pickle')
                candidate.pickle_star(pickle_path)
                processed_tiles.update([pickle_path, ])

                # Need to split the base to tiles again, because it was updated
                # by the deduplicate
                tile_class = {
                    tile.stem: formats.OBJMem(
                    points={},
                    stars={},
                    name=tile.stem
                ) for tile in base_paths}
                for pid,tile in tile_lookup.items():
                    # $$$
                    if tile is None:
                        log.error("tile is None")
                    # $$$
                    if pid is not None:
                        tile_class[tile].points[pid] = base.points[pid]
                        tile_class[tile].stars[pid] = base.stars[pid]
                    else:
                        log.error("pid is None")
                for tile,star_obj in tile_class.items():
                    pickle_path = (dirpath / star_obj.name).with_suffix('.pickle')
                    if len(star_obj.points) != 0:
                        star_obj.pickle_star(pickle_path)
                del tile_class, tile_lookup, star_obj, obj
            # TODO: also need to update the base pickles somehow

    if in_memory:
        base.write(outpath)
    else:
        del base, candidate
        merged = formats.factory.create('objmem', points={}, stars={})
        for tile in processed_tiles:
            with tile.open('rb') as fin:
                tin = pickle.load(fin, encoding='bytes')
                merged.points.update(tin.points)
                merged.stars.update(tin.stars)
        merged.is_valid()
        merged.write(outpath, precision=3)
        for tile in processed_tiles:
            # Delete the pickle
            tile.unlink()


@click.command('import_many')
@click.argument('directory', type=click.Path(exists=True))
@click.pass_context
def import_many_cmd(ctx, directory, outfile):
    """Import TINs in OBJ format to PostgreSQL

    DIRECTORY is the folder that contains the TINs to be merged. All .obj files
    in this folder are merged together into a single file.
    """
    log = ctx.obj['log']
    tin_paths = {} # Will store { morton code : (file path, (range)) }
    dirpath = Path(directory).resolve()
    maxid = 0

    if not dirpath.is_dir():
        raise click.exceptions.ClickException(f"{dirpath} is not a directory")
    for child in dirpath.iterdir():
        if child.suffix == '.obj':
            vertices = formats.OBJ.parse_vertices(child)
            log.debug(f"Computing TIN centroids and Morton-key")
            center = utils.mean_coordinate(vertices)
            morton_key = utils.morton_code(*center)
            tin_paths[morton_key] = Path(child)
    del vertices
    if len(tin_paths) > 1:
        click.echo(f"Found {len(tin_paths)} .obj files in {directory}")
    else:
        raise click.exceptions.ClickException(
            f"You need at least two .obj files in "
            f"{dirpath} in order to merge them")

    # Compute the tilesize. Tiles must be rectangular.
    # TODO: Compute the tilesize from the data itself. See utils for the draft
    tilesize = ctx.obj['cfg']['approximate_tilesize']
    assert isinstance(tilesize[0], float)
    tin_paths = dict(utils.compute_8neighbors(tin_paths, tilesize))

    # write centroids for testing
    log.debug("Writing centroids")
    ctr_path = dirpath / 'centroids.csv'
    with ctr_path.open('w') as cout:
        for i,morton_key in enumerate(sorted(tin_paths)):
            path = tin_paths[morton_key][0]
            x,y = utils.rev_morton_code(morton_key)
            cout.write(f"{i}\t{path.name}\t{x}\t{y}\n")

    # Loop through the TINs in Morton-order
    morton_order = sorted(tin_paths)
    processed_tiles = set() # Will store {pickle_path, ...}
    base = formats.factory.create('objmem')
    candidate = formats.factory.create('objmem')
    for i,candidate_key in enumerate(morton_order[1:], start=1):
        # In case of the first two tiles, the candidate is the second tile
        candidate_path = tin_paths[candidate_key][0]
        log.debug(f"Candidate={candidate_path.name}")

        if candidate_path.name == '37fz2_1.obj':
            a=1
            pass

        # Need to keep track which point belongs to which tile
        tile_lookup = {}  # Will store { point ID : tile name }
        # Collect all the base tiles that have been merged
        if len(processed_tiles) > 0:
            base_paths = processed_tiles.intersection(tin_paths[candidate_key][1])
            log.debug(f"Base={','.join(p.name for p in base_paths)}")
            # FIXME: Why is the Star .points and .stars are not initialized with a dict, but a None?
            base = formats.factory.create('objmem', points={}, stars={})
            for tile in base_paths:
                with tile.open(mode='rb') as fin:
                    try:
                        obj = pickle.load(fin, encoding='bytes')
                        base.points.update(obj.points)
                        base.stars.update(obj.stars)
                        tile_lookup.update((pid, obj.name) for pid in base.points)
                    except pickle.UnpicklingError as e:
                        log.error(
                            f"Cannot load pickle {tile}. Error:\n{e}")
        elif i == 1:
            # The base is the very first tile
            base_key = morton_order[0]
            base_path = tin_paths[base_key][0]
            log.debug(f"Base={base_path.name}")
            base.read(base_path)
            # FIXME: points need to be dicts
            _ = {v:point for v,point in enumerate(base.points)}
            base.points = _
            # D_base_valid = base.is_valid()
            # if not D_base_valid:
            #     log.warning(f"Base is not valid")
        else:
            raise click.exceptions.ClickException("Not finding the base tiles")

        candidate.read(candidate_path)
        # FIXME: points need to be dicts
        _ = {v: point for v, point in enumerate(candidate.points)}
        candidate.points = _

        # $ check validity
        # D_candidate_valid = candidate.is_valid()
        # if not D_candidate_valid:
        #     log.warning(f"Candidate is not valid")

        _m = max(base.stars)
        if maxid < _m:
            start_id = _m
            maxid = _m
        else:
            start_id = maxid
        # Reindex the candidate so that the point indices begin after
        # the base
        candidate.reindex(start_id=start_id)
        # Remove the co-located points from the candidate, merge the
        # stars of the co-located points and keep these stars only in the
        # base
        base.deduplicate(candidate=candidate, precision=3)




main.add_command(import_cmd)
main.add_command(export_cmd)
main.add_command(delete_cmd)
main.add_command(subset_cmd)
main.add_command(merge_cmd)