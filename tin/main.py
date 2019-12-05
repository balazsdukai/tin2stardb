# -*- coding: utf-8 -*-

"""Main app."""

from os import listdir
from pathlib import Path
from sys import stdout
import logging
import json

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
            fmt = suffix.lower().split('.')[1]
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
@click.pass_context
def merge_cmd(ctx, directory, outfile):
    """Merge TINs in OBJ format.

    DIRECTORY is the folder that contains the TINs to be merged. All .obj files
    in this folder are merged together into a single file.

    OUTFILE is the name of the output file. The OUTFILE is written into
    DIRECTORY.
    """
    log = ctx.obj['log']
    tin_paths = {} # Will store (morton code : file path)
    dirpath = Path(directory).resolve()
    outpath = Path(outfile).resolve()

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

    # write centroids for testing
    log.debug("Writing centroids")
    ctr_path = dirpath / 'centroids.csv'
    with ctr_path.open('w') as cout:
        for i,morton_key in enumerate(sorted(tin_paths)):
            path = tin_paths[morton_key]
            x,y = utils.rev_morton_code(morton_key)
            cout.write(f"{i}\t{path.name}\t{x}\t{y}\n")

    # Loop through the TINs in Morton-order
    morton_order = sorted(tin_paths)
    for i in range(len(morton_order)):
        base_key = morton_order[i]
        try:
            candidate_key = morton_order[i+1]
        except IndexError:
            # Reached the last TIN
            candidate_key = None
        base_path = tin_paths[base_key]
        log.debug(f"Processing base={base_path.name}")
        base = formats.factory.create('objmem')
        base.read(base_path)
        if candidate_key:
            candidate_path = tin_paths[candidate_key]
            log.debug(f"Processing candidate={candidate_path.name}")
            candidate = formats.factory.create('objmem')
            candidate.read(candidate_path)
            # The two TINs (base, candidate) are expected to touch along some
            # path, which we call *seam*. Find the points of the seam.
            duplicates = base.find_duplicate_points(candidate)
            # Extract the stars of the seam from the base
            seam = base.get_seam(duplicates)
            # Add the seam to the candidate
            base_maxid = max(base.stars)
            # Remove the seam from the base

        # We only write out the base
        base.write_star(outpath, mode='a')


    # outpath = Path(outfile).resolve()
    # base = formats.factory.create('objmem')
    # candidate = formats.factory.create('objmem')
    # base.read(inpaths[0])
    # candidate.read(inpaths[1])
    # base.merge(candidate, strategy='deduplicate')
    # base.write(outpath)


main.add_command(import_cmd)
main.add_command(export_cmd)
main.add_command(delete_cmd)
main.add_command(subset_cmd)
main.add_command(merge_cmd)