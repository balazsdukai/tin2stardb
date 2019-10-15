# -*- coding: utf-8 -*-

"""Main class for the Star TIN structure."""

import logging
import re
from typing import List, Tuple, Mapping

from psycopg2 import extras, errors
import fiona
from psycopg2 import sql
from shapely import geos
from shapely.geometry import shape, Polygon


log = logging.getLogger(__name__)


class Star(object):
    """Main class for operating on a Star-TIN structure in memory."""

    def __init__(self, points: List[Tuple[float]] = None,
                 stars: Mapping[int, Tuple[int]] = None):
        self.stars = stars
        self.points = points

    def triangles(self):
        """Generate triangles from the stars.

        :returns: A generator over the triangles
        """
        for vtx,link in self.stars:
            # To avoid duplicate triangles:
            # If the vertex ID is smaller than the current neighbor in the link,
            # then skip the triangle and proceed to the next. The skipped
            # triangle will be attached to another star.
            for pt,pid in enumerate(link):
                if (vtx > link[pt-1]) or (vtx > link[pt]):
                    pass
                else:
                    yield (vtx, link[pt-1], link[pt])


class StarDb(object):
    """Main class for operating on a Star-TIN structure in a database."""

    def __init__(self, conn, schema):
        self.conn = conn
        self.schema = schema

    def insert(self, path, epsg, bbox):
        """To be overwritten by format-specific inserters."""
        pass

    def export(self, output):
        """To be overwritten by format-specific exporters."""
        pass

    def _max_vid(self):
        """Get the largest vertex ID from the TIN table."""
        query_params = {
            'id': self.schema.field.id.sqlid,
            'tin': self.schema.schema + self.schema.table
        }
        query = sql.SQL("""
        SELECT max({id}) FROM {tin}
        """).format(**query_params)
        log.debug(self.conn.print_query(query))
        resultset = self.conn.get_query(query)
        if resultset[0][0] is None:
            return 0
        else:
            return resultset[0][0]

    def _insert_generator(self, vertices, stars, epsg):
        """Prepare data for table insert.

        Updates the numbering of the vertex IDs to start from the last IDs in the
        TIN table.
        """
        log.info("Generating insert VALUES")
        x, y, z = [0, 1, 2]
        maxid = self._max_vid()
        for vid, star in stars:
            yield maxid + vid, f"SRID={epsg};POINT({vertices[vid][x]} {vertices[vid][y]} {vertices[vid][z]})", [
                maxid + v for v in star]

    def _insert_query(self, insert_generator):
        log.info("Inserting values...")
        query_params = {
            'tin': self.schema.schema + self.schema.table
        }
        try:
            with self.conn.conn as c:
                with c.cursor() as cur:
                    extras.execute_values(
                        cur,
                        sql.SQL("INSERT INTO {tin} VALUES %s").format(
                            **query_params),
                        insert_generator
                    )
        except BaseException as e:
            raise e

    def create_index(self):
        """Create the necessary indices on the TIN table."""
        query_params = {
            'id': self.schema.field.id.sqlid,
            'tin': self.schema.schema + self.schema.table
        }
        query = sql.SQL("""
        ALTER TABLE IF EXISTS {tin} ADD PRIMARY KEY ({id})
        """).format(**query_params)
        log.info(f"Creating primary key on {self.schema.table.string}")
        log.debug(self.conn.print_query(query))
        try:
            self.conn.send_query(query)
        except errors.InvalidTableDefinition as e:
            if 'multiple primary keys' in str(e):
                log.info(
                    f"Primary key already exists on {self.schema.table.string}")
            else:
                raise e

    def triangles(self):
        query_params = {
            'id': self.schema.field.id.sqlid,
            'star': self.schema.field.star.sqlid,
            'tin': self.schema.schema + self.schema.table
        }
        query = sql.SQL("""
        SELECT {id},{star} FROM {tin}
        """).format(**query_params)
        for vtx,star in self.conn.get_query(query):
            # To avoid duplicate triangles:
            # If the vertex ID is smaller than the current neighbor in the star,
            # then skip the triangle and proceed to the next. The skipped
            # triangle will be attached to another star.
            for pt,pid in enumerate(star):
                if (vtx > star[pt-1]) or (vtx > star[pt]):
                    pass
                else:
                    yield (vtx, star[pt-1], star[pt])

    def points(self):
        query_params = {
            'id': self.schema.field.id.sqlid,
            'geom': self.schema.field.geom.sqlid,
            'tin': self.schema.schema + self.schema.table
        }
        query = sql.SQL("""
        SELECT st_x({geom}), st_y({geom}), st_z({geom}) 
        FROM {tin}
        ORDER BY {id}
        """).format(**query_params)
        yield self.conn.get_query(query)

    @staticmethod
    def read_extent(extent: str) -> Tuple[Polygon, str]:
        """Reads a polygon from a file and returns it as Shapely polygon and
        EWKB.

        .. todo:: Not Implemented, just copied code

        :param extent: Path to a file (eg GeoJSON), contiaining a single polygon
        :return: A tuple of (Shapely Polygon, EWKB). If the extent doesn't have
            a CRS, then a WKB is returned instead of the EWKB.
        """
        raise NotImplementedError
        pattern = re.compile(r'(?<=epsg:)\d+', re.IGNORECASE)
        # Get clip polygon and set SRID
        with fiona.open(extent, 'r') as src:
            epsg = pattern.search(src.crs['init']).group()
            poly = shape(src[0]['geometry'])
            if epsg is None:
                log.warning(f"Did not find the EPSG code of the CRS: {src.crs}")
                wkb = poly.wkb_hex
                return poly, wkb
            else:
                # Change a the default mode to add this, if SRID is set
                geos.WKBWriter.defaults['include_srid'] = True
                # set SRID for polygon
                geos.lgeos.GEOSSetSRID(poly._geom, int(epsg))
                ewkb = poly.wkb_hex
                return poly, ewkb

    def _within_extent(self, ewkb: str) -> sql.Composed:
        """Return a query for the triangles that the extent. If
        features are not provided (as `feature_schema`), then select the tiles
        that intersect the extent.

        .. todo:: Not Implemented, just copied code
        """
        raise NotImplementedError
        query_params = {
            'tin': self.schema.schema + self.schema.table,
            'id': self.schema.field.id.sqlid,
            'geom': self.schema.table + self.schema.field.geom,
            'ewkb': sql.Literal(ewkb)
        }
        query = sql.SQL("""
        SELECT *
        FROM {tin}
        WHERE st_within({geom}, {ewkb}::geometry)
        """).format(**query_params)
        return query