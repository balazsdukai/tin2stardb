# -*- coding: utf-8 -*-

"""Main class for the Star TIN structure."""

from __future__ import annotations
import logging
import re
import random
from typing import List, Tuple, Mapping

from psycopg2 import extras, errors
import fiona
from psycopg2 import sql
from shapely import geos
from shapely.geometry import shape, Polygon

from tin import utils

log = logging.getLogger(__name__)


class Star(object):
    """Main class for operating on a Star-TIN structure in memory."""

    def __init__(self, points: List[Tuple[float, float, float]] = None,
                 stars: Mapping[int, Tuple[int, ...]] = None):
        self.stars = stars
        self.points = points

    def triangles(self):
        """Generate triangles from the stars.

        :return: A generator over the triangles
        """
        for vtx,link in self.stars.items():
            # To avoid duplicate triangles:
            # If the vertex ID is smaller than the current neighbor in the link,
            # then skip the triangle and proceed to the next. The skipped
            # triangle will be attached to another star.
            for pt,pid in enumerate(link):
                if (vtx > link[pt-1]) or (vtx > link[pt]):
                    pass
                else:
                    yield (vtx, link[pt-1], link[pt])

    def add(self, neighbor: Star) -> None:
        """Adds a TIN into the current one by combining their Stars.

        .. warning:: This operation modifies the current TIN.

        .. note:: This operation does not remove duplicate points in case the
            two TINs overlap.
        """
        maxid = max(self.stars)
        self.points += neighbor.points[1:] # because OBJ has 1-based indexing so the first value is None
        stars_nbr = ((star+maxid, [v+maxid for v in link])
                     for star, link in neighbor.stars.items())
        self.stars.update(stars_nbr)

    def merge(self, neighbor: Star) -> None:
        """Merge a TIN into the current one."""
        side, segment = utils.find_side(self.points[1:],
                                        neighbor.points[1:], abs_tol=0.1)

    def pointlocation(self, point: Tuple[float, float]) -> Tuple[int, int, int]:
        """Return the triangle in which the given point is located

        :param point: The point in question, given as (x,y) coordinates
        :return: The triangle in which the point is located
        """
        # -- find closest starting point
        n = 5
        mindist = 1e99
        closeid = 1
        for i in range(n):
            rid = random.randint(1, max(self.stars))
            temp = self.points[rid]
            tempdist = utils.distance(point, temp)
            if tempdist < mindist:
                mindist = tempdist
                closeid = rid
        # print('closest point is:', closeid)
        # -- find first triangle
        rv = closeid
        lk = self.stars[closeid]
        for i in range(len(lk)):
            rv1 = lk[i]
            if utils.orientation(self.points[rv], self.points[rv1], point) == -1:
                # -- find previous in star of closeid
                ii = lk.index(rv1)
                if ii == 0:
                    other = lk[len(lk) - 1]
                else:
                    other = lk[ii - 1]
                tr = (closeid, other, rv1)
                break
        # print('first triangle found:', tr[0], tr[1], tr[2])
        while 1:
            v0 = tr[0]
            v1 = tr[1]
            v2 = tr[2]
            if utils.orientation(self.points[v0], self.points[v1], point) != -1:
                if utils.orientation(self.points[v1], self.points[v2], point) != -1:
                    found = (v0, v1, v2)
                    break
                else:
                    # -- jump-jump-jump!
                    lk = self.stars[v1]
                    ii = lk.index(v2)
                    if ii == 0:
                        other = lk[len(lk) - 1]
                    else:
                        other = lk[ii - 1]
                    tr = (v1, other, v2)  # -- [1, prev, 2]
            else:
                # -- jump-jump-jump!
                lk = self.stars[v0]
                ii = lk.index(v1)
                if ii == 0:
                    other = lk[len(lk) - 1]
                else:
                    other = lk[ii - 1]
                tr = (v0, other, v1)
            # print('tr', tr[0], tr[1], tr[2])
        return tr

    def straight_walk(self, start: Tuple[float, float],
                      dest: Tuple[float, float]) -> List[Tuple[int, int, int], ...]:
        """Return a list of triangles along a straight line through the TIN.

        :param start: First point of the line, given as (x,y) coordinates
        :param dest: Second point of the line, given as (x,y) coordinates
        :return: List of a triangles along the line
        """
        start_tri = self.pointlocation(start)
        trail = []
        trail.append(start_tri)
        tr = start_tri
        #print('first triangle found:', tr[0], tr[1], tr[2])
        while 1:
            v0 = tr[0]
            v1 = tr[1]
            v2 = tr[2]
            if utils.orientation(self.points[v0], self.points[v1], dest) != -1:
                if utils.orientation(self.points[v1], self.points[v2], dest) != -1:
                    found = (v0, v1, v2)
                    break
                else:
                    #-- jump-jump-jump!
                    lk = self.stars[v1]
                    ii = lk.index(v2)
                    if ii == 0:
                        other = lk[len(lk)-1]
                    else:
                        other = lk[ii-1]
                    if utils.orientation(start, dest, self.points[other]) != -1:
                        tr = (v1, other, v2) #-- [1, prev, 2]
                    else:
                        tr = (other, v2, v1)
                    trail.append(tr)
            else:
                #-- jump-jump-jump!
                lk = self.stars[v0]
                ii = lk.index(v1)
                if ii == 0:
                    other = lk[len(lk)-1]
                else:
                    other = lk[ii-1]
                if utils.orientation(start, dest, self.points[other]) != -1:
                    tr = (v0, other, v1)
                else:
                    tr = (other, v1, v0)
                trail.append(tr)
            #print('tr', tr[0], tr[1], tr[2])
        return trail

    def sew(self):
        """.. todo:: Sew two TINs into a topologically valid TIN by traversing along a
                straight line that is the where the two TINs touch.

        The algorithm is an adaptation of the 'straight walk' algorithm. As
        input it requires a *single TIN* and a *line segment*. The single TIN is
        is a combination of two TINs by combining their points into a single
        array and reindexing the vertices of the second. The line segment is
        the line along which the two TINs touch.

        The algorithm progresses along the *line segment* and removes the
        duplicate vertices, by replacing the vertex from the second TIN with
        one from the first TIN in the triangle stars, thereby creating a
        topologically connected TIN.
        """

    def deduplicate(self):
        """Remove duplicate points from a TIN.
        1) Combine the two TINs into a single Star
        2) Cast vertices into strings using the given precision, or just keep
            the precision as it is
        3) Create hash-table (dict) and keep track of the vertex indices,
            not loosing the duplicates
        4) Replace the duplicate with the original vertex
        :return:
        """



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
        # Internally we use a 0-based index, but Postgres uses a 1-based primary key
        maxid = self._max_vid() + 1
        for vid, star in stars:
            yield maxid + vid, \
                  f"SRID={epsg};POINT({vertices[vid][x]} {vertices[vid][y]} {vertices[vid][z]})",\
                  [maxid + v for v in star]

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