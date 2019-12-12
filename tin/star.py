# -*- coding: utf-8 -*-

"""Main class for the Star TIN structure."""

from __future__ import annotations
import logging
import re
import random
from typing import Tuple, Mapping, Optional, Sequence, List
from copy import deepcopy
from statistics import mean, variance
from pathlib import Path
import pickle

from psycopg2 import extras, errors
import fiona
from psycopg2 import sql
from shapely import geos
from shapely.geometry import shape, Polygon

from tin import utils

log = logging.getLogger(__name__)


class Star(object):
    """Main class for operating on a Star-TIN structure in memory."""

    def __init__(self, points: Sequence[Tuple[float, float, float]] = None,
                 stars: Mapping[int, Sequence[int, ...]] = None,
                 points_dict = None, name=None):
        self.stars = stars
        self.points = points
        self.points_dict = points_dict
        self.name = name

    def triangles(self):
        """Generate triangles from the stars.

        :return: A generator over the triangles
        """
        for star, link in self.stars.items():
            # To avoid duplicate triangles:
            # If the vertex ID is smaller than the current neighbor in the link,
            # then skip the triangle and proceed to the next. The skipped
            # triangle will be attached to another star.
            for pt,pid in enumerate(link):
                if pt is None or pid is None:
                    log.error(f"Encountered a None in {star}:{link}")
                    break
                if pt:
                    if (star > link[pt-1]) or (star > link[pt]):
                        pass
                    else:
                        yield (star, link[pt-1], link[pt])
                else:
                    pass
                    # FIXME:   File "/home/balazs/Development/tin/tin/main.py", line 387, in merge_cmd
                    #     merged.write(outpath, precision=3)
                    #   File "/home/balazs/Development/tin/tin/formats.py", line 163, in write
                    #     for tri in self.triangles():
                    #   File "/home/balazs/Development/tin/tin/star.py", line 48, in triangles
                    #     if (star > link[pt-1]) or (star > link[pt]):
                    # TypeError: '>' not supported between instances of 'int' and 'NoneType'




    def pointlocation(self, point: Tuple[float, float]) -> Tuple[int, int, int]:
        """Return the triangle in which the given point is located

        Reference: pgTIN from Hugo Ledoux

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
                      dest: Tuple[float, float]) -> Sequence[Tuple[int, int, int], ...]:
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

    def add(self, candidate: Star) -> None:
        """Adds a TIN into the current one by combining their Stars.

        .. warning:: This operation modifies the TIN.

        .. note:: This operation does not remove duplicate points in case the
            two TINs overlap.
        """
        log.info("Adding a neighboring TIN to the Stars")
        maxid = max(self.stars)
        candidate_start_id = maxid+1
        self.points.update(candidate.points)
        stars_nbr = ((candidate_start_id+star, [candidate_start_id+v for v in link])
                     for star, link in candidate.stars.items())
        self.stars.update(stars_nbr)

    def add_seam(self, start_id: int, seam: Tuple) -> None:
        """Adds the seam to the Star.

        See :py:meth:`~.Star.get_seam` for an explanation on the *seam*.
        """
        seam_stars,seam_points = seam
        new_start_id = start_id+1
        stars_new = {new_start_id+star: [new_start_id+v for v in link]
                     for star,link in self.stars.items()}
        self.stars = stars_new
        # add the stars from the seam
        self.stars.update(seam_stars)
        # add the points from the seam
        self.points_dict = deepcopy(seam_points)
        points_new = ((new_start_id+i, vtx)
                      for i,vtx in enumerate(self.points))
        self.points_dict.update(points_new)
        del self.points

    def reindex(self, start_id: int) -> None:
        """Updates the vertex indices in the stars with a given staring ID."""
        start_id += 1
        stars_new = ((start_id + star, [start_id + v for v in link])
                     for star,link in self.stars.items())
        self.stars = dict(stars_new)
        # FIXME: points should be a dict, not a list
        points_new = ((start_id + i, vtx)
                      for i, vtx in self.points.items())
        self.points = dict(points_new)

    def find_duplicate_points(self, candidate: Star, precision: int = 3) -> List[int]:
        """Finds the points that are duplicated in the candidate Star.

        This function compares the points by casting their coordinates into
        strings using the given precision. Comparison is done based on the
        xy-coordinates only (not z).

        :param candidate: A candidate Star structure which is compared to *this*
            Star.
        :param precision: The number of decimals to use for comparing
            coordinates.
        :returns: A list of point IDs from *this* Star, that have a co-located
            pair in the candidate Star.
        """
        log.info(f"Precision is set to {precision} decimal digits")
        duplicates = []
        pt_hash_tbl_candidate = {}

        for vtx,pt in enumerate(candidate.points):
            pt_str = f"{pt[0]:.{precision}f},{pt[1]:.{precision}f}"
            pt_hash_tbl_candidate[pt_str] = vtx

        for vtx,pt in enumerate(self.points):
            pt_str = f"{pt[0]:.{precision}f},{pt[1]:.{precision}f}"
            if pt_str in pt_hash_tbl_candidate:
                # co-located point found, so we store the vertex ID
                duplicates.append(vtx)

        return duplicates

    def get_seam(self, vertices: List[int]) -> Tuple[Mapping, Mapping]:
        """Get a subset of the TIN by passing in a list of vertex IDs.

        The *seam* is the set of points along which two TINs are touching
        (`vertices`). These points exist in both TINs. This function returns
        the star and coordinates of `vertices`, *without* reindexing the vertex
        IDs.

        The *seam* is used in the streaming merge strategy, in which we need to
        keep the seam in memory, while removing the rest of the base TIN.

        .. note:: The function returns invalid stars, since the link of the
            stars points to vertices that not present in the returned pointlist,
            but only in the original Star. Also, the vertex IDs of the returned
            stars refer to the ordering of the original Star's pointlist and
            not to the returned pointlist.

        :returns: A subset of the Star without reindexing the vertex IDs. The
            subset is returned as a tuple of stars (dict) and coordinates
            (vertex ID : coordinates).
        """
        subset_stars = {star:link for star,link in self.stars.items()
                        if star in vertices}
        subset_vtx = {vtx:self.points[vtx] for vtx in subset_stars.keys()}
        return subset_stars, subset_vtx

    def deduplicate(self, candidate: Star = None, precision: int = 3, return_quality: bool = False) -> Tuple:
        """Remove duplicate points from a TIN.

        The general method for removing duplicates:

            1) Cast vertices into strings using the given precision
            2) Create hash-table (dict) and keep track of the vertex indices,
                not loosing the duplicates
            3) If two points are co-located, take their mean z coordinate as
                the new z
            4) Replace the duplicate with the original vertex

        If `candidate` is not provided, the function modifies the Star in place
        by combining the links of the stars where the star-points are
        co-located. This usually happens when you are merging two adjacent TINs
        into a single TIN.

        If `candidate` is provided, the function modifies both the Star
        instance on which this method is called (the *base*) and the Star
        instance that has co-located points with the base (the *candidate*).
        The principle for this operation is to remove the duplicated points
        from the candidate and add them to the base, while updating the stars
        in both base and candidate.

        Two stars are combined by adding the link from the candidate to the link
        of the base, thus the star of the base becomes the only
        star for co-located points. Then the combined link is sorted in CCW, since we still have
        all the vertices in memory and all the vertex indices are valid in
        both the base and candidate. Finally the co-located points are removed
        from the candidate, the vertex indices are updated, the stars are removed
        and the candidate links are updated to include the stars from the candidate.

        .. note:: Ignores the z coordinate for determining point co-location.

        .. warning:: The function requires that there are no common vertex
            indices in `base` and `candidate`, thus it is necessary to reindex
            the `candidate` beforehand.

        .. warning:: This operation modifies the TIN.

        :return: Tuple of None if return_quality is False; Else a tuple of:
            1) The list of co-located points on the boundary of the two TINs.
                The returned points are tuples of (x,y,z). The z coordinate is
                the average of the z coordinates of the co-located points.
            2) A tuple of (Minimum, Maximum, Mean, Variance) of the absolute
                differences of the z coordinates of the co-located points.
        """
        # base : the Star instance on which this method is called
        # candidate : the Star instance that shares some co-located points with
        # the base (adjacent to the base)
        log.info(f"Removing duplicate points from the TIN. "
                 f"Precision is set to {precision} decimal digits")
        pt_hash_tbl_base = {}
        z_differences = []
        common_points = []

        if candidate is None:
            for vtx2,pt in self.points.items():
                pt_str = f"{pt[0]:.{precision}f},{pt[1]:.{precision}f}"
                # if co-located points found
                if pt_str not in pt_hash_tbl_base:
                    pt_hash_tbl_base[pt_str] = vtx2
                else:
                    vtx1 = pt_hash_tbl_base[pt_str]
                    _z = self.points[vtx1][2]
                    z_differences.append(abs(pt[2] - _z))
                    # average z in case of co-located points
                    new_z = (pt[2] + _z) / 2
                    # keep the point that was found first and replace it's z-coordinate
                    new_pt = (self.points[vtx1][0], self.points[vtx1][1], new_z)
                    common_points.append(new_pt)
                    self.points[vtx1] = new_pt
                    del _z, new_pt, new_z
                    # go through the stars and replace the references of vtx2 to
                    # vtx1 in the links
                    for star in self.stars[vtx2]:
                        try:
                            i = self.stars[star].index(vtx2)
                            self.stars[star][i] = vtx1
                        except ValueError:
                            log.error(f"Expected to find vertex {vtx2} in the stars "
                                      f"of {vtx2}:{self.stars[vtx2]}")
                    # add the link of vtx2 to the star of vtx1
                    link_vtx1 = self.stars[vtx1]
                    link_vtx1.extend(v for v in self.stars[vtx2]
                                     if v not in link_vtx1)
                    # sort the link ccw
                    _d = dict(utils.sort_ccw(self.points, {vtx1: link_vtx1}))
                    # update the link of vtx1 with the merged link of vtx1+vtx2
                    self.stars[vtx1] = _d[vtx1]
                    # delete vtx2 star and point
                    self.points[vtx2] = None
                    del self.stars[vtx2]
                    del _d
            del pt_hash_tbl_base
            # Create a new points list, excluding the deleted points and create a
            # mapping of the old-new vertices
            vid = 0
            new_points = []
            old_new_map = {v:None for v in range(len(self.points))}
            for vtx, point in self.points.items():
                if point is None:
                    old_new_map[vtx] = None
                else:
                    old_new_map[vtx] = vid
                    new_points.append(point)
                    vid += 1
            # Reindex the stars with the new vertices
            new_stars = {}
            for old_v, new_v in old_new_map.items():
                if new_v is not None:
                    old_link = self.stars[old_v]
                    new_stars[new_v] = [old_new_map[v] for v in old_link
                                        if old_new_map[v] is not None]
            # Replace the points and stars
            self.points = deepcopy(new_points)
            self.stars = deepcopy(new_stars)
        else:
            # FIXME: this code below expects that Star.points is a dict, not a list
            # so for now we need to check and cast the Star.points
            if isinstance(self.points, list):
                newpts = ((v, pt) for v,pt in enumerate(self.points))
                self.points = dict(newpts)
                del newpts
            if not self.stars.keys().isdisjoint(candidate.stars.keys()):
                raise ValueError(f"Base {self.name} and candidate {candidate.name} cannot have common vertex indices.")
            pt_hash_tbl_candidate = {}
            updated_stars_base = []

            for vtx, pt in candidate.points.items():
                pt_str = f"{pt[0]:.{precision}f},{pt[1]:.{precision}f}"
                pt_hash_tbl_candidate[pt_str] = vtx
            del vtx, pt

            for star_base, pt in self.points.items():
                pt_str = f"{pt[0]:.{precision}f},{pt[1]:.{precision}f}"
                if pt_str not in pt_hash_tbl_candidate:
                    pt_hash_tbl_base[pt_str] = star_base
                else:
                    # Co-located point in base and candidate.
                    updated_stars_base.append(star_base)
                    # 1) Overwrite base point with the average z-coordinate.
                    star_candidate = pt_hash_tbl_candidate[pt_str]
                    _point = candidate.points[star_candidate]
                    if _point is None:
                        break
                    candidate_z = candidate.points[star_candidate][2]
                    z_differences.append(abs(pt[2] - candidate_z))
                    # Average z of the co-located points
                    new_z = (pt[2] + candidate_z) / 2
                    # Keep the point of the base and replace it's z-coordinate
                    new_pt = (self.points[star_base][0],
                              self.points[star_base][1],
                              new_z)
                    common_points.append(new_pt)
                    self.points[star_base] = new_pt
                    del candidate_z, new_pt, new_z

                    # 2) Take the link of the candidate and merge into the link
                    # of the base, while updating the references in the
                    # link of the candidate star to point to star_base.
                    #
                    # Replace star_candidate to star_base in each link in
                    # candidate that contains star_candidate
                    for star in candidate.stars[star_candidate]:
                        # At this point the link of the candidate already
                        # contains indices from the base
                        if star in candidate.stars:
                            try:
                                i = candidate.stars[star].index(star_candidate)
                                candidate.stars[star][i] = star_base
                            except ValueError:
                                log.error(f"Expected to find point {star_candidate} in "
                                          f"the link of {star_candidate}:{candidate.stars[star_candidate]}")
                    # Append the link of star_candidate to the link of star_base
                    link_base = self.stars[star_base]
                    link_base.extend(v for v in candidate.stars[star_candidate]
                                     if v not in link_base)
                    # Combine the points from the base and candidate into one
                    # container. But at this point we already merged the links,
                    # so we need to check if the link contains indices from the
                    # other Star.
                    combined_points = {v: self.points[v]
                                       for v in self.stars[star_base]
                                       if v in self.points}
                    combined_points[star_base] = self.points[star_base]
                    for v in candidate.stars[star_candidate]:
                        if v in candidate.points:
                            combined_points[v] = candidate.points[v]
                    combined_points[star_candidate] = candidate.points[star_candidate]
                    # FIXME: This is a quickfix for
                    #  7751d13f1d6171d6a79e22fe7eded06176e491d6 so that the
                    #  process runs through. Encountered a reference in a link
                    #  to a missing point (58233) in tile 37fz2_13 when merging
                    #  with 37fz2_12. The missing point is part of 37fz2_7,
                    #  but this tile is not in memory.
                    try:
                        # Sort the link ccw
                        _d = dict(utils.sort_ccw(combined_points, {star_base: link_base}))
                        # Update the link of star_base with the merged link of base+candidate
                        self.stars[star_base] = _d[star_base]
                        del _d
                    except KeyError:
                        log.error("Encountered missing point in sort_ccw()")
                        self.stars[star_base] = link_base
                    # Delete candidate star and point
                    candidate.points[star_candidate] = None
                    del candidate.stars[star_candidate]

            if len(updated_stars_base) == 0:
                log.warning(f"Did not find any co-located points in base {self.name} and candidate {candidate.name}")
                return None, None

            del star_base, pt, star_candidate, pt_hash_tbl_candidate, pt_hash_tbl_base
            # Create a new points list, excluding the deleted points and create a
            # mapping of the old-new vertices
            vid = max(self.points) + 1 # start at where base ends
            new_points = {}
            old_new_map = {v:None for v in candidate.points}
            old_new_map.update((v,v) for v in self.points)
            for vtx, point in candidate.points.items():
                if point is None:
                    old_new_map[vtx] = None
                else:
                    old_new_map[vtx] = vid
                    new_points[vid] = point
                    vid += 1
            del vtx, point
            # Reindex the candidate stars with the new vertices
            new_stars = {}
            for old_vtx, new_vtx in old_new_map.items():
                if new_vtx is not None and old_vtx in candidate.stars:
                    old_link = candidate.stars[old_vtx]
                    new_stars[new_vtx] = [old_new_map[v] for v in old_link
                                          if old_new_map[v] is not None]
            # Replace the points and stars
            candidate.points = deepcopy(new_points)
            candidate.stars = deepcopy(new_stars)
            del new_points, new_stars, old_vtx, new_vtx
            # Also need to update the links of the base that include points
            # from the candidate
            for star_base in updated_stars_base:
                for vtx in self.stars[star_base]:
                    if vtx in old_new_map:
                        # Found a vertex from the candidate
                        try:
                            i = self.stars[star_base].index(vtx)
                            self.stars[star_base][i] = old_new_map[vtx]
                        except ValueError:
                            pass
        del old_new_map
        # TODO: need to update the global max_id

        if len(z_differences) < 2:
            return_quality = False
            log.warning("Only 1 co-located point were found. Cannot compute merge quality")
        else:
            # Compute the quality of the TIN merge
            min_z = round(min(z_differences), 3)
            max_z = round(max(z_differences), 3)
            mean_z = round(mean(z_differences), 3)
            variance_z = round(variance(z_differences), 3)
            log.info(
                f"Difference in z coordinates of co-located points in the two TINs: "
                f"nr_points={len(z_differences)} "
                f"min={min_z}, "
                f"max={max_z}, "
                f"mean={mean_z}, "
                f"variance={variance_z}")
        if return_quality:
            return common_points, (min_z, max_z, mean_z, variance_z)
        else:
            return None, None


    def streaming_deduplicate(self):
        """ """


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

    def merge(self, candidate: Star, strategy: str = 'deduplicate',
              precision: int = 3, return_quality: bool = False) -> Optional[None, Tuple[List[Tuple[float], ...]], Tuple[float]]:
        """Merge a TIN into the current one.

        :param candidate: An adjacent TIN in Star structure that needs to be
            merged into the current TIN.
        :param strategy: The strategy to use for merging. If *deduplicate*, then
            the two TINs are expected to touch, thus have a set of co-located
            points along the edge in which they are touching.
        :param return_quality: If True, return a list common (co-located) points
            (x,y,z) in the two TINs. The z coordinate of the coordinates is the
            mean of the z coordinate of the two co-located points.
        :return: None if return_quality is False; Else a tuple of:
            1) The list of co-located points on the boundary of the two TINs.
                The returned points are tuples of (x,y,z). The z coordinate is
                the average of the z coordinates of the co-located points.
            2) A tuple of (Minimum, Maximum, Mean, Variance) of the absolute
                differences of the z coordinates of the co-located points.
        """
        # side, segment = utils.find_side(self.points[1:],
        #                                 candidate.points[1:], abs_tol=0.1)
        if strategy.lower() == 'deduplicate':
            self.add(candidate)
            common_pts, quality = self.deduplicate(precision=precision,
                                          return_quality=return_quality)
            if return_quality:
                return common_pts, quality
            else:
                return None
        elif strategy.lower() == 'streaming_deduplicate':
            # Find duplicate points between the two TINs. The two TINs are
            # expected to touch along these points.
            duplicates = self.find_duplicate_points(candidate, precision=3)
            #
            subset_stars, subset_points = self.get_seam(duplicates)
        else:
            raise ValueError(f"Unknown merge strategy {strategy}")

    def is_valid(self) -> bool:
        """Validates a Star structure.

        (1) checks if the links are consistent, thus vertex A is also present in
            the link of vertex B, if vertex B is in the link of vertex A
        (2) checks if the link is ordered CounterClockWise
        (3) checks for each triangle that the vertices in adjacent triangles
            are consistent
        """
        # Check that there are no Null points
        null_points = [id for id, point in self.points.items()
                       if point is None]
        if len(null_points) > 0:
            log.warning(f"There are {len(null_points)} NULL points. NULL={null_points}")
        # Check if there is any Null in the links
        null_in_link_of_star = [star
                                for star,link in self.stars.items()
                                for v in link
                                if v is None]
        if len(null_in_link_of_star) > 0:
            log.warning(f"There are NULL in the link of {null_in_link_of_star}. Aborting validation.")
            return False
        # Check that links are consistent
        validation_summary = {star:[] for star in self.stars}
        for star, consistent in utils.link_is_consistent(self.stars):
            validation_summary[star].append(consistent)
            if not consistent:
                log.warning(f"Link of star {star} is not consistent. "
                            f"{self.stars[star]}; "
                            f"{[self.points[v] for v in self.stars[star]]}")
        # Check if link is CCW
        for star, ccw in utils.link_is_ccw(self.points, self.stars):
            validation_summary[star].append(ccw)
            if not ccw:
                log.warning(f"Link of star {star} is not CCW. "
                            f"{self.stars[star]}; "
                            f"{[self.points[v] for v in self.stars[star]]}")
        # Check that triangles are consistent
        _t = []
        for tri, consistent in utils.triangle_is_consistent(self.stars, self.triangles()):
            _t.append(consistent)
            if not consistent:
                log.warning(f"Triangle {tri} is not consistent")

        no_null = len(null_points) == 0
        no_null_in_link = len(null_in_link_of_star) == 0
        links_consistent = all(result[0] for star, result in validation_summary.items())
        ccw = all(result[1] for star, result in validation_summary.items())
        tris_consistent = all(_t)

        return no_null and \
               no_null_in_link and \
               links_consistent and \
               ccw and \
               tris_consistent

    def write_star(self, path: Path, mode='a'):
        """Write the Star to a file, conserving its structure.

        The output format is:
            s star_ID link_ID1 ...
            v point_ID x-coordinate y-coordinate z-coordinate

        Thus the format is similar to Wavefront OBJ, but in this case the
        rows beginning with `s` denote a *star*. In a *star* the `star_` and
        `link_` IDs are the 0-based vertex IDs. The rows beginning with `v`
        denote the vertex coordinates, just like in OBJ.
        """
        log.info(f"Writing Star to {path} in mode={mode}")
        with path.open(mode=mode) as fout:
            for star,link in self.stars.items():
                fout.write(f"s {star} {' '.join(str(i) for i in link)}\n")
            for id,vtx in self.points.items():
                fout.write(f"v {id} {' '.join(str(i) for i in vtx)}\n")

    def pickle_star(self, path: Path):
        """Pickle the Star and write it to a file."""
        log.info(f"Pickling Star to {path}")
        with path.open(mode='wb') as fout:
            try:
                pickle.dump(self, fout, protocol=pickle.HIGHEST_PROTOCOL)
            except pickle.PicklingError as e:
                log.error(f"Cannot pickle the Star instance. Error:\n{e}")


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