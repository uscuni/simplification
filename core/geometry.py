"""Geometry-related helpers"""

from itertools import combinations

import numpy as np
import shapely
from scipy import spatial

__all__ = ["is_within", "voronoi_skeleton"]


def is_within(
    line: shapely.LineString, poly: shapely.Polygon, rtol: float = -1e4
) -> bool:
    """Check if the line is within a polygon with a set relative tolerance.

    Parameters
    ----------
    line : shapely.LineString
        Input line to check relationship.
    poly : shapely.Polygon
        Input polygon to check relationship.
    rtol : float (default -1e4)
        The set relative tolerance.

    Returns
    -------
    bool
        ``True`` if ``line`` is either entirely within ``poly`` or if
        ``line`` is within `poly`` based on a relaxed ``rtol`` relative tolerance.
    """
    if line.within(poly):
        return True

    intersection = line.intersection(poly)
    return abs(intersection.length - line.length) <= rtol


def voronoi_skeleton(lines, poly=None, snap_to=None, distance=2):
    """
    Returns average geometry.


    Parameters
    ----------
    lines : array_like
        LineStrings connected at endpoints
    poly : shapely.geometry.Polygon
        polygon enclosed by `lines`
    distance : float
        distance for interpolation

    Returns list of averaged geometries
    """
    if not poly:
        poly = shapely.box(*lines.total_bounds)
    # get an additional line around the lines to avoid infinity issues with Voronoi
    extended_lines = list(lines) + [poly.buffer(distance).exterior]

    # interpolate lines to represent them as points for Voronoi
    shapely_lines = extended_lines
    points, ids = shapely.get_coordinates(
        shapely.segmentize(shapely_lines, distance), return_index=True
    )

    # remove duplicated coordinates
    unq, count = np.unique(points, axis=0, return_counts=True)
    mask = np.isin(points, unq[count > 1]).all(axis=1)
    points = points[~mask]
    ids = ids[~mask]

    # generate Voronoi diagram
    voronoi_diagram = spatial.Voronoi(points)

    # get all rigdes and filter only those between the two lines
    pts = voronoi_diagram.ridge_points
    mapped = np.take(ids, pts)
    rigde_vertices = np.array(voronoi_diagram.ridge_vertices)

    # iterate over segment-pairs and keep rigdes between input geometries
    edgelines = []
    for a, b in combinations(range(len(lines)), 2):
        mask = (
            np.isin(mapped[:, 0], [a, b])
            & np.isin(mapped[:, 1], [a, b])
            & (mapped[:, 0] != mapped[:, 1])
        )
        verts = rigde_vertices[mask]

        # generate the line in between the lines
        edgeline = shapely.simplify(
            shapely.line_merge(
                shapely.multilinestrings(voronoi_diagram.vertices[verts])
            ),
            distance / 2,
        )
        intersection = shapely_lines[b].intersection(shapely_lines[a])
        if intersection.is_empty:
            edgelines.append(edgeline)
            if snap_to is None:
                edgelines.append(shapely.shortest_line(edgeline, poly.boundary))
            else:
                for target in snap_to:
                    edgelines.append(shapely.shortest_line(edgeline, target))
        else:
            snapped = shapely.snap(
                edgeline, shapely_lines[b].intersection(shapely_lines[a]), distance
            )
            edgelines.append(snapped)
    return edgelines
