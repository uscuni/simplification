"""Geometry-related helpers"""

import warnings

import geopandas as gpd
import numpy as np
import shapely
from libpysal import graph
from scipy import spatial

from .algorithms.nodes import consolidate_nodes

__all__ = ["is_within", "voronoi_skeleton"]


def is_within(
    line: np.ndarray[shapely.Geometry], poly: shapely.Polygon, rtol: float = 1e-4
) -> bool:
    """Check if the line is within a polygon with a set relative tolerance.

    Parameters
    ----------
    line : np.ndarray[shapely.LineString]
        Input line to check relationship.
    poly : shapely.Polygon
        Input polygon to check relationship.
    rtol : float (default -1e4)
        The set relative tolerance.

    Returns
    -------
    np.ndarray[bool]
        ``True`` if ``line`` is either entirely within ``poly`` or if
        ``line`` is within `poly`` based on a relaxed ``rtol`` relative tolerance.
    """
    # if line.within(poly):
    #     return True

    # intersection = line.intersection(poly)
    # return abs(intersection.length - line.length) <= rtol

    within = shapely.within(line, poly)
    if within.all():
        return within

    intersection = shapely.intersection(line, poly)
    return np.abs(shapely.length(intersection) - shapely.length(line)) <= rtol


def voronoi_skeleton(
    lines,
    poly=None,
    snap_to=None,
    max_segment_length=1,
    buffer=None,
    secondary_snap_to=None,
    limit_distance=2,
    consolidation_tolerance=None,
):
    """
    Returns average geometry.


    Parameters
    ----------
    lines : array_like
        LineStrings connected at endpoints
    poly : shapely.geometry.Polygon
        polygon enclosed by `lines`
    snap_to : gpd.GeoSeries
        series of geometries that shall be connected to the skeleton
    distance : float
        distance for interpolation
    buffer : float
        optional custom buffer distance for dealing with Voronoi infinity issues
    consolidation_tolerance : float
        tolerance passed to node consolidation within the resulting skeleton. If None,
        no consolidation happens

    Returns array of averaged geometries
    """
    if buffer is None:
        buffer = max_segment_length * 20
    if not poly:
        poly = shapely.box(*lines.total_bounds)
    # get an additional line around the lines to avoid infinity issues with Voronoi
    extended_lines = list(lines) + [poly.buffer(buffer).boundary]

    # interpolate lines to represent them as points for Voronoi
    shapely_lines = extended_lines
    points, ids = shapely.get_coordinates(
        shapely.segmentize(shapely_lines, max_segment_length), return_index=True
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
    to_add = []
    splitters = []

    # determine the negative buffer distance to avoid overclipping of narrow polygons
    # this can still result in some missing links, but only in rare cases
    dist = min(
        [limit_distance, shapely.ops.polylabel(poly).distance(poly.boundary) * 0.4]
    )
    limit = poly.buffer(-dist)

    # drop ridges that are between points coming from the same line
    selfs = mapped[:, 0] == mapped[:, 1]
    buff = (mapped == mapped.max()).any(axis=1)
    mapped = mapped[~(selfs | buff)]
    rigde_vertices = rigde_vertices[~(selfs | buff)]
    unique = np.unique(np.sort(mapped, axis=1), axis=0)

    for a, b in unique:
        mask = ((mapped[:, 0] == a) | (mapped[:, 0] == b)) & (
            (mapped[:, 1] == a) | (mapped[:, 1] == b)
        )

        verts = rigde_vertices[mask]

        # generate the line in between the lines
        edgeline = shapely.line_merge(
            shapely.multilinestrings(voronoi_diagram.vertices[verts])
        )

        # check if the edgeline is within polygon
        if not edgeline.within(limit):
            # if not, clip it by the polygon with a small negative buffer to keep
            # the gap between edgeline and poly boundary to avoid possible
            # overlapping lines
            edgeline = shapely.intersection(edgeline, limit)

            # in edge cases, this can result in a MultiLineString with one sliver part
            if edgeline.geom_type == "MultiLineString":
                parts = shapely.get_parts(edgeline)
                edgeline = parts[np.argmax(shapely.length(parts))]

        # check if a, b lines share a node
        intersection = shapely_lines[b].intersection(shapely_lines[a])
        # if they do, add shortest line from the edgeline to the shared node and
        # combine it with the edgeline. Also, avoid an inner loop in more complex input
        # that would create connection across
        if not intersection.is_empty and not (
            intersection.geom_type == "MultiPoint"
            and (len(intersection.geoms) == 2 and len(lines) != 2)
        ):
            # we need union of edgeline and shortest because snap is buggy in GEOS
            # and line_merge as well. This results in a MultiLineString but we can
            # deal with those later. For now, we just need this extended edgeline to
            # be a single geometry to ensure the component discovery below works as
            # intended
            # get_parts is needed as in case of voronoi based on two lines, these
            # intersect on both ends, hence both need to be extended
            edgeline = shapely.union(
                edgeline,
                shapely.union_all(
                    shapely.shortest_line(
                        shapely.get_parts(intersection), edgeline.boundary
                    )
                ),
            )
        # add final edgeline to the list
        edgelines.append(edgeline)

    edgelines = np.array(edgelines)[~(shapely.is_empty(edgelines))]

    if edgelines.shape[0] > 0:
        # if there is no explicit snapping target, snap to the boundary of the polygon
        # via the shortest line. That is by definition always within the polygon
        # (Martin thinks)
        if snap_to is not False:
            if snap_to is None:
                sl = shapely.shortest_line(
                    shapely.union_all(edgelines).boundary, poly.boundary
                )
                to_add.append(sl)
                splitters.append(shapely.get_point(sl, -1))

            # if we have some snapping targets, we need to figure out
            # what shall be snapped to what
            else:
                additions, splits = snap_to_targets(
                    edgelines, poly, snap_to, secondary_snap_to
                )
                to_add.extend(additions)
                splitters.extend(splits)

            # concatenate edgelines and their additions snapping to edge
            edgelines = np.concatenate([edgelines, to_add])
        # simplify to avoid unnecessary point density and some wobbliness
        edgelines = shapely.simplify(edgelines, max_segment_length)
    # drop empty
    edgelines = edgelines[edgelines != None]  # noqa: E711

    edgelines = shapely.line_merge(edgelines[shapely.length(edgelines) > 0])
    if np.unique(shapely.get_type_id(edgelines)).shape[0] > 1:
        edgelines = shapely.get_parts(edgelines)

    if consolidation_tolerance and edgelines.shape[0] > 0:
        edgelines = consolidate_nodes(
            edgelines, tolerance=consolidation_tolerance, preserve_ends=True
        ).geometry.to_numpy()

    return edgelines, splitters


def snap_to_targets(edgelines, poly, snap_to, secondary_snap_to=None):
    to_add = []
    to_split = []
    # cast edgelines to gdf
    edgelines_df = gpd.GeoDataFrame(geometry=edgelines)
    # build queen contiguity on edgelines and extract component labels
    comp_labels = graph.Graph.build_contiguity(
        edgelines_df[~(edgelines_df.is_empty | edgelines_df.geometry.isna())],
        rook=False,
    ).component_labels
    # compute size of each component
    comp_counts = comp_labels.value_counts()
    # get MultiLineString geometry per connected component
    components = edgelines_df.dissolve(comp_labels)

    # if there are muliple components, loop over all and treat each
    if len(components) > 1:
        for comp_label, comp in components.geometry.items():
            # if component does not interest the boundary, it needs to be snapped
            # if it does but has only one part, this part interesect only on one
            # side (the node remaining from the removed edge) and needs to be
            # snapped on the other side as well
            if (
                (not comp.intersects(poly.boundary))
                or comp_counts[comp_label] == 1
                or (
                    not comp.intersects(shapely.union_all(snap_to))
                )  # ! this fixes one thing but may break others
            ):
                # add segment composed of the shortest line to the nearest snapping
                # target. We use boundary to snap to endpoints of edgelines only
                sl = shapely.shortest_line(comp.boundary, shapely.union_all(snap_to))
                if is_within(sl, poly):
                    to_add.append(sl)
                    to_split.append(shapely.get_point(sl, -1))
                else:
                    if secondary_snap_to is not None:
                        sl = shapely.shortest_line(
                            comp.boundary, shapely.union_all(secondary_snap_to)
                        )
                        to_split.append(shapely.get_point(sl, -1))
                        to_add.append(sl)
    else:
        # if there is a single component, ensure it gets a shortest line to an
        # endpoint from each snapping target
        for target in snap_to:
            sl = shapely.shortest_line(components.boundary.item(), target)
            if is_within(sl, poly):
                to_split.append(shapely.get_point(sl, -1))
                to_add.append(sl)
            else:
                warnings.warn(
                    "Could not create a connection as it would lead outside "
                    "of the artifact.",
                    UserWarning,
                    stacklevel=2,
                )
    return to_add, to_split
