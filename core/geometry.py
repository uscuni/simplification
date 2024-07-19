"""Geometry-related helpers"""

from itertools import combinations

import geopandas as gpd
import numpy as np
import shapely
from libpysal import graph
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
    # if line.within(poly):
    #     return True

    # intersection = line.intersection(poly)
    # return abs(intersection.length - line.length) <= rtol

    within = shapely.within(line, poly)
    if within.all():
        return within

    intersection = shapely.intersection(line, poly)
    return np.abs(shapely.length(intersection) - shapely.length(line)) <= rtol


def voronoi_skeleton(lines, poly=None, snap_to=None, distance=2, buffer=None):
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
    if buffer is None:
        buffer = distance * 20
    if not poly:
        poly = shapely.box(*lines.total_bounds)
    # get an additional line around the lines to avoid infinity issues with Voronoi
    extended_lines = list(lines) + [poly.buffer(buffer).exterior]

    # interpolate lines to represent them as points for Voronoi
    shapely_lines = extended_lines
    points, ids = shapely.get_coordinates(
        shapely.segmentize(shapely_lines, distance / 2), return_index=True
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
    for a, b in combinations(range(len(lines)), 2):
        mask = (
            np.isin(mapped[:, 0], [a, b])
            & np.isin(mapped[:, 1], [a, b])
            & (mapped[:, 0] != mapped[:, 1])
        )
        verts = rigde_vertices[mask]

        # generate the line in between the lines
        edgeline = shapely.line_merge(
            shapely.multilinestrings(voronoi_diagram.vertices[verts])
        )

        if not edgeline.within(poly):
            edgeline = shapely.intersection(edgeline, poly.buffer(-distance))
        # edgelines.append(edgeline)
        intersection = shapely_lines[b].intersection(shapely_lines[a])
        if not intersection.is_empty:
            # edgeline = shapely.snap(edgeline, intersection, distance * 10)
            edgeline = shapely.union(
                edgeline, shapely.shortest_line(edgeline.boundary, intersection)
            )
            # to_add.append(shapely.shortest_line(edgeline.boundary, intersection))
        edgelines.append(edgeline)

    if snap_to is None:
        to_add.append(
            shapely.shortest_line(shapely.union_all(edgelines).boundary, poly.boundary)
        )
    else:
        union = shapely.union_all(snap_to)
        edgelines_df = gpd.GeoDataFrame(geometry=edgelines)
        comp_labels = graph.Graph.build_contiguity(
            edgelines_df[~edgelines_df.is_empty], rook=False
        ).component_labels
        comp_counts = comp_labels.value_counts()
        components = edgelines_df.dissolve(comp_labels)
        if len(components) > 1:
            for comp_label, comp in components.geometry.items():
                if not comp.intersects(poly.boundary) or comp_counts[comp_label] == 1:
                    to_add.append(shapely.shortest_line(union, comp.boundary))
        else:
            for target in snap_to:
                to_add.append(shapely.shortest_line(target, components.boundary.item()))

    edgelines = np.concatenate([edgelines, to_add])
    edgelines = shapely.simplify(edgelines, distance / 2)

    return edgelines


# connections_to_fix = new_connnections[~connections_within]
# new_connnections = new_connnections[connections_within]
# fixes = []
# for g in connections_to_fix:
#     affected_nodes = relevant_nodes[relevant_nodes.intersects(g)]
#     affected_edge = edges[edges.covers(affected_nodes.union_all())].geometry.item()
#     chull_segments = shapely.get_parts(
#         affected_edge.convex_hull.boundary.intersection(artifact.geometry)
#     )
#     fixes.append(
#         shapely.line_merge(
#             shapely.union_all(
#                 chull_segments[shapely.intersects(chull_segments, artifact.geometry)]
#             )
#         )
#     )
# new_connnections = np.concatenate([new_connnections, fixes])
