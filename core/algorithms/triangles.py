import geopandas as gpd
import pandas as pd
import shapely
from tqdm.auto import tqdm

from ..geometry import is_within


def case_1(edges, *, geom, to_drop, to_add, nodes, distance_threshold=1.05):
    """If one edge is its own contiuity group while the other two are part of larger
    continuity blocks where one continues before and after the artifact

    - If the shortest line from the entry point to the street that continues before and
    after the artifact is singificantly shorter than the existing connection use that
    one; unless it is not fully within the original artifact - in that case use the
    existing connection

    - If the shortest line from the entry point to the street that continues before and
    after the artifact is not singificantly shorter than the existing connection, use
    the existing connection

    Parameters
    ----------
    edges : GeoDataFrame
        geometries forming the artifact
    geom : shapely.Polygon
        polygonal representation of the artifact
    to_drop : list
        list collecting geometries to be dropped
    to_add : list
        list collecting geometries to be added
    nodes : GeoSeries
        nodes forming the artifact
    distance_threshold : float
        use new line if distance_threshold * length is shorter than existing link
    """
    relevant_nodes = nodes.iloc[nodes.sindex.query(geom, predicate="intersects")]
    main_road = edges.geometry[~edges.coins_end].item()
    entry_node = relevant_nodes[relevant_nodes.disjoint(main_road)]
    shortest_line = shapely.shortest_line(entry_node.item(), main_road)

    to_drop.append((edges.coins_count == 1).idxmax())
    edges_fixed = edges[edges.coins_count != 1]
    existing_link = edges_fixed.geometry[edges_fixed.coins_end]

    if (
        shortest_line.length * distance_threshold < existing_link.length.item()
        and is_within(shortest_line, geom)
    ):
        to_drop.append(existing_link.index[0])
        to_add.append(shortest_line)


def case_2(edges, *, to_drop):
    """If one edge is its own contiuity group while the other two are part of larger
    continuity blocks where both continue before and after the artifact

    - drop the edge that is its own continuity group

    Parameters
    ----------
    edges : GeoDataFrame
        geometries forming the artifact
    to_drop : list
        list collecting geometries to be dropped
    """
    to_drop.append((edges.coins_count == 1).idxmax())


def case_5(edges, *, geom, to_drop, to_add, roads, nodes):
    """If all three edges are part of larger continuity blocks & only one of them is
    part of a continuity block that continues before and after the artifact:

    - if one of the non-primary continuity blocks intersects the primary
        - keep the existing node between the two continuity blocks, remove the remaining
          edge
    - if the non-primary continuity blocks do not intersect the primary
        - drop the two and link the entry point to the continuity block that continues
          before and after the artifact if the shortest line is fully wihtin the
          artifact. otherwise keep the shorter of the two links

    Parameters
    ----------
    edges : GeoDataFrame
        geometries forming the artifact
    geom : shapely.Polygon
        polygonal representation of the artifact
    to_drop : list
        list collecting geometries to be dropped
    to_add : list
        list collecting geometries to be added
    roads : GeoDataFrame
        original road network in full extent
    nodes : GeoSeries
        nodes forming the artifact
    """
    edges_extended = roads.iloc[roads.sindex.query(geom, predicate="intersects")]
    main_road = edges[~edges.coins_end]
    main_group = main_road.coins_group.item()
    non_main = edges_extended[edges_extended.coins_group != main_group]

    intersections_per_group = (
        non_main.intersects(main_road.geometry.item())
        .groupby(non_main.coins_group)
        .sum()
        > 1
    )
    intersections_per_group = intersections_per_group[
        intersections_per_group.index.isin(edges.coins_group)
    ]

    if intersections_per_group.any():
        group_to_drop = intersections_per_group.idxmin()
        to_drop.append(edges[edges.coins_group == group_to_drop].index[0])
    else:
        relevant_nodes = nodes.iloc[nodes.sindex.query(geom, predicate="intersects")]
        entry_node = relevant_nodes[relevant_nodes.disjoint(main_road.geometry.item())]
        shortest_line = shapely.shortest_line(
            entry_node.item(), main_road.geometry.item()
        )

        if is_within(shortest_line, geom):
            to_add.append(shortest_line)
            to_drop.extend(edges[edges.coins_group != main_group].index.to_list())
        else:
            to_drop.extend(edges[edges.coins_group != main_group].length.idxmin())


def case_6(edges, *, to_drop):
    """If all three edges are part of larger continuity blocks & two of them are part of
    a continuity block that continues before and after the artifact

    - drop the one that is not part of a continuity block that continues before and
      after he artifact

    Parameters
    ----------
    edges : GeoDataFrame
        geometries forming the artifact
    to_drop : list
        list collecting geometries to be dropped
    """
    to_drop.append(edges[edges.coins_end].index[0])


def case_7(edges, *, geom, to_drop, to_add, nodes):
    """If all three edges are part of larger continuity blocks & none of then continues
    before and after the artifact

    - consider the largest edge a main one and create a shortest path from the remaining
      entry point. Use the shortest path if it is fully within the artifact. Otherwise
      use the shorted of the two edges.

    Parameters
    ----------
    edges : GeoDataFrame
        geometries forming the artifact
    geom : shapely.Polygon
        polygonal representation of the artifact
    to_drop : list
        list collecting geometries to be dropped
    to_add : list
        list collecting geometries to be added
    nodes : GeoSeries
        nodes forming the artifact
    """
    longest = edges.length.idxmax()
    main_road = edges.loc[[longest]]

    relevant_nodes = nodes.iloc[nodes.sindex.query(geom, predicate="intersects")]
    entry_node = relevant_nodes[relevant_nodes.disjoint(main_road.geometry.item())]
    shortest_line = shapely.shortest_line(entry_node.item(), main_road.geometry.item())
    if is_within(shortest_line, geom):
        to_add.append(shortest_line)
        to_drop.extend(edges.index.drop([longest]).to_list())
    else:
        to_drop.extend([edges.length.idxmin()])


def case_8(edges, *, geom, to_drop, to_add, nodes):
    """If all three edges are part of the same continuiry block (roundabout):

    - drop all of them and link the entry points to the centroid of the triangle

    Parameters
    ----------
    edges : GeoDataFrame
        geometries forming the artifact
    geom : shapely.Polygon
        polygonal representation of the artifact
    to_drop : list
        list collecting geometries to be dropped
    to_add : list
        list collecting geometries to be added
    nodes : GeoSeries
        nodes forming the artifact
    """
    centroid = geom.centroid
    relevant_nodes = nodes.iloc[nodes.sindex.query(geom, predicate="intersects")]

    to_drop.extend(edges.index.to_list())
    to_add.extend(shapely.shortest_line(relevant_nodes, centroid).to_list())


def resolve(roads, artifacts):
    """Resolve triangular artifacts

    Parameters
    ----------
    roads : GeoDataFrame
        road networks
    artifacts : GeoDataFrame
        face artifacts

    Returns
    -------
    GeoSeries
    """
    first = shapely.get_point(roads.geometry, 0)
    last = shapely.get_point(roads.geometry, -1)
    possible_nodes = gpd.GeoSeries(
        pd.concat([first, last], ignore_index=True), crs=roads.crs
    )
    nodes = possible_nodes[~possible_nodes.duplicated()]

    to_drop = []
    to_add = []

    for geom in tqdm(artifacts.geometry, total=len(artifacts)):
        edges = roads.iloc[roads.sindex.query(geom, predicate="covers")]

        if len(edges) != 3:  # ensure triangle
            continue

        if edges.coins_end.sum() == 2 and (edges.coins_count == 1).sum() == 1:
            case_1(
                edges,
                geom=geom,
                to_drop=to_drop,
                to_add=to_add,
                nodes=nodes,
                distance_threshold=1.05,
            )

        elif edges.coins_end.sum() == 1 and (edges.coins_count == 1).sum() == 1:
            case_2(edges, to_drop=to_drop)

        elif edges.coins_end.sum() == 2 and (edges.coins_count == 1).sum() == 0:
            case_5(
                edges,
                geom=geom,
                to_drop=to_drop,
                to_add=to_add,
                roads=roads,
                nodes=nodes,
            )

        elif edges.coins_end.sum() == 1 and (edges.coins_count == 1).sum() == 0:
            case_6(edges, to_drop=to_drop)
        elif edges.coins_end.sum() == 3 and (edges.coins_count == 1).sum() == 0:
            case_7(edges, geom=geom, to_drop=to_drop, to_add=to_add, nodes=nodes)
        elif (
            edges.coins_end.sum() == 0
            and (edges.coins_count == 1).sum() == 0
            and edges.coins_group.nunique() == 1
        ):
            case_8(edges, geom=geom, to_drop=to_drop, to_add=to_add, nodes=nodes)

    return pd.concat(
        [roads.geometry.drop(to_drop), gpd.GeoSeries(to_add, crs=roads.crs)]
    )
