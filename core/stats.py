from collections import Counter

import geopandas
import networkx
import shapely

__all__ = [
    "add_node_degree",
    "get_edge_stats",
    "get_node_stats",
]


def add_node_degree(
    node_gdf: geopandas.GeoDataFrame, graph: networkx.Graph
) -> geopandas.GeoDataFrame:
    """Add a node degree column to a nodes dataframe.

    Parameters
    ----------
    node_gdf : geopandas.GeoDataFrame
        Node data.
    graph : networkx.Graph
        Graph representation of a street network.

    Returns
    -------
    node_gdf : geopandas.GeoDataFrame
        Updated node data.
    """
    node_gdf["degree"] = node_gdf["nodeID"].map(lambda x: networkx.degree(graph, x))

    return node_gdf


def get_edge_stats(
    edge_gdf: geopandas.GeoDataFrame, geom: shapely.Polygon
) -> tuple[int, float]:
    """Calculate basic edge stats for a single H3 hex cell.

    Parameters
    ----------
    edge_gdf : geopandas.GeoDataFrame
        Edge data.
    geom: shapely.Polygon
        Single H3 hex cell.

    Returns
    -------
    edge_count : int
        Number of edges within the cell.
    edge_length : float
        Cumulative length of edges within the cell.
    """
    cell = geopandas.clip(edge_gdf, geom)
    edge_count = len(cell)
    edge_length = cell.length.sum()
    return edge_count, edge_length


def _avg_degree(histdict: dict) -> int | float:
    """Helper to calculate average node degree."""
    if len(histdict) > 0:
        return sum([k * v for k, v in histdict.items()]) / sum(list(histdict.values()))
    return 0


def get_node_stats(
    node_gdf: geopandas.GeoDataFrame, geom: shapely.Polygon
) -> tuple[int, dict, float]:
    """Calculate basic node stats for a single H3 hex cell.

    Parameters
    ----------
    node_gdf : geopandas.GeoDataFrame
        Node data.
    geom: shapely.Polygon
        Single H3 hex cell.

    Returns
    -------
    node_count : int
        Number of nodes within the cell.
    degree_distr : dict
        Distribution of node degrees within the cell.
    avg_degree : float
        Mean node degree within the cell.
    """
    cell = geopandas.clip(node_gdf, geom)
    node_count = len(cell)
    degree_distr = dict(Counter(cell["degree"]))
    avg_degree = _avg_degree(degree_distr)
    return node_count, degree_distr, avg_degree
