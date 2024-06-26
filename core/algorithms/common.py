import geopandas
import momepy
from libpysal import graph


def continuity(roads: geopandas.GeoDataFrame) -> geopandas.GeoDataFrame:
    """Assign COINS-based information to roads.

    Parameters
    ----------
    roads :  geopandas.GeoDataFrame
        Road network.

    Returns
    -------
    roads : geopandas.GeoDataFrame
        The input ``roads`` with additional columns.
    """
    roads = roads.copy()
    # Measure continuity of street network
    coins = momepy.COINS(roads)

    # Assing continuity group
    roads["coins_group"] = coins.stroke_attribute()

    # Assign length of each continuity group and a number of segments within the group.
    coins_grouped = roads.length.groupby(roads.coins_group)
    roads["coins_len"] = coins_grouped.sum()[roads.coins_group].values
    roads["coins_count"] = coins_grouped.size()[roads.coins_group].values

    # Figure out which segments are on the ends of their continuity groups.
    roads_contiguity = graph.Graph.build_contiguity(
        roads, rook=False
    ).assign_self_weight()
    grouper = roads.coins_group.take(
        roads_contiguity._adjacency.index.codes[1]
    ).groupby(roads_contiguity._adjacency.index.codes[0])
    coins_end = []
    for i, g in grouper:
        coins_end.append((g == g.loc[i]).sum() in (1, 2))
    roads["coins_end"] = coins_end

    return roads
