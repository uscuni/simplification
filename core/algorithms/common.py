import momepy


def continuity(roads, angle_threshold=120):
    """Assign COINS-based information to roads

    Parameters
    ----------
    roads : GeoDataFrame
        Road network

    Returns
    -------
    GeoDataFrame
        original roads with additional columns
    """
    roads = roads.copy()
    # Measure continuity of street network
    coins = momepy.COINS(roads, angle_threshold=angle_threshold, flow_mode=True)

    # Assing continuity group
    group, end = coins.stroke_attribute(True)
    roads["coins_group"] = group
    roads["coins_end"] = end

    # Assign length of each continuity group and a number of segments within the group.
    coins_grouped = roads.length.groupby(roads.coins_group)
    roads["coins_len"] = coins_grouped.sum()[roads.coins_group].values
    roads["coins_count"] = coins_grouped.size()[roads.coins_group].values

    return roads, coins
