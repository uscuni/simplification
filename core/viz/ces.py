import geopandas as gpd
import pandas as pd
import shapely


def plot_ces(roads, geom, ax):
    edges = roads.iloc[roads.sindex.query(geom, predicate="covers")].copy()
    edges_touching = roads.iloc[
        roads.sindex.query(edges.union_all(), predicate="touches")
    ].copy()

    edges["ces"] = None
    edges.loc[edges.coins_count == 1, "ces"] = "S"
    ecg = edges.coins_group
    all_ends = edges[edges.coins_end]
    ae_cg = all_ends.coins_group
    edges.loc[~ecg.isin(ae_cg), "ces"] = "C"
    edges.loc[edges.ces.isna(), "ces"] = "E"

    ces = []
    for _, row in edges_touching.iterrows():
        if row.coins_group in edges.coins_group.values:
            ces.append(edges[edges.coins_group == row.coins_group].ces.iloc[0])
        else:
            ces.append("-")
    edges_touching["ces"] = ces

    combined = pd.concat([edges, edges_touching])
    colors = {
        "C": "#11A579",
        "E": "#E73F74",
        "S": "#3969AC",
        "-": "#bbbbbb",
    }

    # Extract vertices
    _endpoint = lambda x, gdf=edges: shapely.get_point(  # noqa: E731
        gdf.geometry.explode(ignore_index=True), x
    )

    # _endpoint = get_endpoint(gdf)
    _endpoints = pd.concat([_endpoint(0), _endpoint(-1)]).drop_duplicates()
    vertices = gpd.GeoDataFrame(geometry=_endpoints, crs=edges.crs)

    for e in ["C", "E", "S", "-"]:
        sub = combined[combined.ces == e]
        if not sub.empty:
            sub.clip(geom.buffer(15)).plot(
                ax=ax,
                color=colors[e],
            )
    gpd.GeoSeries([geom]).plot(
        ax=ax, edgecolor="lightgray", hatch="..", facecolor="none"
    )
    vertices.plot(ax=ax, facecolor="lightgrey", edgecolor="k", markersize=25, zorder=2)

    # ensure we always plot a rectangle
    bounds = combined.clip(geom.buffer(15)).total_bounds
    center_x = (bounds[0] + bounds[2]) / 2
    center_y = (bounds[1] + bounds[3]) / 2
    center_point = shapely.Point(center_x, center_y)
    dist_x = abs(center_x - bounds[2])
    dist_y = abs(center_y - bounds[3])
    max_dist = max(dist_x, dist_y)
    buffer_geometry = center_point.buffer(max_dist, cap_style="square")
    gpd.GeoSeries([buffer_geometry], crs=combined.crs).plot(
        ax=ax,
        edgecolor="none",
        facecolor="none",
    )


def plot_simplified(roads, geom, ax):
    edges = roads.iloc[roads.sindex.query(geom, predicate="intersects")]
    edges.clip(geom.buffer(15)).plot(ax=ax, color="#F2B701")
    gpd.GeoSeries([geom]).plot(
        ax=ax, edgecolor="lightgray", hatch="..", facecolor="none", linewidth=0
    )

    # Extract vertices
    _endpoint = lambda x, gdf=edges: shapely.get_point(  # noqa: E731
        gdf.geometry.explode(ignore_index=True), x
    )

    # _endpoint = get_endpoint(gdf)
    _endpoints = pd.concat([_endpoint(0), _endpoint(-1)]).drop_duplicates()
    vertices = gpd.GeoDataFrame(geometry=_endpoints, crs=edges.crs).clip(
        geom.buffer(15)
    )
    vertices.plot(ax=ax, facecolor="lightgrey", edgecolor="k", markersize=25, zorder=2)
    # ensure we always plot a rectangle
    bounds = (
        pd.concat([gpd.GeoSeries([geom]), edges]).clip(geom.buffer(15)).total_bounds
    )
    center_x = (bounds[0] + bounds[2]) / 2
    center_y = (bounds[1] + bounds[3]) / 2
    center_point = shapely.Point(center_x, center_y)
    dist_x = abs(center_x - bounds[2])
    dist_y = abs(center_y - bounds[3])
    max_dist = max(dist_x, dist_y)
    buffer_geometry = center_point.buffer(max_dist, cap_style="square")
    gpd.GeoSeries([buffer_geometry], crs=edges.crs).plot(
        ax=ax,
        edgecolor="none",
        facecolor="none",
    )
