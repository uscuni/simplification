from collections import Counter

import contextily as cx
import geopandas
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import networkx
import tobler

##################
# h3 funcs for gridcell-wise evaluation
##################


def make_grid(fua, res, proj_crs):
    """
    create a hex grid for given FUA (using original OSM dataset) and resolution.
    fua: int, ID of city (any of use cases: 809, 869, 1133, 1656, 4617, 4881, 8989)
    res: int, hex grid resolution
    (suggested: 8 - with cell side length ~530m or 9 - cell side length ~200m)
    """

    # CRS for h3 must be 4326
    orig_crs = "EPSG:4326"

    # get city geom and name
    meta = read_sample_data()  # noqa: F821 - now ``core.utils.read_sample_data()``
    geom = meta.loc[meta.eFUA_ID == fua, "geometry"]

    # read in OSM data
    orig = read_parquet_roads(fua)  # noqa: F821 - now ``core.utils.read_parquet_roads()``
    orig = orig.to_crs(orig_crs)

    assert meta.crs == orig.crs

    grid = tobler.util.h3fy(
        source=geom, resolution=res, clip=False, buffer=False, return_geoms=True
    )

    # have hex id as column rather than index
    grid["hex_id"] = grid.index
    grid = grid.reset_index(drop=True)

    # keep only the grid cells that contain some piece of the orig data
    mytree = orig.sindex
    q = mytree.query(grid.geometry, predicate="intersects", sort=True)
    grid = grid.loc[list(set(q[0]))]
    grid = grid.reset_index(drop=True)
    grid = grid.to_crs(proj_crs)
    return grid


def add_node_degree(gdf, graph):
    gdf["degree"] = gdf.apply(lambda x: networkx.degree(graph, x.nodeID), axis=1)
    return gdf


def get_geom_stats(gdf, geom):
    """
    clip gdf to geom (hex cell)
    return1: geometry count within the cell
    return2: total length of geometries within the cell
    """
    cell = geopandas.clip(gdf, geom)
    geom_count = len(cell)
    geom_length = cell.length.sum()
    return geom_count, geom_length


def _avg_degree(histdict):
    if len(histdict) > 0:
        return sum([k * v for k, v in histdict.items()]) / sum(list(histdict.values()))
    return 0


def get_node_stats(node_gdf, geom):
    cell = geopandas.clip(node_gdf, geom)
    nodecount = len(cell)
    degree_distr = dict(Counter(cell["degree"]))
    avg_degree = _avg_degree(degree_distr)
    return nodecount, degree_distr, avg_degree


def read_manual(fua, proj_crs):
    gdf = geopandas.read_parquet(f"../data/{fua}/manual/{fua}.parquet")
    gdf = gdf[["geometry"]]
    gdf = gdf.explode(ignore_index=True)
    gdf = gdf.to_crs(proj_crs)
    return gdf


def read_parenx(fua, option, proj_crs):
    gdf = geopandas.read_parquet(f"../data/{fua}/parenx/{option}.parquet")
    gdf = gdf.explode(ingore_index=True)
    gdf = gdf.to_crs(proj_crs)
    return gdf


# single-cell plot
def plot_cell(grid_id, grid, orig, base, comp):
    geom = grid.loc[grid_id, "geometry"]

    fig = plt.figure(figsize=(10, 14))
    spec = gridspec.GridSpec(3, 3, figure=fig)

    ax1 = fig.add_subplot(spec[0, 0])

    ax2 = fig.add_subplot(spec[0, 1])
    ax2.sharex(ax1)
    ax2.sharey(ax1)

    ax3 = fig.add_subplot(spec[0, 2])
    ax3.sharex(ax1)
    ax3.sharey(ax1)

    ax4 = fig.add_subplot(spec[1:, 0:4])
    ax4.sharex(ax1)
    ax4.sharey(ax1)

    # small plot: original network
    ax = ax1
    geopandas.clip(orig, geom).plot(ax=ax, color="black")
    ax.set_title("orig (osm)")

    ax = ax2
    geopandas.clip(base, geom).plot(ax=ax, color="red")
    ax.set_title("base")

    ax = ax3
    geopandas.clip(comp, geom).plot(ax=ax, color="blue")
    ax.set_title("comp")

    ax = ax4
    geopandas.clip(base, geom).plot(
        ax=ax, color="red", alpha=0.5, lw=2, zorder=1, label="base (manual)"
    )
    geopandas.clip(comp, geom).plot(
        ax=ax, color="blue", alpha=1, lw=2, zorder=2, label="simplified (auto)"
    )
    geopandas.clip(grid, geom).plot(
        ax=ax, facecolor=(0, 0, 0, 0), lw=0.5, edgecolor=(0, 0, 0, 1), zorder=3
    )
    cx.add_basemap(ax=ax, crs=grid.crs, source=cx.providers.CartoDB.Voyager, zorder=0)
    ax.set_title("Overlay")
    ax.legend()

    axs = [ax1, ax2, ax3, ax4]
    for ax in axs:
        ax.set_axis_off()

    return fig, ax
