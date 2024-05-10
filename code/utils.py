import json
import pathlib

import contextily as cx
import cv2
import geopandas
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import networkx
import pyproj
import tobler
from shapely.geometry import Point

__all__ = [
    "city_fua",
    "fua_city",
    "read_sample_data",
    "read_parquet_roads",
    "graph_size",
    "load_usecases",
    "viz_class_path",
    "viz_class_location",
    "viz_class_param_plot",
    "viz_class_video",
]


ddir = pathlib.Path("..", "data")
parq = "parquet"
samp_parq = ddir / f"sample.{parq}"

# dict of fua ID: cityname
fua_city = {
    1133: "Aleppo",
    869: "Auckland",
    4617: "Bucaramanga",
    809: "Douala",
    1656: "Liège",
    4881: "Salt Lake City",
    8989: "Wuhan",
}

# dict of cityname: fua ID
city_fua = {c: f for f, c in fua_city.items()}


def read_sample_data() -> geopandas.GeoDataFrame:
    return geopandas.read_parquet(samp_parq)


def read_parquet_roads(fua: int) -> geopandas.GeoDataFrame:
    """Read OSM roads from parquet format; return bare columns."""
    return geopandas.read_parquet(pathlib.Path(ddir, f"{fua}", f"roads_osm.{parq}"))[
        ["highway", "geometry"]
    ].reset_index(drop=True)


def graph_size(info: str, g: networkx.Graph) -> str:
    return f"{info}\n\t* {g}"


def load_usecases(city: str) -> tuple[dict, pathlib.Path]:
    """Load use cases for a city and also return the directory path."""
    path_base = pathlib.Path("..", "usecases", str(city_fua[city]))
    with open(path_base / "points.json") as f:
        points = json.load(f)
    return points, path_base


##################
# Viz funcs
##################


def viz_class_path(myclass: str, fpath_pack: pathlib.Path) -> pathlib.Path:
    """make subfolder for plot saving"""
    fpath_class = fpath_pack / myclass
    fpath_class.mkdir(parents=True, exist_ok=True)
    return fpath_class


def viz_class_location(
    mypoint: tuple[float, float], crs: pyproj.CRS, buffer: int = 250
) -> geopandas.GeoSeries:
    """get center frame for clipping"""
    return (
        geopandas.GeoDataFrame(geometry=[Point(mypoint)], crs="epsg:4326")
        .to_crs(crs)
        .buffer(buffer, cap_style=3)
    )


def viz_class_param_plot(
    n: geopandas.GeoDataFrame,
    e: geopandas.GeoDataFrame,
    loc: geopandas.GeoSeries,
    fmt_title: str,
    fmt_fname: str,
    fpath: pathlib.Path,
    crs: pyproj.CRS,
):
    """make 1 class-param plot."""
    # make a plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    # clip geometries to box
    n_clipped = n.copy().clip(loc)
    e_clipped = e.copy().clip(loc)

    # plot
    e_clipped.plot(ax=ax, zorder=1, color="black", linewidth=2)
    n_clipped.plot(ax=ax, zorder=2, color="red", markersize=15, alpha=0.9)
    cx.add_basemap(ax=ax, source=cx.providers.CartoDB.Voyager, crs=crs)
    ax.set_axis_off()
    ax.set_title(fmt_title)
    plt.tight_layout()

    # save to subfolder
    fig.savefig(fpath / f"{fmt_fname}.png", dpi=300)
    plt.close()


def viz_class_video(fpath: pathlib.Path):
    """make video class-param set"""
    fps = 1
    images = sorted(fpath.glob("*.png"))
    video_name = pathlib.Path(f"{fpath}.mp4")
    frame = cv2.imread(str(images[0]))
    height, width, layers = frame.shape
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")
    video = cv2.VideoWriter(str(video_name), fourcc, fps, (width, height))
    for image in images:
        video.write(cv2.resize(cv2.imread(str(image)), (width, height)))
    cv2.destroyAllWindows()
    video.release()


##################
# h3 funcs
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
    meta = read_sample_data()
    geom = meta.loc[meta.eFUA_ID == fua, "geometry"]

    # read in OSM data
    orig = read_parquet_roads(fua)
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
