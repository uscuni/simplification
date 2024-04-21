import json
import pathlib

import contextily as cx
import cv2
import geopandas
import matplotlib.pyplot as plt
import networkx
import pyproj
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
