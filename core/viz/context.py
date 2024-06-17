import pathlib

import contextily as cx
import cv2
import geopandas
import matplotlib.pyplot as plt
import pyproj
import shapely

__all__ = [
    "path",
    "location",
    "param_plot",
    "video",
]


def path(myclass: str, fpath_pack: pathlib.Path) -> pathlib.Path:
    """make subfolder for plot saving"""
    fpath_class = fpath_pack / myclass
    fpath_class.mkdir(parents=True, exist_ok=True)
    return fpath_class


def location(
    mypoint: tuple[float, float], crs: str | int | pyproj.CRS, buffer: int = 250
) -> geopandas.GeoSeries:
    """get center frame for clipping"""
    return (
        geopandas.GeoDataFrame(geometry=[shapely.Point(mypoint)], crs="epsg:4326")
        .to_crs(crs)
        .buffer(buffer, cap_style=3)
    )


def param_plot(
    n: geopandas.GeoDataFrame,
    e: geopandas.GeoDataFrame,
    loc: geopandas.GeoSeries,
    fmt_title: str,
    fmt_fname: str,
    fpath: pathlib.Path,
    crs: str | int | pyproj.CRS,
    close: bool = True,  # use False for image testing
):
    """make 1 context-param plot."""
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
    if close:
        plt.close()


def video(fpath: pathlib.Path):
    """make video context-param set"""
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
