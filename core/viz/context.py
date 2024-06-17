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


def path(context: str, fpath_pack: pathlib.Path) -> pathlib.Path:
    """Make subfolder for plot saving.

    Parameters
    ----------
    context : str
        Name for the context of the plot. For example, ``'intersection'``.
    fpath_pack: pathlib.Path
        Directory and name for the package being tested. For example, ``'osmnx'``.

    Returns
    -------
    fpath_context : pathlib.Path
        ``fpath_pack`` with ``context`` as an appended subdirectory.
    """
    fpath_context = fpath_pack / context
    fpath_context.mkdir(parents=True, exist_ok=True)
    return fpath_context


def location(
    point: tuple[float, float], crs: str | int | pyproj.CRS, buffer: int = 250
) -> geopandas.GeoSeries:
    """Get center frame for clipping.

    Parameters
    ----------
    point : tuple[float, float]
        Center point for zoomed in AOI.
    crs : str | int | pyproj.CRS
        Coordinate system.
    buffer : int (default 250)
        Desired buffer around ``point``.

    Returns
    -------
    geopandas.GeoSeries
        Buffered area around ``point``.
    """
    return (
        geopandas.GeoDataFrame(geometry=[shapely.Point(point)], crs="epsg:4326")
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
    close: bool = True,
) -> None:
    """Make 1 context-param plot.

    Parameters
    ----------
    n : geopandas.GeoDataFrame
        Nodes data.
    e : geopandas.GeoDataFrame
        Edges data.
    loc : geopandas.GeoSeries
        Zoomed in AOI. See ``context.location()``.
    fmt_title : str
        Plot title.
    fmt_fname : str
        File name to save
    fpath : pathlib.Path
        Directory for to save ``fmt_fname``.
    crs : str | int | pyproj.CRS
        Coordinate system.
    close : bool (default True)
        Close the plot -- **Use ``False`` for image testing**.
    """
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


def video(fpath: pathlib.Path) -> None:
    """Make video context-param set.

    Parameters
    ----------
    fpath : pathlib.Path
        Directory to save video.
    """
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
