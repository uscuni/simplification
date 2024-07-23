import pathlib

import geopandas
import matplotlib

__all__ = [
    "plot_case",
]


def plot_case(
    lines: geopandas.GeoDataFrame,
    vertices: geopandas.GeoDataFrame,
    title: str,
    fpath: str | pathlib.Path,
    **savefig_kwargs: dict,
) -> None:
    """Plot a single Case-Type example of our simplification protocol.
    For example, 'Case 01 - Parallel Roads - Original.'

    Parameters
    ----------
    lines : geopandas.GeoDataFrame
        Line data to plot. See ``core.protocol.generate_case()``.
    vertices : geopandas.GeoDataFrame
        Vertex data to plot. See ``core.protocol.generate_case()``.
    title : str
        Fully formatted title of the resultant plot.
    fpath : str | pathlib.Path
        Full path and file name to save out resultant plot as a ``.png``.
    **savefig_kwargs : dict
        Keyword arguments passed into ``matplotlib.pyplot.savefig()``.
    """

    # Ensure directory exists
    if isinstance(fpath, str):
        fpath = pathlib.Path(fpath)
    fpath.parent.mkdir(parents=True, exist_ok=True)

    ax = lines.plot(figsize=(5, 5), color="black", linewidth=0.5, linestyle="-")
    ax.set_axis_off()
    vertices.plot(ax=ax, color="k", markersize=2)
    matplotlib.pyplot.title(title)
    matplotlib.pyplot.savefig(fpath, **savefig_kwargs)
