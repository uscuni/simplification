import contextily
import geopandas
import matplotlib

__all__ = [
    "plot_aoi",
    "plot_analysis",
    "plot_cell",
]


def plot_aoi(
    cell_gdf: geopandas.GeoDataFrame,
    city: str,
    res: int,
):
    """Plot the H3 cell grid over the AOI.

    Parameters
    ----------
    cell_gdf : geopandas.GeoDataFrame
        H3 hexbins covering the FUA that contain data.
    city : str
        City name.
    res : int
        H3 hex grid resolution.
    """

    ax = cell_gdf.boundary.plot(linewidth=0.5, alpha=1, color="black", figsize=(5, 8))
    ax.set_axis_off()
    ax.set_title(f"Grid cell analysis for {city}, H3 resolution {res}")
    contextily.add_basemap(
        ax=ax, crs=cell_gdf.crs, source=contextily.providers.CartoDB.Voyager
    )


def plot_analysis(
    cell_gdf: geopandas.GeoDataFrame, var: str, info: None | str = None
) -> None:
    """Plot the H3 cell grid symbolized by analysis values.

    Parameters
    ----------
    cell_gdf : geopandas.GeoDataFrame
        H3 hexbins covering the FUA that contain data.
    var : str
        Analysis values column name in ``cell_gdf``.
    info : None | str (default None)
        Addition information for the title.
    """

    cmap = matplotlib.cm.coolwarm
    ax = cell_gdf.plot(
        column=var,
        norm=matplotlib.colors.CenteredNorm(vcenter=1),
        cmap=cmap,
        legend=True,
        figsize=(10, 10),
    )
    ax.set_axis_off()
    ax.set_title(f"{var} {info}" if info else f"{var}")


def plot_cell(
    grid_id: int | str,
    grid: geopandas.GeoDataFrame,
    orig: geopandas.GeoDataFrame,
    orig_name: str,
    base: geopandas.GeoDataFrame,
    base_name: str,
    comp: geopandas.GeoDataFrame,
    comp_name: str,
) -> tuple[matplotlib.figure.Figure, dict]:
    """Generate a mosaic plot comparing original OSM streets to a baseline
    method and secondary method.

    Parameters
    ----------
    grid_id : int
        If ``int`` index location of grid cell.
        if ``str`` the H3 hex ID of grid cell.
    grid : geopandas.GeoDataFrame
        H3 hexbins covering the FUA that contain data.
    orig : geopandas.GeoDataFrame
        Original (OSM) streets data.
    base : geopandas.GeoDataFrame
        Baseline street data -- can be any (other than ``orig``).
    comp : geopandas.GeoDataFrame
        Comparative street data -- can be any (other than ``orig`` or ``base``).

    Returns
    -------
    fig : matplotlib.figure.Figure
        Plot figure.
    axs : dict
        Subplot mosaic container.
    """

    if isinstance(grid_id, str):
        geom = grid.loc[(grid["hex_id"] == grid_id)].geometry.squeeze()
    else:
        geom = grid.loc[grid_id].geometry

    fig, axs = matplotlib.pyplot.subplot_mosaic(
        [["orig", "base", "comp"], ["agg", "agg", "agg"]],
        sharex=True,
        sharey=True,
        figsize=(10, 14),
    )

    _orig = geopandas.clip(orig, geom).copy()
    _base = geopandas.clip(base, geom).copy()
    _comp = geopandas.clip(comp, geom).copy()
    _grid = geopandas.clip(grid, geom).copy()

    for data, color, label, title in [
        [_orig, "black", "orig", orig_name],
        [_base, "red", "base", base_name],
        [_comp, "blue", "comp", comp_name],
    ]:
        data.plot(ax=axs[label], color=color)
        axs[label].set_title(title)

    agg = axs["agg"]
    _base.plot(ax=agg, color="red", alpha=0.5, lw=2, zorder=1, label=base_name)
    _comp.plot(ax=agg, color="blue", alpha=1, lw=2, zorder=2, label=comp_name)
    _grid.plot(ax=agg, facecolor="None", lw=0.5, ec="black", zorder=3)
    contextily.add_basemap(
        ax=agg, crs=_grid.crs, source=contextily.providers.CartoDB.Voyager, zorder=0
    )
    agg.set_title("Overlay")
    agg.legend()

    for ax in axs:
        axs[ax].set_axis_off()

    matplotlib.pyplot.subplots_adjust(wspace=0, hspace=-0.25)

    return fig, axs
