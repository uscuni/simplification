import json
import pathlib

import geopandas
import momepy
import networkx
import pyproj
import tobler

__all__ = [
    "city_fua",
    "fua_city",
    "read_sample_data",
    "read_original",
    "read_no_degree_2_roads",
    "read_manual",
    "read_parenx",
    "graph_size",
    "load_usecases",
    "make_grid",
    "remove_degree_2_nodes",
]

parq = "parquet"

# testing hack... MUST BE EDITABLE INSTALL < probably need to rethink...
curr_path = pathlib.Path(__file__).resolve()
if curr_path.parts[-1] == "utils.py":
    top_dir = curr_path.parents[1]
else:
    top_dir = curr_path / "simplification"
data_dir = top_dir / "data"
usecase_dir = top_dir / "usecases"


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
    return geopandas.read_parquet(data_dir / f"sample.{parq}")


def read_original(fua: int | str, geom_only: bool = True) -> geopandas.GeoDataFrame:
    """Read OSM roads from parquet format; return bare columns."""
    if isinstance(fua, str):
        fua = city_fua[fua]
    _fua_path = data_dir / f"{fua}" / f"roads_osm.{parq}"
    cols = ["highway", "geometry"] if not geom_only else ["geometry"]
    return geopandas.read_parquet(_fua_path, columns=cols).reset_index(drop=True)


def read_no_degree_2_roads(fua: int | str) -> geopandas.GeoDataFrame:
    """Read OSM roads from parquet format; return bare columns."""
    if isinstance(fua, str):
        fua = city_fua[fua]
    _fua_path = data_dir / pathlib.Path(f"{fua}", "no_degree_2", f"{fua}.{parq}")
    return geopandas.read_parquet(_fua_path)


def read_manual(fua: int, proj_crs: str | int | pyproj.CRS) -> geopandas.GeoDataFrame:
    """Read in manually prepared simplified road data."""
    if isinstance(fua, str):
        fua = city_fua[fua]
    _fua_path = data_dir / pathlib.Path(f"{fua}", "manual", f"{fua}.{parq}")
    return (
        geopandas.read_parquet(_fua_path)[["geometry"]]
        .explode(ignore_index=True, index_parts=False)
        .to_crs(proj_crs)
    )


def read_parenx(
    fua: int, option: str, proj_crs: str | int | pyproj.CRS
) -> geopandas.GeoDataFrame:
    """Read in prepared parenx data."""
    if isinstance(fua, str):
        fua = city_fua[fua]
    _fua_path = data_dir / pathlib.Path(f"{fua}", "parenx", f"{option}.{parq}")
    return (
        geopandas.read_parquet(_fua_path)
        .explode(ingore_index=True, index_parts=False)
        .to_crs(proj_crs)
    )


def graph_size(info: str, g: networkx.Graph) -> str:
    return f"{info}\n\t* {g}"


def load_usecases(city: str) -> tuple[dict, pathlib.Path]:
    """Load use cases for a city and also return the directory path."""
    path_base = usecase_dir / str(city_fua[city])
    with open(path_base / "points.json") as f:
        points = json.load(f)
    return points, path_base


def make_grid(
    fua: int,
    res: int,
    proj_crs: str | int | pyproj.CRS,
    base_crs: str | int | pyproj.CRS = "EPSG:4326",
) -> geopandas.GeoDataFrame:
    """Create a hex grid for given FUA (using original OSM dataset) and resolution.

    Parameters
    ----------
    fua : int
        ID of city (any of use cases: 809, 869, 1133, 1656, 4617, 4881, 8989)
    res : int
        H3 hex grid resolution. Suggested resolutions:
            * 8 - with cell side length ~530m
            * 9 - cell side length ~200m
    proj_crs :  str | int | pyproj.CRS
        Target projection.
    base_crs :  str | int | pyproj.CRS (default EPSG:4326)
        Default system needed for H3 operability.

    Returns
    -------
    geopandas.GeoDataFrame
        H3 hexbins covering the FUA that contain data.
    """

    # get city geom and name
    meta = read_sample_data()
    geom = meta.loc[meta.eFUA_ID == fua, "geometry"].copy()

    # read in OSM data
    orig = read_original(fua).to_crs(base_crs)
    assert meta.crs == orig.crs, "CRS of 'meta' and 'orig' are not equal."

    grid = tobler.util.h3fy(
        source=geom, resolution=res, clip=False, buffer=False, return_geoms=True
    )

    # set hex id in columns & return only the grid cells that contain data
    return (
        grid.assign(hex_id=grid.index)
        .reset_index(drop=True)
        .loc[
            list(
                set(
                    orig.sindex.query(grid.geometry, predicate="intersects", sort=True)[
                        0
                    ]
                )
            )
        ]
        .sort_values(by="hex_id")
        .to_crs(proj_crs)
        .reset_index(drop=True)
    )


def remove_degree_2_nodes(fua: int | str) -> geopandas.GeoDataFrame:
    """Remove [interstitial / non-articulation / degree 2] nodes from road network."""
    return momepy.remove_false_nodes(read_original(fua))
