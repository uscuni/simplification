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
    "read_no_degree_2",
    "read_results",
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

# dict of method: color in plot
method_color = {
    "original": "#999999",
    "cityseer": "#11A579",
    "osmnx": "#7F3C8D",
    "parenx-voronoi": "#E73F74",
    "parenx-skeletonize": "#3969AC",
    "neatnet": "#F2B701",
    "neatnet-buildings": "#008695",
    "manual": "#80BA5A",
}

# dict of use cases
# dict of use cases
use_cases = {
    0: {"coords": (9.81981, 4.00249), "buffer": 220, "situation": "Parallel roads"},
    1: {"coords": (9.82132, 4.00986), "buffer": 100, "situation": "Roundabouts"},
    2: {"coords": (37.06837, 36.31436), "buffer": 200, "situation": "Diverging roads"},
    3: {"coords": (174.7512, -36.88972), "buffer": 100, "situation": "T-junction"},
    4: {
        "coords": (37.06599, 36.31055),
        "buffer": 120,
        "situation": "Simple intersection",
    },
    5: {
        "coords": (9.74106, 4.09687),
        "buffer": 80,
        "situation": "Cross-shaped intersection",
    },
    6: {"coords": (37.1427, 36.23658), "buffer": 30, "situation": "Intersection"},
    7: {"coords": (37.13986, 36.24012), "buffer": 220, "situation": "Side roads"},
    8: {"coords": (5.631484, 50.605893), "buffer": 110, "situation": "Cul-de-sac"},
    9: {"coords": (37.16807, 36.19395), "buffer": 130, "situation": "Ovalabout"},
    10: {
        "coords": (37.22222, 36.19417),
        "buffer": 200,
        "situation": "Cloverleaf interchange",
    },
    11: {
        "coords": (174.84, -36.91861),
        "buffer": 200,
        "situation": "Multi-level carriageway",
    },
    12: {
        "coords": (5.61551, 50.67645),
        "buffer": 100,
        "situation": "Special case roundabouts",
    },
    13: {
        "coords": (-73.16139, 7.06444),
        "buffer": 170,
        "situation": "Parallel and linking roads",
    },
    14: {"coords": (-73.1275, 7.12528), "buffer": 90, "situation": "Outliers"},
    15: {
        "coords": (9.66117, 4.08809),
        "buffer": 200,
        "situation": "Parallel road leading to different levels",
    },
    16: {
        "coords": (37.17978, 36.20862),
        "buffer": 100,
        "situation": "Multi-level roundabout",
    },
    17: {
        "coords": (36.99121, 36.09677),
        "buffer": 560,
        "situation": "Partial cloverleaf interchange",
    },
    18: {
        "coords": (174.7616667, -36.7991667),
        "buffer": 550,
        "situation": "Complicated freeway intersection",
    },
}


def read_sample_data() -> geopandas.GeoDataFrame:
    return geopandas.read_parquet(data_dir / f"sample.{parq}")


def _fua_code(fua: int | str) -> int:
    """Helper for converting a city name into a FUA code."""
    if isinstance(fua, str):
        fua = city_fua[fua]
    return fua


def _fua_path(fua: int, dataset: str) -> pathlib.Path:
    """Helper for parsing input dataset paths."""
    fua = _fua_code(fua)
    dset = data_dir / pathlib.Path(f"{fua}", dataset)
    return dset / f"{fua}.{parq}"


def read_original(fua: int | str, geom_only: bool = True) -> geopandas.GeoDataFrame:
    """Read OSM roads from parquet format; return bare columns."""
    return geopandas.read_parquet(
        _fua_path(fua, "original"),
        columns=["highway", "geometry"] if not geom_only else ["geometry"],
    ).reset_index(drop=True)


def read_no_degree_2(fua: int | str) -> geopandas.GeoDataFrame:
    """Read OSM roads from parquet format; return bare columns."""
    return geopandas.read_parquet(_fua_path(fua, "no_degree_2"))


def read_results(
    fua: int | str, method: str, proj_crs: str | int | pyproj.CRS
) -> geopandas.GeoDataFrame:
    """Read results parquet format.

    Parameters
    -----------
    fua : int | str
        Code or name of FUA.
    method : str
        Simplification method used.
    proj_crs : str | int | pyproj.CRS
        Projected CRS information.

    Returns
    --------
    geopandas.GeoDataFrame
        Simplified results from ``fua`` by ``method``.
    """

    return (
        geopandas.read_parquet(_fua_path(fua, method))
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
    df = momepy.remove_false_nodes(read_original(fua))
    df = df[~df.geometry.duplicated()].reset_index(drop=True)
    return momepy.remove_false_nodes(df)
