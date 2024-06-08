import json
import pathlib

import geopandas
import networkx
import pyproj
import shapely

__all__ = [
    "city_fua",
    "fua_city",
    "read_sample_data",
    "read_parquet_roads",
    "graph_size",
    "load_usecases",
    "viz_class_path",
    "viz_class_location",
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


##################
# basic utils
##################

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


def read_parquet_roads(fua: int) -> geopandas.GeoDataFrame:
    """Read OSM roads from parquet format; return bare columns."""
    return geopandas.read_parquet(data_dir / f"{fua}" / f"roads_osm.{parq}")[
        ["highway", "geometry"]
    ].reset_index(drop=True)


def graph_size(info: str, g: networkx.Graph) -> str:
    return f"{info}\n\t* {g}"


def load_usecases(city: str) -> tuple[dict, pathlib.Path]:
    """Load use cases for a city and also return the directory path."""
    path_base = usecase_dir / str(city_fua[city])
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
        geopandas.GeoDataFrame(geometry=[shapely.Point(mypoint)], crs="epsg:4326")
        .to_crs(crs)
        .buffer(buffer, cap_style=3)
    )
