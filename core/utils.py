import json
import pathlib

import geopandas
import networkx

__all__ = [
    "city_fua",
    "fua_city",
    "read_sample_data",
    "read_parquet_roads",
    "graph_size",
    "load_usecases",
]

parq = "parquet"

top_dir = pathlib.Path(__file__).absolute().parents[1]
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
