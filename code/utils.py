import pathlib

import geopandas
import networkx

__all__ = [
    "city_fua",
    "fua_city",
    "read_sample_data",
    "read_parquet_roads",
    "graph_size",
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
