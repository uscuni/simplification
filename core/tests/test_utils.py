import geopandas.testing
import networkx
import pytest

import core


def test_read_sample_data():
    df = core.utils.read_sample_data()
    assert isinstance(df, geopandas.GeoDataFrame)
    assert df.shape == (131, 15)
    assert df.crs == "EPSG:4326"


cities = core.utils.city_fua.keys()
records = [86_101, 109_976, 110_074, 88_746, 105_174, 110_211, 107_551]


@pytest.mark.parametrize("city, n_records", zip(cities, records, strict=True))
def test_read_parquet_roads(city, n_records):
    fua = core.utils.city_fua[city]

    gdf_1 = core.utils.read_parquet_roads(fua)
    gdf_2 = core.utils.read_parquet_roads(core.utils.city_fua[core.utils.fua_city[fua]])

    geopandas.testing.assert_geodataframe_equal(gdf_1, gdf_2)

    assert gdf_1.shape[0] == gdf_2.shape[0] == n_records


def test_graph_size():
    known = "test graph\n\t* Graph with 0 nodes and 0 edges"
    observed = core.utils.graph_size("test graph", networkx.Graph())
    assert known == observed


def test_load_usecases():
    usecases, path = core.utils.load_usecases("Liège")

    known = {
        "coords": [5.592738, 50.616253],
        "address": "Carrefour Quai des Vennes",
        "comments": "...",
    }
    assert usecases["parkinglot"] == known

    known = ("usecases", "1656")
    assert known == path.parts[-2:]
