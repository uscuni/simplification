import geopandas.testing
import networkx
import numpy
import pytest

import core


def test_read_sample_data():
    df = core.utils.read_sample_data()
    assert isinstance(df, geopandas.GeoDataFrame)
    assert df.shape == (131, 15)
    assert df.crs == pytest.epsg_4326


cities = core.utils.city_fua.keys()
records = [86_101, 109_976, 110_074, 88_746, 105_174, 110_211, 107_551]


@pytest.mark.parametrize("city, n_records", zip(cities, records, strict=True))
def test_read_parquet_roads(city, n_records):
    fua = core.utils.city_fua[city]

    gdf_1 = core.utils.read_parquet_roads(fua)
    gdf_2 = core.utils.read_parquet_roads(core.utils.city_fua[core.utils.fua_city[fua]])

    geopandas.testing.assert_geodataframe_equal(gdf_1, gdf_2)

    assert gdf_1.shape[0] == gdf_2.shape[0] == n_records


def test_read_manual():
    fua = pytest.auckland

    observed = core.utils.read_manual(fua, pytest.epsg_4326)

    known_records = 35650
    assert observed.shape[0] == known_records
    assert observed.crs == pytest.epsg_4326


@pytest.mark.parametrize(
    "option, n_records", [["skeletonize", 32_029], ["voronoi", 31_251]]
)
def test_read_parenex(option, n_records):
    fua = pytest.auckland

    observed = core.utils.read_parenx(fua, option, pytest.epsg_4326)

    known_records = n_records
    assert observed.shape[0] == known_records
    assert observed.crs == pytest.epsg_4326


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


def test_make_grid():
    fua = pytest.auckland
    resolution = 8

    observed = core.utils.make_grid(fua, resolution, pytest.epsg_4326)

    assert observed.crs == pytest.epsg_4326
    known_records = 294
    assert observed.shape[0] == known_records
    known_bounds = [
        174.6123711149229,
        -37.01764879430995,
        174.90988694906858,
        -36.74330035902645,
    ]
    numpy.testing.assert_array_equal(observed.total_bounds, known_bounds)

    known_hexid_len = 15
    assert observed["hex_id"].map(lambda x: len(x)).unique()[0] == known_hexid_len
