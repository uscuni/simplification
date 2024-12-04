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


cities = list(core.utils.city_fua.keys())
osm_records = [78_908, 60_364, 79_317, 84_819, 79_907, 50_917, 92_667]
xnd2_records = [43_233, 12_439, 16_302, 30_552, 16_554, 13_468, 29_314]
man_records = [38_772, 8_640, 14_170, 29_252, 13_508, 11_032, numpy.nan]
parenx_voronoi_records = [42_302, 9_569, 16_844, 29_446, 15_516, 14_242, 35_694]
parenx_skeletonize_records = [44_294, 9_561, 17_469, 30_641, 16_075, 14_784, 37_557]
osmnx_records = [43_451, 15_770, 24_025, 37_954, 21_755, 17_704, 30_156]
neatnet_records = [39_489, 8_229, 13_587, 29_028, 12_395, 11_059, 16_298]


@pytest.mark.parametrize("city, n_records", zip(cities, osm_records, strict=True))
def test_read_original(city, n_records):
    fua = core.utils.city_fua[city]

    gdf_1 = core.utils.read_original(fua)
    gdf_2 = core.utils.read_original(core.utils.fua_city[fua])

    geopandas.testing.assert_geodataframe_equal(gdf_1, gdf_2)

    assert gdf_1.shape[0] == gdf_2.shape[0] == n_records


@pytest.mark.parametrize("city, n_records", zip(cities, xnd2_records, strict=True))
def test_read_no_degree_2(city, n_records):
    fua = core.utils.city_fua[city]

    gdf_1 = core.utils.read_no_degree_2(fua)
    gdf_2 = core.utils.read_no_degree_2(core.utils.fua_city[fua])

    geopandas.testing.assert_geodataframe_equal(gdf_1, gdf_2)

    assert gdf_1.shape[0] == gdf_2.shape[0] == n_records


@pytest.mark.parametrize(
    "city, n_records", zip(cities[:-1], man_records[:-1], strict=True)
)
def test_read_manual(city, n_records):
    fua = core.utils.city_fua[city]

    gdf_1 = core.utils.read_manual(fua, pytest.epsg_4326)
    gdf_2 = core.utils.read_manual(core.utils.fua_city[fua], pytest.epsg_4326)

    geopandas.testing.assert_geodataframe_equal(gdf_1, gdf_2)

    assert gdf_1.shape[0] == gdf_2.shape[0] == n_records
    assert gdf_1.crs == gdf_2.crs == pytest.epsg_4326


@pytest.mark.parametrize("city, n_records", zip(cities, osmnx_records, strict=True))
def test_read_osmnx(city, n_records):
    fua = core.utils.city_fua[city]

    gdf_1 = core.utils.read_osmnx(fua, pytest.epsg_4326)
    gdf_2 = core.utils.read_osmnx(core.utils.fua_city[fua], pytest.epsg_4326)

    geopandas.testing.assert_geodataframe_equal(gdf_1, gdf_2)

    assert gdf_1.shape[0] == gdf_2.shape[0] == n_records
    assert gdf_1.crs == gdf_2.crs == pytest.epsg_4326


@pytest.mark.parametrize(
    "city, n_records", zip(cities, parenx_voronoi_records, strict=True)
)
def test_read_parenx_voronoi(city, n_records):
    fua = core.utils.city_fua[city]
    option = "voronoi"
    gdf_1 = core.utils.read_parenx(fua, proj_crs=pytest.epsg_4326, option=option)
    gdf_2 = core.utils.read_parenx(
        core.utils.fua_city[fua], proj_crs=pytest.epsg_4326, option=option
    )

    geopandas.testing.assert_geodataframe_equal(gdf_1, gdf_2)

    assert gdf_1.shape[0] == gdf_2.shape[0] == n_records
    assert gdf_1.crs == gdf_2.crs == pytest.epsg_4326


@pytest.mark.parametrize(
    "city, n_records", zip(cities, parenx_skeletonize_records, strict=True)
)
def test_read_parenx_skeletonize(city, n_records):
    fua = core.utils.city_fua[city]
    option = "skeletonize"
    gdf_1 = core.utils.read_parenx(fua, proj_crs=pytest.epsg_4326, option=option)
    gdf_2 = core.utils.read_parenx(
        core.utils.fua_city[fua], proj_crs=pytest.epsg_4326, option=option
    )

    geopandas.testing.assert_geodataframe_equal(gdf_1, gdf_2)

    assert gdf_1.shape[0] == gdf_2.shape[0] == n_records
    assert gdf_1.crs == gdf_2.crs == pytest.epsg_4326


@pytest.mark.parametrize("city, n_records", zip(cities, neatnet_records, strict=True))
def test_read_neatnet(city, n_records):
    fua = core.utils.city_fua[city]

    gdf_1 = core.utils.read_neatnet(fua, pytest.epsg_4326)
    gdf_2 = core.utils.read_neatnet(core.utils.fua_city[fua], pytest.epsg_4326)

    geopandas.testing.assert_geodataframe_equal(gdf_1, gdf_2)

    assert gdf_1.shape[0] == gdf_2.shape[0] == n_records
    assert gdf_1.crs == gdf_2.crs == pytest.epsg_4326


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
    known_records = 285
    assert observed.shape[0] == known_records
    known_bounds = [
        174.6123711149229,
        -37.01764879430995,
        174.90988694906858,
        -36.74330035902645,
    ]
    numpy.testing.assert_array_almost_equal(observed.total_bounds, known_bounds)

    known_hexid_len = 15
    assert observed["hex_id"].map(lambda x: len(x)).unique()[0] == known_hexid_len


def test_remove_degree_2_nodes():
    fua = core.utils.city_fua["Aleppo"]

    known = core.utils.read_no_degree_2(fua)
    observed = core.utils.remove_degree_2_nodes(fua)

    geopandas.testing.assert_geodataframe_equal(known, observed)
