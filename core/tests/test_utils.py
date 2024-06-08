import pathlib
import shutil

import geopandas.testing
import networkx
import pytest
import shapely.testing

import core


def test_read_sample_data():
    df = core.utils.read_sample_data()
    assert isinstance(df, geopandas.GeoDataFrame)
    assert df.shape == (131, 15)
    assert df.crs == "EPSG:4326"


cities = core.utils.city_fua.keys()
records = [86_101, 109_976, 110_074, 88_746, 105_174, 110_211, 107_551]


# basic utils tests ------------------------------


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


# Viz funcs tests ------------------------------
def test_viz_class_path():
    class_for_test = "some_class"
    dir_for_test_str = "test_viz_class_path"
    dir_for_test_path = pathlib.Path(dir_for_test_str)

    known = pathlib.Path(dir_for_test_str, class_for_test)
    observed = core.utils.viz_class_path(class_for_test, dir_for_test_path)

    assert known == observed
    assert observed.exists()
    shutil.rmtree(dir_for_test_path)
    assert not observed.exists()


def test_viz_class_location():
    known_crs = "EPSG:3857"
    known_area = 10000.0
    known_x, known_y = -82.3347, 27.214658
    known_point = shapely.Point(known_x, known_y)

    observed = core.utils.viz_class_location(known_point, crs=known_crs, buffer=50)

    assert known_crs == observed.crs
    assert known_area == observed.squeeze().area
    shapely.testing.assert_geometries_equal(
        known_point, observed.centroid.to_crs("EPSG:4326").squeeze()
    )
