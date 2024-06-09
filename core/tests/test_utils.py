import pathlib
import shutil

import geopandas.testing
import networkx
import pytest
import shapely.testing
from matplotlib.testing.decorators import image_comparison

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


@image_comparison(baseline_images=["test_plot_movie.png"], style="mpl20")
def test_plot_movie():
    # declare AOI
    city = "Liège"
    roads_bare = core.utils.read_parquet_roads(core.utils.city_fua[city])
    orig_crs = roads_bare.crs

    # declare singular param set
    merge_midline = True
    contains_buffer_dist = 17

    points, fpath_base = core.utils.load_usecases(city)
    package = "cityseer"
    fpath_pack = fpath_base / package
    myclass = "parallel_edges_1"

    # isolate central AOI
    mypoint = points[myclass]["coords"]
    center = core.utils.viz_class_location(mypoint, orig_crs)

    # read in prepared subset for testing
    gpkg = pathlib.Path("core", "tests", "data", "test_data_liege.gpkg")
    _layer = "cityseer_parallel_edges_1_midline_True_17_"
    nodes = geopandas.read_file(gpkg, layer=f"{_layer}nodes")
    edges = geopandas.read_file(gpkg, layer=f"{_layer}edges")

    # make subfolder for plot saving
    _myclass = f"TESTING_{myclass}_midline_{merge_midline}"
    fpath_class = core.utils.viz_class_path(_myclass, fpath_pack)

    # make class-param plot
    buff_pad_3 = f"{contains_buffer_dist:03d}"
    fmt_title = (
        f"Contains Buffer Dist: {buff_pad_3}m\n"
        f"Merge Edges by Midline: {merge_midline}"
    )
    fmt_fname = f"{buff_pad_3}_{merge_midline}"
    core.utils.viz_class_param_plot(
        nodes,
        edges,
        center,
        fmt_title,
        fmt_fname,
        fpath_class,
        orig_crs,
        close=False,
    )

    # make video
    core.utils.viz_class_video(fpath_class)

    # plot existence & removal
    out_png_path = fpath_class / f"{fmt_fname}.png"
    assert out_png_path.exists()
    shutil.rmtree(fpath_class)
    assert not fpath_class.exists()

    # movie existence & removal
    out_mp4_path = fpath_class.with_suffix(".mp4")
    assert out_mp4_path.exists()
    out_mp4_path.unlink()
    assert not out_mp4_path.exists()
