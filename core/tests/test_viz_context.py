import pathlib
import shutil

import geopandas.testing
import pytest
import shapely.testing
from matplotlib.testing.decorators import image_comparison

import core


def test_viz_context_path():
    context_for_test = "some_context"
    dir_for_test_str = "test_viz_context_path"
    dir_for_test_path = pathlib.Path(dir_for_test_str)

    known = pathlib.Path(dir_for_test_str, context_for_test)
    observed = core.viz.context.path(context_for_test, dir_for_test_path)

    assert known == observed
    assert observed.exists()
    shutil.rmtree(dir_for_test_path)
    assert not observed.exists()


def test_viz_context_location():
    known_crs = pytest.epsg_3857
    known_area = 10000.0
    known_x, known_y = -82.3347, 27.214658
    known_point = shapely.Point(known_x, known_y)

    observed = core.viz.context.location(known_point, crs=known_crs, buffer=50)

    assert known_crs == observed.crs
    assert known_area == observed.squeeze().area
    shapely.testing.assert_geometries_equal(
        known_point, observed.centroid.to_crs(pytest.epsg_4326).squeeze()
    )


# RMS (12.174 // 255) for Ubuntu DEV --> MPL via pypi (very small diff)
tol = 12.2 if pytest.ENV_TYPE == "dev" else 1


@image_comparison(
    baseline_images=["test_viz_context_param_plot_video.png"], style="mpl20", tol=tol
)
def test_viz_context_param_plot_video():
    # declare AOI
    city = "Liège"
    roads_bare = core.utils.read_original(core.utils.city_fua[city])
    orig_crs = roads_bare.crs

    # declare singular param set
    merge_midline = True
    contains_buffer_dist = 17

    points, fpath_base = core.utils.load_usecases(city)
    package = "cityseer"
    fpath_pack = fpath_base / package
    _context = "parallel_edges_1"

    # isolate central AOI
    point = points[_context]["coords"]
    center = core.viz.context.location(point, orig_crs)

    # read in prepared subset for testing
    gpkg = pathlib.Path("core", "tests", "data", "test_data_liege.gpkg")
    _layer = "cityseer_parallel_edges_1_midline_True_17_"
    nodes = geopandas.read_file(gpkg, layer=f"{_layer}nodes")
    edges = geopandas.read_file(gpkg, layer=f"{_layer}edges")

    # make subfolder for plot saving
    __context = f"TESTING_{_context}_midline_{merge_midline}"
    fpath_context = core.viz.context.path(__context, fpath_pack)

    # make class-param plot
    buff_pad_3 = f"{contains_buffer_dist:03d}"
    fmt_title = (
        f"Contains Buffer Dist: {buff_pad_3}m\nMerge Edges by Midline: {merge_midline}"
    )
    fmt_fname = f"{buff_pad_3}_{merge_midline}"
    core.viz.context.param_plot(
        nodes,
        edges,
        center,
        fmt_title,
        fmt_fname,
        fpath_context,
        orig_crs,
        close=False,
    )

    # make video
    core.viz.context.video(fpath_context)

    # plot existence & removal
    out_png_path = fpath_context / f"{fmt_fname}.png"
    assert out_png_path.exists()
    shutil.rmtree(fpath_context)
    assert not fpath_context.exists()

    # movie existence & removal
    out_mp4_path = fpath_context.with_suffix(".mp4")
    assert out_mp4_path.exists()
    out_mp4_path.unlink()
    assert not out_mp4_path.exists()
