import momepy
import pytest

import core


def pytest_addoption(parser):
    """Add custom command line arguments."""

    # flag for determining CI environment
    parser.addoption(
        "--env_type",
        action="store",
        default="latest",
        help="Testing environment type label",
        type=str,
    )


def pytest_configure(config):
    """PyTest session attributes, methods, etc."""

    env_type = config.getoption("env_type")

    # set ``ENV_TYPE``
    valid_env_types = ["oldest", "latest", "dev"]
    pytest.ENV_TYPE = env_type.split("-")[-1].split(".")[0]
    assert pytest.ENV_TYPE in valid_env_types

    pytest.epsg_4326 = "EPSG:4326"
    pytest.epsg_3857 = "EPSG:3857"

    pytest.auckland = core.utils.city_fua["Auckland"]


def _osm_auckland():
    return core.utils.read_parquet_roads(pytest.auckland)


@pytest.fixture
def osm_auckland():
    """Helper for reading in original OSM Auckland network."""
    return _osm_auckland()


@pytest.fixture
def manual_auckland():
    """Helper for setup of 'manual' Auckland network."""
    roads = core.utils.read_manual(pytest.auckland, _osm_auckland().crs)
    graph = momepy.gdf_to_nx(roads, length="length", integer_labels=True)
    nodes, edges = momepy.nx_to_gdf(graph)
    return graph, nodes, edges


@pytest.fixture
def parenx_auckland():
    """Helper for setup of parenx voronoi Auckland network."""
    roads = core.utils.read_parenx(pytest.auckland, "voronoi", _osm_auckland().crs)
    graph = momepy.gdf_to_nx(roads, length="length", integer_labels=True)
    nodes, edges = momepy.nx_to_gdf(graph)
    return graph, nodes, edges


@pytest.fixture
def grid_7_auckland():
    """Helper for generation of Auckland H3 grid at resolution 7."""
    return core.utils.make_grid(pytest.auckland, 7, _osm_auckland().crs)


@pytest.fixture
def grid_8_auckland():
    """Helper for generation of Auckland H3 grid at resolution 8."""
    return core.utils.make_grid(pytest.auckland, 8, _osm_auckland().crs)


@pytest.fixture
def grid_9_auckland():
    """Helper for generation of Auckland H3 grid at resolution 9."""
    return core.utils.make_grid(pytest.auckland, 9, _osm_auckland().crs)
