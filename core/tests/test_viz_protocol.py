import pathlib

import geopandas
import pytest
import shapely
from matplotlib.testing.decorators import image_comparison

import core


@image_comparison(
    baseline_images=["test_plot_case.png"],
    style="mpl20",
    tol=6.1 if pytest.ENV_TYPE == "dev" else 1,
)
def test_plot_case():
    p1, p2 = shapely.Point(20, 20), shapely.Point(21, 21)
    lines = geopandas.GeoDataFrame(geometry=[shapely.LineString((p1, p2))])
    verts = geopandas.GeoDataFrame(geometry=[p1, p2])

    protocol_case, protocol_type = "test_case", "test_type"
    image_fpath = pathlib.Path(f"{protocol_case}_{protocol_type}.png")
    core.viz.protocol.plot_case(lines, verts, image_fpath)

    image_fpath.unlink()
