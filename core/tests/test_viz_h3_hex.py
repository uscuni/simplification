import pytest
from matplotlib.testing.decorators import image_comparison

import core


@image_comparison(
    baseline_images=["test_viz_h3_hex_plot_aoi.png"],
    style="mpl20",
    tol=7 if pytest.ENV_TYPE == "dev" else 1,  # (RMS 6.914)
)
def test_viz_h3_hex_plot_aoi(grid_8_auckland):
    core.viz.h3_hex.plot_aoi(grid_8_auckland, "Auckland", 8)


@image_comparison(
    baseline_images=["test_viz_h3_hex_plot_analysis.png"],
    style="mpl20",
    tol=35 if pytest.ENV_TYPE == "dev" else 1,  # (RMS 34.017)
)
def test_viz_h3_hex_plot_analysis(manual_auckland, parenx_auckland, grid_7_auckland):
    _, _, edges_manual = manual_auckland
    _, _, edges_parenx = parenx_auckland
    grid = grid_7_auckland

    _m = "manual"
    _p = "parenx"
    var = "edge_count"
    for dname, df in [[_m, edges_manual], [_p, edges_parenx]]:
        grid[[f"{var}_{dname}", "IGNORE", "ALSO_IGNORE"]] = grid.apply(
            lambda x: core.stats.get_edge_stats(df, x.geometry),  # noqa: B023
            axis=1,
            result_type="expand",
        )

    grid[f"{var}_ratio"] = grid[f"{var}_{_p}"] / grid[f"{var}_{_m}"]

    core.viz.h3_hex.plot_analysis(grid, f"{var}_ratio", info=f"({_p}/{_m})")


tol = 5 if pytest.ENV_TYPE == "dev" else 1  # (RMS 4.030)


class TestVizH3HexPlotCell:
    @pytest.fixture(autouse=True)
    def setup_method(
        self, osm_auckland, manual_auckland, parenx_auckland, grid_9_auckland
    ):
        self.edges_osm = osm_auckland
        _, _, self.edges_manual = manual_auckland
        _, _, self.edges_parenx = parenx_auckland
        self.grid = grid_9_auckland

        self._o = "osm"
        self._m = "manual"
        self._p = "parenx-voronoi"

    @image_comparison(
        baseline_images=["test_viz_h3_hex_plot_cell_loc_id.png"], style="mpl20", tol=tol
    )
    def test_loc_id(self):
        core.viz.h3_hex.plot_cell(
            84,
            self.grid,
            self.edges_osm,
            self._o,
            self.edges_manual,
            self._m,
            self.edges_parenx,
            self._p,
        )

    @image_comparison(
        baseline_images=["test_viz_h3_hex_plot_cell_hex_id.png"],
        style="mpl20",
        tol=tol,
    )
    def test_hex_id(self):
        core.viz.h3_hex.plot_cell(
            "89bb5000563ffff",
            self.grid,
            self.edges_osm,
            self._o,
            self.edges_manual,
            self._m,
            self.edges_parenx,
            self._p,
        )
