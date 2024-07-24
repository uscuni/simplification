import geopandas
import pandas
import pytest

import core


class TestContinuity:
    def setup_method(self):
        self.city = "Liège"
        self.roads = core.utils.read_no_degree_2(self.city)

    def test_basic(self):
        roads = core.algorithms.common.continuity(self.roads)

        assert isinstance(roads, geopandas.GeoDataFrame)
        assert roads.geom_type.unique()[0] == "LineString"

        assert roads.columns.tolist() == [
            "geometry",
            "coins_group",
            "coins_len",
            "coins_count",
            "coins_end",
        ]

        known_counts_coins_end = pandas.DataFrame(
            {"coins_end": [True, False], "count": [15205, 13401]}
        )
        observed_counts_coins_end = (
            roads["coins_end"].value_counts().to_frame().reset_index()
        )
        pandas.testing.assert_frame_equal(
            known_counts_coins_end, observed_counts_coins_end
        )

    def test_basic_non_depuded(self):
        with pytest.raises(
            ValueError,
            match=(
                "Lines are identical. Please revise input data to "
                "ensure no lines are identical or overlapping."
            ),
        ):
            core.algorithms.common.continuity(self.roads, dedup=False)
