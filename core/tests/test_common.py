import geopandas

import core


class TestContinuity:
    def setup_method(self):
        self.city = "Liège"
        self.roads = core.utils.read_no_degree_2_roads(self.city)

    def test_basic(self):
        roads = core.algorithms.common.continuity(self.roads)

        assert isinstance(roads, geopandas.GeoDataFrame)
        assert roads.geom_type.unique()[0] == "LineString"
