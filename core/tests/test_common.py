import geopandas
import momepy

import core


class TestContinuity:
    def setup_method(self):
        self.city = "Liège"
        self.roads = momepy.remove_false_nodes(core.utils.read_parquet_roads(self.city))

    def test_basic(self):
        roads = core.algorithms.common.continuity(self.roads)

        assert isinstance(roads, geopandas.GeoDataFrame)
        assert roads.geom_type.unique()[0] == "LineString"
