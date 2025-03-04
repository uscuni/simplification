import pandas
import pytest

import core


def test_add_node_degree(manual_auckland):
    manual_graph, nodes_manual, _ = manual_auckland

    observed = core.stats.add_node_degree(nodes_manual, manual_graph)["degree"]

    assert isinstance(observed, pandas.Series)
    assert observed.shape[0] == 6755
    assert observed.sum() == 17280


def test_get_edge_stats(manual_auckland, grid_9_auckland):
    _, _, edges_manual = manual_auckland
    grid_cell = grid_9_auckland[grid_9_auckland["hex_id"] == "89bb50031c7ffff"].geometry

    known_count, known_length, known_coord_count = (5, 551.6152216584356, 15)
    observed_count, observed_length, observed_coord_count = core.stats.get_edge_stats(
        edges_manual, grid_cell
    )

    assert observed_count == known_count
    assert known_length == pytest.approx(observed_length)
    assert known_coord_count == observed_coord_count


def test_avg_degree():
    assert not core.stats._avg_degree({})
    assert core.stats._avg_degree({1: 1, 2: 1}) == 1.5


def test_get_node_stats(manual_auckland, grid_9_auckland):
    manual_graph, nodes_manual, _ = manual_auckland
    grid_cell = grid_9_auckland[grid_9_auckland["hex_id"] == "89bb50031c7ffff"].geometry

    known_count, known_distr, known_avg = (4, {1: 3, 3: 1}, 1.5)
    observed_count, observed_distr, observed_avg = core.stats.get_node_stats(
        core.stats.add_node_degree(nodes_manual, manual_graph), grid_cell
    )

    assert observed_count == known_count
    assert observed_distr == known_distr
    assert observed_avg == known_avg
