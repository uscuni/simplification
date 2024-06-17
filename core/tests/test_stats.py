import pandas

import core


def test_add_node_degree(manual_auckland):
    manual_graph, nodes_manual, _ = manual_auckland

    observed = core.stats.add_node_degree(nodes_manual, manual_graph)["degree"]

    assert isinstance(observed, pandas.Series)
    assert observed.shape[0] == 31050
    assert observed.sum() == 71300


def test_get_edge_stats(manual_auckland, grid_9_auckland):
    _, _, edges_manual = manual_auckland
    grid_cell = grid_9_auckland[grid_9_auckland["hex_id"] == "89bb50031c7ffff"].geometry

    known_count, known_length = (5, 578.75495776233)
    observed_count, observed_length = core.stats.get_edge_stats(edges_manual, grid_cell)

    assert observed_count == known_count
    assert known_length == observed_length


def test_avg_degree():
    assert not core.stats._avg_degree({})
    assert core.stats._avg_degree({1: 1, 2: 1}) == 1.5


def test_get_node_stats(manual_auckland, grid_9_auckland):
    manual_graph, nodes_manual, _ = manual_auckland
    grid_cell = grid_9_auckland[grid_9_auckland["hex_id"] == "89bb50031c7ffff"].geometry

    known_count, known_distr, known_avg = (3, {1: 2, 3: 1}, 1.6666666666666667)
    observed_count, observed_distr, observed_avg = core.stats.get_node_stats(
        core.stats.add_node_degree(nodes_manual, manual_graph), grid_cell
    )

    assert observed_count == known_count
    assert observed_distr == known_distr
    assert observed_avg == known_avg
