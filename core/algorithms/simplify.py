import geopandas as gpd
import momepy
import numpy as np
import pandas as pd
import shapely
from libpysal import graph
from scipy import sparse

from ..geometry import is_within, snap_to_targets, voronoi_skeleton

__all__ = [
    "ccss_special_case",
    "filter_connections",
    "avoid_forks",
    "reconnect",
    "remove_dangles",
    "one_remaining",
    "multiple_remaining",
    "one_remaining_c",
    "loop",
]


def ccss_special_case(
    primes, conts_groups, highest_hierarchy, relevant_nodes, split_points
):
    # If there are primes on both Cs, connect them. If there's prime on one C,
    # connect it to the other C. If there are no primes, get midpoints on Cs and
    # connect those.
    if primes.empty:
        # midpoints solution
        c0 = conts_groups.geometry.iloc[0]
        c1 = conts_groups.geometry.iloc[1]
        p0 = shapely.line_interpolate_point(c0, 0.5, normalized=True)
        p1 = shapely.line_interpolate_point(c1, 0.5, normalized=True)
        new_connections = [shapely.LineString([p0, p1])]
        split_points.append(p0)
        split_points.append(p1)

    # one prime, get shortest line to the other C
    elif primes.shape[0] == 1:
        no_prime_c = conts_groups[conts_groups.disjoint(primes.geometry.item())]
        sl = shapely.shortest_line(primes.geometry.item(), no_prime_c.geometry.item())
        new_connections = [sl]
        split_points.append(shapely.get_point(sl, -1))

    # two primes, connect them
    elif primes.shape[0] == 2:
        new_connections = [
            shapely.shortest_line(primes.geometry.iloc[0], primes.geometry.iloc[1])
        ]

    # multiple primes, connect two nearest on distinct Cs
    else:
        primes_on_c0 = primes[primes.interesects(conts_groups.geometry.iloc[0])]
        primes_on_c1 = primes[primes.interesects(conts_groups.geometry.iloc[1])]
        new_connections = [
            shapely.shortest_line(primes_on_c0.union_all(), primes_on_c1.union_all())
        ]

    # some nodes may have ended unconnected. Find them and reconnect them.
    combined_linework = pd.concat(
        [
            highest_hierarchy,
            gpd.GeoSeries(new_connections, crs=highest_hierarchy.crs),
        ]
    ).union_all()
    missing = relevant_nodes[relevant_nodes.disjoint(combined_linework)]
    new_connections.extend(
        shapely.shortest_line(missing.geometry, combined_linework).tolist()
    )

    return new_connections


def filter_connections(primes, conts_groups, new_connections):
    # The skeleton returns connections to all the nodes. We need to keep only
    # some, if there are multiple connections to a single C. We don't touch
    # the other.

    unwanted = []
    keeping = []
    for c in conts_groups.geometry:
        int_mask = shapely.intersects(new_connections, c)
        connections_intersecting_c = new_connections[int_mask]
        if len(connections_intersecting_c) > 1:
            prime_mask = shapely.intersects(
                connections_intersecting_c, primes.union_all()
            )
            connections_intersecting_primes = connections_intersecting_c[prime_mask]
            # if there are multiple connections to a single C, drop them and keep only
            # the shortest one leading to prime
            if (
                len(connections_intersecting_c) > 1
                and len(connections_intersecting_primes) > 0
            ):
                lens = shapely.length(connections_intersecting_primes)
                unwanted.append(connections_intersecting_c)
                keeping.append(connections_intersecting_primes[[np.argmin(lens)]])

            # fork on two nodes on C
            elif len(connections_intersecting_c) > 1:
                lens = shapely.length(connections_intersecting_c)
                unwanted.append(connections_intersecting_c)
                keeping.append(connections_intersecting_c[[np.argmin(lens)]])

    if len(unwanted) > 0:
        if len(keeping) > 0:
            new_connections = np.concatenate(
                [
                    new_connections[
                        ~np.isin(new_connections, np.concatenate(unwanted))
                    ],
                    np.concatenate(keeping),
                ]
            )
        else:
            new_connections = new_connections[
                ~np.isin(new_connections, np.concatenate(unwanted))
            ]
    return new_connections, connections_intersecting_c, connections_intersecting_primes


def avoid_forks(
    highest_hierarchy,
    new_connections,
    relevant_targets,
    artifact,
    split_points,
):
    # mutliple Cs that are not intersecting. Avoid forks on the ends of Voronoi. If
    # one goes to relevant node, keep it. If not, remove both and replace with
    # a new shortest connection
    int_mask = shapely.intersects(new_connections, highest_hierarchy.union_all())
    targets_mask = shapely.intersects(new_connections, relevant_targets.union_all())
    new_connections = new_connections[(int_mask * targets_mask) | np.invert(int_mask)]
    cont_diss = highest_hierarchy.dissolve(highest_hierarchy.coins_group).geometry
    addition, splitters = snap_to_targets(
        new_connections,
        artifact.geometry,
        cont_diss[cont_diss.disjoint(shapely.union_all(new_connections))],
    )
    split_points.extend(splitters)
    new_connections = np.concatenate([new_connections, addition])

    return new_connections


def reconnect(conts_groups, new_connections, artifact, split_points):
    # check for disconnected Cs and reconnect
    new_connections_comps = graph.Graph.build_contiguity(
        gpd.GeoSeries(new_connections), rook=False
    ).component_labels
    new_components = gpd.GeoDataFrame(geometry=new_connections).dissolve(
        new_connections_comps
    )
    additions = []
    for c in conts_groups.geometry:
        mask = new_components.intersects(c)
        if not mask.all():
            adds, splitters = snap_to_targets(
                new_components[~mask].geometry, artifact.geometry, [c]
            )
            additions.extend(adds)
            split_points.extend(splitters)
    if len(additions) > 0:
        new_connections = np.concatenate([new_connections, additions])

    return new_connections


def remove_dangles(new_connections, artifact):
    # the drop above could've introduced a dangling edges. Remove those.

    new_connections = shapely.line_merge(new_connections)
    pts0 = shapely.get_point(new_connections, 0)
    pts1 = shapely.get_point(new_connections, -1)
    pts = shapely.buffer(np.concatenate([pts0, pts1]), 1e-6)
    all_idx, pts_idx = shapely.STRtree(pts).query(
        np.concatenate([pts, [artifact.geometry.boundary]]),
        predicate="intersects",
    )
    data = [True] * len(all_idx)
    sp = sparse.coo_array((data, (pts_idx, all_idx)), shape=(len(pts), len(pts) + 1))
    dangles = pts[sp.sum(axis=1) == 1]
    new_connections = new_connections[
        shapely.disjoint(new_connections, shapely.union_all(dangles))
    ]
    return new_connections


def one_remaining(
    relevant_targets, remaining_nodes, artifact, edges, es_mask, distance, split_points
):
    # find the nearest relevant target
    remaining_nearest, target_nearest = relevant_targets.sindex.nearest(
        remaining_nodes.geometry, return_all=False
    )
    # create a new connection as the shortest straight line
    new_connections = shapely.shortest_line(
        remaining_nodes.geometry.iloc[remaining_nearest].values,
        relevant_targets.geometry.iloc[target_nearest].values,
    )
    # check if the new connection is within the artifact
    connections_within = is_within(new_connections, artifact.geometry, 0.1)
    # if it is not within, discard it and use the skeleton instead
    if not connections_within.all():
        print("CONDITION is_within False")

        new_connections, splitters = voronoi_skeleton(
            edges[es_mask].geometry,  # use edges that are being dropped
            poly=artifact.geometry,
            snap_to=relevant_targets.geometry,  # snap to relevant node targets
            distance=distance,
            buffer=distance,  # TODO: figure out if we need this
        )
        split_points.extend(splitters)

    return new_connections


def multiple_remaining(
    edges,
    es_mask,
    artifact,
    distance,
    highest_hierarchy,
    split_points,
    snap_to,
):
    # use skeleton to ensure all nodes are naturally connected
    new_connections, splitters = voronoi_skeleton(
        edges[es_mask].geometry,  # use edges that are being dropped
        poly=artifact.geometry,
        snap_to=snap_to,  # snap to relevant node targets
        distance=distance,
        secondary_snap_to=highest_hierarchy.geometry,
        # buffer = highest_hierarchy.length.sum() * 1.2
    )
    split_points.extend(splitters)

    return remove_dangles(new_connections, artifact)


def one_remaining_c(
    remaining_nodes, highest_hierarchy, artifact, edges, es_mask, distance, split_points
):
    # create a new connection as the shortest straight line to any C
    new_connections = shapely.shortest_line(
        remaining_nodes.geometry.values,
        highest_hierarchy.union_all(),
    )
    splitters = shapely.get_point(new_connections, -1)
    # check if the new connection is within the artifact
    connections_within = is_within(new_connections, artifact.geometry, 0.1)
    # if it is not within, discard it and use the skeleton instead
    if not connections_within.all():
        print("CONDITION is_within False")

        new_connections, splitters = voronoi_skeleton(
            edges[es_mask].geometry,  # use edges that are being dropped
            poly=artifact.geometry,
            snap_to=highest_hierarchy.dissolve("coins_group").geometry,  # snap to Cs
            distance=distance,
            # buffer = highest_hierarchy.length.sum() * 1.2
        )
    split_points.extend(splitters)

    return new_connections


def loop(edges, es_mask, highest_hierarchy, artifact, distance, split_points):
    # check if we need to add a deadend to represent the space
    to_add = []
    dropped = edges[es_mask].geometry.item()
    segments = list(
        map(
            shapely.LineString,
            zip(dropped.coords[:-1], dropped.coords[1:], strict=True),
        )
    )  # TODO: vectorize this shit

    # figure out if there's a snapping node
    # Get nodes on Cs
    # bd_points = highest_hierarchy.boundary.explode()
    # Identify nodes on primes
    # TODO: what about this if/else?
    # primes = bd_points[bd_points.duplicated()]
    # # if primes.empty:
    snap_to = highest_hierarchy.dissolve("coins_group").geometry
    # else:
    #     snap_to = [primes.union_all()]

    possible_dangle, splitters = voronoi_skeleton(
        segments,  # use edges that are being dropped
        poly=artifact.geometry,
        snap_to=snap_to,
        distance=distance,
        # buffer = highest_hierarchy.length.sum() * 1.2
    )
    split_points.extend(splitters)

    possible_dangle = possible_dangle[shapely.disjoint(possible_dangle, dropped)]
    if (
        graph.Graph.build_contiguity(
            gpd.GeoSeries(possible_dangle), rook=False
        ).n_components
        == 1
    ):
        dangle_coins = momepy.COINS(
            gpd.GeoSeries(shapely.line_merge(possible_dangle)).explode(),
            flow_mode=True,
        ).stroke_gdf()
        candidate = dangle_coins.loc[dangle_coins.length.idxmax()].geometry
        if candidate.intersects(snap_to.union_all().buffer(1e-6)):
            to_add.append(candidate)

    return to_add


def split(split_points, cleaned_roads, roads):
    # split lines on new nodes
    split_points = gpd.GeoSeries(split_points)
    for split in split_points.drop_duplicates():
        _, ix = cleaned_roads.sindex.nearest(split, max_distance=1e-6)
        edge = cleaned_roads.geometry.iloc[ix]
        if edge.shape[0] == 1:
            snapped = shapely.snap(edge.item(), split, tolerance=1e-6)
            lines_split = shapely.get_parts(shapely.ops.split(snapped, split))
            cleaned_roads = pd.concat(
                [
                    cleaned_roads.drop(edge.index[0]),
                    gpd.GeoSeries(lines_split, crs=roads.crs),
                ],
                ignore_index=True,
            )

    return cleaned_roads
