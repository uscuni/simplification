import collections
import logging
import math

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
        logging.debug("CONDITION is_within False")

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
        logging.debug("CONDITION is_within False")

        new_connections, splitters = voronoi_skeleton(
            edges[es_mask].geometry,  # use edges that are being dropped
            poly=artifact.geometry,
            snap_to=highest_hierarchy.dissolve("coins_group").geometry,  # snap to Cs
            distance=distance,
            # buffer = highest_hierarchy.length.sum() * 1.2
        )
    split_points.extend(splitters)

    return new_connections


def loop(
    edges,
    es_mask,
    highest_hierarchy,
    artifact,
    distance,
    split_points,
    min_dangle_length,
):
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
        if candidate.intersects(snap_to.union_all().buffer(1e-6)) and (
            candidate.length > min_dangle_length
        ):
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


def n1_g1_identical(edges, *, to_drop, to_add, geom, distance=2, min_dangle_length=10):
    """If there is only 1 continuity group {C, E, S} and only 1 node

    - drop the edge

    Parameters
    ----------
    edges : GeoDataFrame
        geometries forming the artifact
    to_drop : list
        list collecting geometries to be dropped
    """
    to_drop.append(edges.index[0])
    dropped = edges.geometry.item()

    segments = list(
        map(
            shapely.LineString,
            zip(dropped.coords[:-1], dropped.coords[1:], strict=True),
        )
    )

    snap_to = shapely.get_point(dropped, 0)

    possible_dangle, _ = voronoi_skeleton(
        segments,  # use edges that are being dropped
        poly=geom,
        snap_to=[snap_to],
        distance=distance,
        # buffer = highest_hierarchy.length.sum() * 1.2
    )
    disjoint = shapely.disjoint(possible_dangle, dropped)
    connecting = shapely.intersects(possible_dangle, snap_to)
    dangle = possible_dangle[disjoint | connecting]

    dangle_geoms = gpd.GeoSeries(shapely.line_merge(dangle)).explode()
    dangle_coins = momepy.COINS(
        dangle_geoms, flow_mode=True, angle_threshold=120
    ).stroke_attribute()
    strokes = gpd.GeoDataFrame({"coin": dangle_coins}, geometry=dangle_geoms).dissolve(
        "coin"
    )
    entry = strokes.geometry[strokes.intersects(snap_to)].item()
    if entry.length > min_dangle_length:
        to_add.append(entry)


def nx_gx_identical(edges, *, geom, to_drop, to_add, nodes, angle):
    """If there are  1+ identical continuity groups, and more than 1 node (n>=2)

    - drop all of them and link the entry points to the centroid

    Parameters
    ----------
    edges : GeoDataFrame
        geometries forming the artifact
    geom : shapely.Polygon
        polygonal representation of the artifact
    to_drop : list
        list collecting geometries to be dropped
    to_add : list
        list collecting geometries to be added
    nodes : GeoSeries
        nodes forming the artifact
    """
    centroid = geom.centroid
    relevant_nodes = nodes.geometry.iloc[
        nodes.sindex.query(geom, predicate="intersects")
    ]

    to_drop.extend(edges.index.to_list())
    lines = shapely.shortest_line(relevant_nodes, centroid)

    # if the angle between two lines is too sharp, replace with a direct connection
    # between the nodes
    if len(lines) == 2:
        if angle_between_two_lines(lines.iloc[0], lines.iloc[1]) < angle:
            to_add.append(
                shapely.LineString([relevant_nodes.iloc[0], relevant_nodes.iloc[1]])
            )
        else:
            to_add.extend(lines.tolist())
    else:
        to_add.extend(lines.tolist())


def nx_gx(
    edges,
    *,
    artifact,
    to_drop,
    to_add,
    split_points,
    nodes,
    distance=2,
    min_dangle_length=10,
):
    """
    Drop all but highest hierarchy. If there are unconnected nodes after drop, connect
    to nearest remaining edge or nearest intersection if there are more remaining edges.
    If there three or more of the highest hierarchy, use roundabout solution.

    If, after dropping, we end up with more than one connected components based on
    remaining edges, create a connection either as a shortest line between the two or
    using skeleton if that is not inside or there are 3 or more components.

    Connection point should ideally be an existing nearest node with degree 4 or above.
    """

    # filter ends
    all_ends = edges[edges.coins_end]

    # determine if we have C present or not. Based on that, ensure that we correctly
    # pick-up the highest hierarchy and drop all lower
    if artifact.C > 0:
        # define mask for E and S strokes
        es_mask = edges.coins_group.isin(all_ends.coins_group)
        # filter Cs
        highest_hierarchy = edges[~es_mask]
    else:
        singles = set()
        visited = []
        for coins_count, group in zip(
            all_ends.coins_count, all_ends.coins_group, strict=True
        ):
            if (group not in visited) and (
                coins_count == (edges.coins_group == group).sum()
            ):
                singles.add(group)
                visited.append(group)
        # filter ends
        all_ends = edges[edges.coins_group.isin(singles)]
        # define mask for E and S strokes
        es_mask = edges.coins_group.isin(all_ends.coins_group)
        # filter Cs
        highest_hierarchy = edges[~es_mask]

    # get nodes forming the artifact
    relevant_nodes = nodes.iloc[
        nodes.sindex.query(artifact.geometry, predicate="dwithin", distance=1e-6)
    ]
    # filter nodes that lie on Cs (possibly primes)
    nodes_on_cont = relevant_nodes.index[
        relevant_nodes.sindex.query(
            highest_hierarchy.geometry.union_all(), predicate="dwithin", distance=1e-6
        )
    ]
    # get nodes that are not on Cs
    remaining_nodes = relevant_nodes.drop(nodes_on_cont)

    # get all remaining geometries and determine if they are all connected or new
    # connections need to happen
    remaining_geoms = pd.concat([remaining_nodes.geometry, highest_hierarchy.geometry])
    heads_ix, tails_ix = remaining_geoms.sindex.query(
        remaining_geoms, predicate="intersects"
    )
    n_comps = graph.Graph.from_arrays(heads_ix, tails_ix, 1).n_components

    # add list of existing edges to be removed from the network
    to_drop.extend(edges[es_mask].index.tolist())

    # more than one component in the remaining geometries
    # (either highest_hierarchy or remaining nodes)
    if n_comps > 1:
        logging.debug("CONDITION n_comps > 1 True")

        # get nodes that are relevant snapping targets (degree 4+)
        relevant_targets = relevant_nodes.loc[nodes_on_cont].query("degree > 3")

        cont_comp_labels = graph.Graph.build_contiguity(
            highest_hierarchy, rook=False
        ).component_labels
        conts_groups = highest_hierarchy.dissolve(cont_comp_labels)

        # BRANCH 1 - multiple Cs
        if len(highest_hierarchy) > 1:
            logging.debug("CONDITION len(highest_hierarchy) > 1 True")

            # Get nodes on Cs
            bd_points = highest_hierarchy.boundary.explode()
            # Identify nodes on primes
            primes = bd_points[bd_points.duplicated()]

            # For CCSS we need a special case solution if the lenght of S is
            # significantly shorter than the lenght of C. In that case, Voronoi does not
            # create shortest connections but a line that is parallel to Cs.
            if (
                highest_hierarchy.coins_group.nunique() == 2
                and artifact.S == 2
                and (highest_hierarchy.length.sum() > all_ends.length.sum())
            ):
                logging.debug("CONDITION for CCSS special case True")

                # this also appends to split_points
                new_connections = ccss_special_case(
                    primes,
                    conts_groups,
                    highest_hierarchy,
                    relevant_nodes,
                    split_points,
                )

            else:
                logging.debug("CONDITION for CCSS special case False")

                # Get new connections via skeleton
                new_connections, splitters = voronoi_skeleton(
                    edges.geometry,  # use all edges as an input
                    poly=artifact.geometry,
                    snap_to=relevant_targets.geometry,  # snap to nodes
                    distance=distance,
                    # buffer = highest_hierarchy.length.sum() * 1.2
                )
                split_points.extend(splitters)

                # The skeleton returns connections to all the nodes. We need to keep
                # only some, if there are multiple connections to a single C. We don't
                # touch the other.

                (
                    new_connections,
                    connections_intersecting_c,
                    connections_intersecting_primes,
                ) = filter_connections(primes, conts_groups, new_connections)

                # mutliple Cs that are not intersecting. Avoid forks on the ends of
                # Voronoi. If one goes to relevant node, keep it. If not, remove both
                # and replace with a new shortest connection
                if (
                    len(connections_intersecting_c) > 1
                    and len(connections_intersecting_primes) == 0
                ):
                    # this also appends to split_points
                    new_connections = avoid_forks(
                        highest_hierarchy,
                        new_connections,
                        relevant_targets,
                        artifact,
                        split_points,
                    )

                # check for disconnected Cs and reconnect
                new_connections = reconnect(
                    conts_groups, new_connections, artifact, split_points
                )

                # the drop above could've introduced a dangling edges. Remove those.
                new_connections = remove_dangles(new_connections, artifact)

        # BRANCH 2 - relevant node targets exist
        elif relevant_targets.shape[0] > 0:
            logging.debug("CONDITION relevant_targets.shape[0] > 0 True")

            # SUB BRANCH - only one remaining node
            if remaining_nodes.shape[0] < 2:
                logging.debug("CONDITION remaining_nodes.shape[0] < 2 True")

                # this also appends to split_points
                new_connections = one_remaining(
                    relevant_targets,
                    remaining_nodes,
                    artifact,
                    edges,
                    es_mask,
                    distance,
                    split_points,
                )

            # SUB BRANCH - more than one remaining node
            else:
                logging.debug("CONDITION remaining_nodes.shape[0] < 2 False")

                # this also appends to split_points
                new_connections = multiple_remaining(
                    edges,
                    es_mask,
                    artifact,
                    distance,
                    highest_hierarchy,
                    split_points,
                    relevant_targets.geometry,
                )

        # BRANCH 3 - no target nodes - snapping to C
        else:
            logging.debug(
                "CONDITION relevant_targets.shape[0] > 0 False, snapping to C"
            )

            # SUB BRANCH - only one remaining node
            if remaining_nodes.shape[0] < 2:
                logging.debug("CONDITION remaining_nodes.shape[0] < 2 True")

                # this also appends to split_points
                new_connections = one_remaining_c(
                    remaining_nodes,
                    highest_hierarchy,
                    artifact,
                    edges,
                    es_mask,
                    distance,
                    split_points,
                )

            # SUB BRANCH - more than one remaining node
            else:
                logging.debug("CONDITION remaining_nodes.shape[0] < 2 False")

                # this also appends to split_points
                new_connections = multiple_remaining(
                    edges,
                    es_mask,
                    artifact,
                    distance,
                    highest_hierarchy,
                    split_points,
                    highest_hierarchy.dissolve("coins_group").geometry,
                )

            new_connections = reconnect(
                conts_groups, new_connections, artifact, split_points
            )

        # add new connections to a list of features to be added to the network
        to_add.extend(list(new_connections))

    # there may be loops or half-loops we are dropping. If they are protruding enough
    # we want to replace them by a deadend representing their space
    elif artifact.C == 1 and (artifact.E + artifact.S) == 1:
        logging.debug("CONDITION is_loop True")

        sl = shapely.shortest_line(
            relevant_nodes.geometry.iloc[0], relevant_nodes.geometry.iloc[1]
        )
        if (
            is_within(sl, artifact.geometry)
            and sl.length < highest_hierarchy.length.sum()
        ):
            logging.debug("DEVIATION replacing with shortest")
            to_add.append(sl)
            to_drop.append(highest_hierarchy.index[0])

        dangles = loop(
            edges,
            es_mask,
            highest_hierarchy,
            artifact,
            distance,
            split_points,
            min_dangle_length,
        )
        if len(dangles) > 0:
            to_add.extend(dangles)

    elif artifact.node_count == 2 and artifact.stroke_count == 2:
        logging.debug("CONDITION is_sausage True")

        sl = shapely.shortest_line(
            relevant_nodes.geometry.iloc[0], relevant_nodes.geometry.iloc[1]
        )
        if (
            is_within(sl, artifact.geometry)
            and sl.length < highest_hierarchy.length.sum()
        ):
            logging.debug("DEVIATION replacing with shortest")
            to_add.append(sl)
            to_drop.append(highest_hierarchy.index[0])


def angle_between_two_lines(line1, line2):
    # based on momepy.coins but adapted to shapely lines
    # extract points
    a, b, c, d = shapely.get_coordinates([line1, line2]).tolist()
    a, b, c, d = tuple(a), tuple(b), tuple(c), tuple(d)

    # assertion: we expect exactly 2 of the 4 points to be identical
    # (lines touch in this point)
    points = collections.Counter([a, b, c, d])

    # points where line touch = "origin" (for vector-based angle calculation)
    origin = [k for k, v in points.items() if v == 2][0]
    # other 2 unique points (one on each line)
    point1, point2 = (k for k, v in points.items() if v == 1)

    # translate lines into vectors (numpy arrays)
    v1 = [point1[0] - origin[0], point1[1] - origin[1]]
    v2 = [point2[0] - origin[0], point2[1] - origin[1]]

    # compute angle between 2 vectors in degrees
    dot_product = v1[0] * v2[0] + v1[1] * v2[1]
    norm_v1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2)
    norm_v2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2)
    cos_theta = round(dot_product / (norm_v1 * norm_v2), 6)  # precision issues fix
    angle = math.degrees(math.acos(cos_theta))

    return angle


def simplify_singletons(artifacts, roads, nodes, distance=2):
    # collect changes
    to_drop = []
    to_add = []
    split_points = []

    planar = artifacts[~artifacts.non_planar]

    for artifact in planar.itertuples():
        # get edges relevant for an artifact
        edges = roads.iloc[roads.sindex.query(artifact.geometry, predicate="covers")]

        if (artifact.node_count == 1) and (artifact.stroke_count == 1):
            logging.debug("FUNCTION n1_g1_identical")
            n1_g1_identical(
                edges, to_drop=to_drop, to_add=to_add, geom=artifact.geometry
            )

        elif (artifact.node_count > 1) and (len(set(artifact.ces_type[1:])) == 1):
            logging.debug("FUNCTION nx_gx_identical")
            nx_gx_identical(
                edges,
                geom=artifact.geometry,
                to_add=to_add,
                to_drop=to_drop,
                nodes=nodes,
                angle=75,
            )

        elif (artifact.node_count > 1) and (len(artifact.ces_type) > 2):
            logging.debug("FUNCTION nx_gx")
            nx_gx(
                edges,
                artifact=artifact,
                to_drop=to_drop,
                to_add=to_add,
                split_points=split_points,
                nodes=nodes,
                distance=distance,
            )
        else:
            logging.debug("NON PLANAR")

    cleaned_roads = roads.geometry.drop(to_drop)
    # split lines on new nodes
    cleaned_roads = split(split_points, cleaned_roads, roads)

    # create new roads with fixed geometry. Note that to_add and to_drop lists shall be
    # global and this step should happen only once, not for every artifact
    new_roads = pd.concat(
        [
            cleaned_roads,
            gpd.GeoSeries(to_add, crs=roads.crs).line_merge().simplify(distance),
        ],
        ignore_index=True,
    )
    new_roads = momepy.remove_false_nodes(new_roads[~new_roads.is_empty])

    return new_roads
