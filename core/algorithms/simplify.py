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

from ..geometry import is_within, remove_false_nodes, snap_to_targets, voronoi_skeleton
from .common import continuity, get_stroke_info

logger = logging.getLogger(__name__)

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
    conn_c = []
    conn_p = []
    for c in conts_groups.geometry:
        int_mask = shapely.intersects(new_connections, c)
        connections_intersecting_c = new_connections[int_mask]
        conn_c.append(connections_intersecting_c)
        if len(connections_intersecting_c) > 1:
            prime_mask = shapely.intersects(
                connections_intersecting_c, primes.union_all()
            )
            connections_intersecting_primes = connections_intersecting_c[prime_mask]
            conn_p.append(connections_intersecting_primes)
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
    return (
        new_connections,
        np.concatenate(conn_c) if len(conn_c) > 0 else np.array([]),
        np.concatenate(conn_p) if len(conn_p) > 0 else np.array([]),
    )


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


def remove_dangles(new_connections, artifact, eps=1e-6):
    # the drop above could've introduced a dangling edges. Remove those.

    new_connections = shapely.line_merge(new_connections)
    pts0 = shapely.get_point(new_connections, 0)
    pts1 = shapely.get_point(new_connections, -1)
    pts = shapely.buffer(np.concatenate([pts0, pts1]), eps)
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
        logger.debug("CONDITION is_within False")

        new_connections, splitters = voronoi_skeleton(
            edges[es_mask].geometry,  # use edges that are being dropped
            poly=artifact.geometry,
            snap_to=relevant_targets.geometry.iloc[target_nearest],  # snap to nearest
            distance=distance,
            buffer=distance,  # TODO: figure out if we need this
        )
        split_points.extend(splitters)

    return remove_dangles(new_connections, artifact)


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
        logger.debug("CONDITION is_within False")

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
    eps=1e-6,
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
    bd_points = highest_hierarchy.boundary.explode()
    # Identify nodes on primes
    primes = bd_points[bd_points.duplicated()]
    if primes.empty:
        logger.debug("SNAP TO highest_hierarchy")
        snap_to = highest_hierarchy.dissolve("coins_group").geometry
    else:
        logger.debug("SNAP TO primes")
        snap_to = primes

    possible_dangle, splitters = voronoi_skeleton(
        segments,  # use edges that are being dropped
        poly=artifact.geometry,
        snap_to=snap_to,
        distance=distance,
        # buffer = highest_hierarchy.length.sum() * 1.2
    )
    split_points.extend(splitters)

    possible_dangle = possible_dangle[shapely.disjoint(possible_dangle, dropped)]
    n_comps = graph.Graph.build_contiguity(
        gpd.GeoSeries(possible_dangle), rook=False
    ).n_components
    if n_comps == 1:
        logger.debug("LOOP components 1")
        dangle_coins = momepy.COINS(
            gpd.GeoSeries(shapely.line_merge(possible_dangle)).explode(),
            flow_mode=True,
        ).stroke_gdf()
        candidate = dangle_coins.loc[dangle_coins.length.idxmax()].geometry
        if candidate.intersects(snap_to.union_all().buffer(eps)) and (
            candidate.length > min_dangle_length
        ):
            logger.debug("LOOP intersects and length > min_dangle_length")
            if not primes.empty:
                points = [
                    shapely.get_point(candidate, 0),
                    shapely.get_point(candidate, -1),
                ]
                distances = shapely.distance(points, highest_hierarchy.union_all())
                if distances.max() > min_dangle_length:
                    logger.debug("LOOP prime check passed")
                    to_add.append(candidate)
            else:
                to_add.append(candidate)
    elif n_comps > 1:
        # NOTE: it is unclear to me what exactly should happen here. I believe that
        # there will be cases when we may want to keep multiple dangles. Now keeping
        # only one.
        logger.debug("LOOP components many")
        dangle_coins = momepy.COINS(
            gpd.GeoSeries(shapely.line_merge(possible_dangle)).explode(),
            flow_mode=True,
        ).stroke_gdf()
        candidate = dangle_coins.loc[dangle_coins.length.idxmax()].geometry
        if candidate.intersects(snap_to.union_all().buffer(eps)) and (
            candidate.length > min_dangle_length
        ):
            logger.debug("LOOP intersects and length > min_dangle_length")
            if not primes.empty:
                points = [
                    shapely.get_point(candidate, 0),
                    shapely.get_point(candidate, -1),
                ]
                distances = shapely.distance(points, highest_hierarchy.union_all())
                if distances.max() > min_dangle_length:
                    logger.debug("LOOP prime check passed")
                    to_add.append(candidate)
            else:
                to_add.append(candidate)

    return to_add


def split(split_points, cleaned_roads, roads, eps=1e-6):
    # split lines on new nodes
    split_points = gpd.GeoSeries(split_points)
    for split in split_points.drop_duplicates():
        _, ix = cleaned_roads.sindex.nearest(split, max_distance=eps)
        edge = cleaned_roads.geometry.iloc[ix]
        if edge.shape[0] == 1:
            snapped = shapely.snap(edge.item(), split, tolerance=eps)
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


def nx_gx_identical(edges, *, geom, to_drop, to_add, nodes, angle, distance, eps=1e-6):
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
        nodes.sindex.query(geom, predicate="dwithin", distance=eps)
    ]

    to_drop.extend(edges.index.to_list())
    lines = shapely.shortest_line(relevant_nodes, centroid)

    if not is_within(lines, geom).all():
        logger.debug("NOT WITHIN replacing with skeleton")
        lines, _ = voronoi_skeleton(
            edges.geometry,  # use edges that are being dropped
            poly=geom,
            distance=distance,
            snap_to=relevant_nodes,
        )
        to_add.extend(lines.tolist())
    # if the angle between two lines is too sharp, replace with a direct connection
    # between the nodes
    elif len(lines) == 2:
        if angle_between_two_lines(lines.iloc[0], lines.iloc[1]) < angle:
            logger.debug(
                "TWO LINES WITH SHARP ANGLE replacing with straight connection"
            )
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
    eps=1e-6,
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
        logger.debug("HIGHEST C")
        # define mask for E and S strokes
        es_mask = edges.coins_group.isin(all_ends.coins_group)
        # filter Cs
        highest_hierarchy = edges[~es_mask]
    else:
        logger.debug("HIGHEST E")
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
        nodes.sindex.query(artifact.geometry, predicate="dwithin", distance=eps)
    ]
    # filter nodes that lie on Cs (possibly primes)
    nodes_on_cont = relevant_nodes.index[
        relevant_nodes.sindex.query(
            highest_hierarchy.geometry.union_all(), predicate="dwithin", distance=eps
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
        logger.debug("CONDITION n_comps > 1 True")

        # get nodes that are relevant snapping targets (degree 4+)
        relevant_targets = relevant_nodes.loc[nodes_on_cont].query("degree > 3")

        cont_comp_labels = graph.Graph.build_contiguity(
            highest_hierarchy, rook=False
        ).component_labels
        conts_groups = highest_hierarchy.dissolve(cont_comp_labels)

        # BRANCH 1 - multiple Cs
        if len(highest_hierarchy) > 1:
            logger.debug("CONDITION len(highest_hierarchy) > 1 True")

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
                and artifact.E == 0
                and (highest_hierarchy.length.sum() > all_ends.length.sum())
            ):
                logger.debug("CONDITION for CCSS special case True")

                # this also appends to split_points
                new_connections = ccss_special_case(
                    primes,
                    conts_groups,
                    highest_hierarchy,
                    relevant_nodes,
                    split_points,
                )

            else:
                logger.debug("CONDITION for CCSS special case False")

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
            logger.debug("CONDITION relevant_targets.shape[0] > 0 True")

            # SUB BRANCH - only one remaining node
            if remaining_nodes.shape[0] < 2:
                logger.debug("CONDITION remaining_nodes.shape[0] < 2 True")

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
                logger.debug("CONDITION remaining_nodes.shape[0] < 2 False")

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
            logger.debug("CONDITION relevant_targets.shape[0] > 0 False, snapping to C")

            # SUB BRANCH - only one remaining node
            if remaining_nodes.shape[0] < 2:
                logger.debug("CONDITION remaining_nodes.shape[0] < 2 True")

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
                logger.debug("CONDITION remaining_nodes.shape[0] < 2 False")

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
        logger.debug("CONDITION is_loop True")

        sl = shapely.shortest_line(
            relevant_nodes.geometry.iloc[0], relevant_nodes.geometry.iloc[1]
        )

        if (
            (artifact.interstitial_nodes == 0)
            and is_within(sl, artifact.geometry)
            and (sl.length * 1.1) < highest_hierarchy.length.sum()
        ):
            logger.debug("DEVIATION replacing with shortest")
            to_add.append(sl)
            to_drop.append(highest_hierarchy.index[0])

        else:
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
        logger.debug("CONDITION is_sausage True")

        sl = shapely.shortest_line(
            relevant_nodes.geometry.iloc[0], relevant_nodes.geometry.iloc[1]
        )
        if (
            is_within(sl, artifact.geometry)
            and (sl.length * 1.1) < highest_hierarchy.length.sum()
        ):
            logger.debug("DEVIATION replacing with shortest")
            to_add.append(sl)
            to_drop.append(highest_hierarchy.index[0])
    else:
        logger.debug("DROP ONLY")


def nx_gx_cluster(edges, *, cluster_geom, nodes, to_drop, to_add, distance=2, eps=1e-6):
    """treat an n-artifact cluster: merge all artifact polygons; drop
    all lines fully within the merged polygon; skeletonize and keep only
    skeletonized edges and connecting nodes"""

    # get edges on boundary
    edges_on_boundary = edges.intersection(cluster_geom.boundary.buffer(eps)).explode(
        ignore_index=True
    )
    edges_on_boundary = edges_on_boundary[
        (~edges_on_boundary.is_empty)
        & (edges_on_boundary.geom_type.str.contains("Line"))
        & (edges_on_boundary.length > 10 * eps)
    ]  # keeping only (multi)linestrings of length>>eps
    edges_on_boundary = edges_on_boundary.to_frame("geometry")

    # find nodes ON the cluster polygon boundary (to be partially kept)
    nodes_on_boundary = nodes.iloc[
        nodes.sindex.query(cluster_geom.boundary.buffer(eps), predicate="intersects")
    ].copy()

    # find edges that cross but do not lie within
    edges_crossing = edges.iloc[
        edges.sindex.query(cluster_geom.buffer(eps), predicate="crosses")
    ]

    # the nodes to keep are those that intersect with these crossing edges
    nodes_to_keep = nodes_on_boundary.iloc[
        nodes_on_boundary.sindex.query(
            edges_crossing.union_all(), predicate="intersects"
        )
    ].copy()

    # merging lines between nodes to keep:
    buffered_nodes_to_keep = nodes_to_keep.buffer(eps).union_all()

    # make queen contiguity graph on MINUSBUFFERED outline road segments,
    # and copy component labels into edges_on_boundary gdf
    edges_on_boundary = edges_on_boundary.explode(ignore_index=True)
    queen = graph.Graph.build_fuzzy_contiguity(
        edges_on_boundary.difference(buffered_nodes_to_keep)
    )
    edges_on_boundary["comp"] = queen.component_labels

    # skeletonize
    skel, _ = voronoi_skeleton(
        edges_on_boundary.dissolve(by="comp").geometry,
        cluster_geom,
        snap_to=False,
        distance=distance,
    )

    lines_to_drop = edges.iloc[
        edges.sindex.query(cluster_geom.buffer(eps), predicate="contains")
    ].index.to_list()
    lines_to_add = list(skel)

    to_add.extend(lines_to_add)
    to_drop.extend(lines_to_drop)

    ### RECONNECTING NON-PLANAR INTRUDING EDGES TO SKELETON

    # considering only edges that are kept
    edges_kept = edges.copy().drop(lines_to_drop, axis=0)

    to_reconnect = []

    skel_merged = shapely.line_merge(skel)
    skel_merged = gpd.GeoSeries(skel_merged, crs=edges.crs)

    skel_nodes = list(shapely.get_point(skel_merged, 0))
    skel_nodes.extend(list(shapely.get_point(skel_merged, -1)))
    skel_nodes = gpd.GeoSeries(skel_nodes, crs=edges.crs).union_all()

    # loop through endpoints of kept edges...
    for i in [0, -1]:
        # do the same for "end" points
        endpoints = gpd.GeoSeries(
            shapely.get_point(edges_kept.geometry, i), crs=edges.crs
        )

        # which are contained by artifact...
        endpoints = endpoints.iloc[
            endpoints.sindex.query(cluster_geom, predicate="contains")
        ]

        # ...but NOT on skeleton
        endpoints = endpoints.difference(skel_merged.union_all())

        to_reconnect.extend(endpoints.geometry)

    # to_reconnect now contains a list of points which need to be connected to the
    # nearest skel node: from those nodes, we need to add shapely shortest lines between
    # those edges_kept.endpoints and
    non_planar_connections = shapely.shortest_line(skel_nodes, to_reconnect)

    ### extend our list "to_add" with this artifact clusters' contribution:
    to_add.extend(non_planar_connections)


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


def simplify_singletons(
    artifacts, roads, distance=2, compute_coins=True, min_dangle_length=10, eps=1e-6
):
    # Get nodes from the network.
    nodes = momepy.nx_to_gdf(momepy.node_degree(momepy.gdf_to_nx(roads)), lines=False)

    # Link nodes to artifacts
    node_idx, artifact_idx = artifacts.sindex.query(
        nodes.geometry, predicate="dwithin", distance=eps
    )
    intersects = sparse.coo_array(
        ([True] * len(node_idx), (node_idx, artifact_idx)),
        shape=(len(nodes), len(artifacts)),
        dtype=np.bool_,
    )

    # Compute number of nodes per artifact
    artifacts["node_count"] = intersects.sum(axis=0)

    # Compute number of stroke groups per artifact
    if compute_coins:
        roads, _ = continuity(roads)
    strokes, c_, e_, s_ = get_stroke_info(artifacts, roads)

    artifacts["stroke_count"] = strokes
    artifacts["C"] = c_
    artifacts["E"] = e_
    artifacts["S"] = s_

    # Filter artifacts caused by non-planar intersections. (TODO: Note that this is not
    # perfect and some 3CC artifacts were non-planar but not captured here).
    artifacts["non_planar"] = artifacts["stroke_count"] > artifacts["node_count"]
    a_idx, r_idx = roads.sindex.query(artifacts.geometry.boundary, predicate="overlaps")
    artifacts.iloc[np.unique(a_idx), -1] = True

    # Count intersititial nodes (primes).
    artifacts["interstitial_nodes"] = artifacts.node_count - artifacts[
        ["C", "E", "S"]
    ].sum(axis=1)

    # Define the type label.
    ces_type = []
    for x in artifacts[["node_count", "C", "E", "S"]].itertuples():
        ces_type.append(f"{x.node_count}{'C' * x.C}{'E' * x.E}{'S' * x.S}")
    artifacts["ces_type"] = ces_type

    # collect changes
    to_drop = []
    to_add = []
    split_points = []

    planar = artifacts[~artifacts.non_planar]
    planar["buffered"] = planar.buffer(eps)
    if artifacts.non_planar.any():
        logger.debug(f"IGNORING {artifacts.non_planar.sum()} non planar artifacts")

    for artifact in planar.itertuples():
        # get edges relevant for an artifact
        edges = roads.iloc[roads.sindex.query(artifact.buffered, predicate="covers")]

        if (artifact.node_count == 1) and (artifact.stroke_count == 1):
            logger.debug("FUNCTION n1_g1_identical")
            n1_g1_identical(
                edges, to_drop=to_drop, to_add=to_add, geom=artifact.geometry
            )

        elif (artifact.node_count > 1) and (len(set(artifact.ces_type[1:])) == 1):
            logger.debug("FUNCTION nx_gx_identical")
            nx_gx_identical(
                edges,
                geom=artifact.geometry,
                to_add=to_add,
                to_drop=to_drop,
                nodes=nodes,
                angle=75,
                distance=distance,
            )

        elif (artifact.node_count > 1) and (len(artifact.ces_type) > 2):
            logger.debug("FUNCTION nx_gx")
            nx_gx(
                edges,
                artifact=artifact,
                to_drop=to_drop,
                to_add=to_add,
                split_points=split_points,
                nodes=nodes,
                distance=distance,
                min_dangle_length=min_dangle_length,
            )
        else:
            logger.debug("NON PLANAR")

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
    new_roads = remove_false_nodes(new_roads[~new_roads.is_empty].to_frame("geometry"))

    return new_roads


def simplify_clusters(artifacts, roads, distance=2, eps=1e-6):
    # Get nodes from the network.
    nodes = momepy.nx_to_gdf(momepy.node_degree(momepy.gdf_to_nx(roads)), lines=False)

    # collect changes
    to_drop = []
    to_add = []

    for _, artifact in artifacts.groupby("comp"):
        # get artifact cluster polygon
        cluster_geom = artifact.union_all(method="coverage")
        # get edges relevant for an artifact
        edges = roads.iloc[
            roads.sindex.query(cluster_geom, predicate="intersects")
        ].copy()

        nx_gx_cluster(
            edges=edges,
            cluster_geom=cluster_geom,
            nodes=nodes,
            to_drop=to_drop,
            to_add=to_add,
            eps=eps,
        )

    cleaned_roads = roads.geometry.drop(to_drop)

    # create new roads with fixed geometry. Note that to_add and to_drop lists shall be
    # global and this step should happen only once, not for every artifact
    new_roads = pd.concat(
        [
            cleaned_roads,
            gpd.GeoSeries(to_add, crs=roads.crs).line_merge().simplify(distance),
        ],
        ignore_index=True,
    ).explode()
    new_roads = remove_false_nodes(new_roads[~new_roads.is_empty].to_frame("geometry"))

    return new_roads


def consolidate_nodes(gdf, tolerance=2):
    """Return geoemtry with consolidated nodes.

    Replace clusters of nodes with a single node (weighted centroid
    of a cluster) and snap linestring geometry to it. Cluster is
    defined using DBSCAN on coordinates with ``tolerance``==``eps`.

    Does not preserve any attributes, function is purely geometric.

    Parameters
    ----------
    gdf : GeoDataFrame
        GeoDataFrame with LineStrings (usually representing street network)
    tolerance : float
        The maximum distance between two nodes for one to be considered
        as in the neighborhood of the other. Nodes within tolerance are
        considered a part of a single cluster and will be consolidated.

    Returns
    -------
    GeoSeries
    """
    # TODO: this should not dumbly merge all nodes within the cluster to a single
    # TODO: centroid but iteratively - do the two nearest and add other only if the
    # TODO: distance is still below the tolerance

    # TODO: make it work on GeoDataFrames preserving attributes
    from sklearn.cluster import DBSCAN

    nodes = momepy.nx_to_gdf((momepy.gdf_to_nx(gdf)), lines=False)

    # get clusters of nodes which should be consolidated
    db = DBSCAN(eps=tolerance, min_samples=2).fit(nodes.get_coordinates())
    nodes["lab"] = db.labels_
    change = nodes[nodes.lab > -1].set_index("lab")

    # get pygeos geometry
    geom = gdf.geometry.copy()

    # loop over clusters, cut out geometry within tolerance / 2 and replace it
    # with spider-like geometry to the weighted centroid of a cluster
    spiders = []
    midpoints = []
    for cl in change.index.unique():
        cluster = change.loc[cl]
        cookie = cluster.buffer(tolerance / 2).union_all()
        inds = geom.sindex.query(cookie, predicate="intersects")
        pts = geom.iloc[inds].intersection(cookie.boundary).get_coordinates()
        pts = shapely.get_coordinates(geom.iloc[inds].intersection(cookie.boundary))
        geom.iloc[inds] = geom.iloc[inds].difference(cookie)
        if pts.shape[0] > 0:
            midpoint = np.mean(cluster.get_coordinates(), axis=0)
            midpoints.append(midpoint)
            mids = np.array(
                [
                    midpoint,
                ]
                * len(pts)
            )
            spider = shapely.linestrings(
                np.array([pts[:, 0], mids[:, 0]]).T,
                y=np.array([pts[:, 1], mids[:, 1]]).T,
            )
            spiders.append(spider)

    # combine geometries
    geometry = pd.concat([geom, gpd.GeoSeries(np.hstack(spiders), crs=geom.crs)])

    return remove_false_nodes(geometry[~geometry.is_empty].to_frame("geometry"))